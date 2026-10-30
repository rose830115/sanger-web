#!/usr/bin/env python3
"""Run NCBI BLAST for Sanger consensus FASTA and summarize taxonomic hits.

參數來自 config.yaml 的 blast 區塊（email 預設留空，線上送件前須填實驗室信箱）。
用法：python3 -m pipeline.sanger_ncbi_blast FASTA -o OUTPUT_DIR [--assembly-summary CSV] [--reuse-xml]
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import logging
import re
import threading
import time
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

from pipeline.config import blast_settings, load_config
from pipeline.stage_guard import atomic_write_text

log = logging.getLogger("sanger.blast")

# 逾時預設值（config.yaml blast.* 可覆蓋；內部暫定，可調）
DEFAULT_REQUEST_TIMEOUT_SEC = 900
DEFAULT_RESEND_TIMEOUT_SEC = 1800
DEFAULT_MAX_RESENDS = 1
QUEUE_BUSY_MSG = "NCBI 佇列壅塞，已自動重送，請勿關閉網頁"


class BlastTimeout(Exception):
    """BLAST request 逾時且重發後仍逾時（與網路失敗區分：type 輪旗標 TIMEOUT_ABANDONED）。"""


def call_with_timeout_resend(
    call,
    timeout_sec: float,
    resend_timeout_sec: float,
    max_resends: int = DEFAULT_MAX_RESENDS,
    label: str = "BLAST",
    on_event=None,
):
    """在背景 thread 跑 call()（NCBIWWW.qblast 是黑盒輪詢，無法中途打斷），超過時限就棄掉該 thread 的結果、同參數重發。

    第一次等 timeout_sec；每次重發等 resend_timeout_sec；重發 max_resends 次後仍逾時 -> raise BlastTimeout。
    call() 自己丟的例外（網路失敗等）原樣往外丟，不重發，維持既有失敗語意。timeout_sec<=0 = 不設逾時（直接呼叫）。
    被棄掉的 thread 為 daemon：不阻擋程式結束；它可能繼續輪詢到 NCBI 回應為止，結果不採用。
    """
    if not timeout_sec or timeout_sec <= 0:
        return call()
    attempts = 1 + max(0, int(max_resends))
    for attempt in range(attempts):
        wait = timeout_sec if attempt == 0 else resend_timeout_sec
        box: dict = {}
        done = threading.Event()

        def worker(box=box, done=done):
            try:
                box["value"] = call()
            except BaseException as exc:  # noqa: BLE001
                box["error"] = exc
            finally:
                done.set()

        threading.Thread(target=worker, daemon=True, name=f"qblast-{label}-{attempt}").start()
        if done.wait(wait):
            if "error" in box:
                raise box["error"]
            return box["value"]
        if attempt + 1 < attempts:
            reason = f"request 超過 {wait:g} 秒未回應"
            log.warning("BLAST_RESENT %s：%s，放棄該 request、同參數重發（第 %d 次）", label, reason, attempt + 1)
            print(f"  BLAST_RESENT {label}：{reason}，自動重發")
            if on_event:
                on_event(QUEUE_BUSY_MSG)
    raise BlastTimeout(
        f"{label}：重發 {attempts - 1} 次後仍逾時（最後一個時窗 {wait:g} 秒），放棄 [TIMEOUT_ABANDONED]"
    )


GENUS_STOPWORDS = {
    "Uncultured",
    "Unidentified",
    "Unknown",
    "Fungal",
    "Fungus",
    "Clone",
    "Culture",
    "Isolate",
    "Strain",
    "Endophyte",
    "Ascomycota",
    "Basidiomycota",
    "Synthetic",
    "Environmental",
    "Mixed",
}
# 目、科、綱、門等高階名稱不是屬名（🟢-1：科目名被當屬名）
_HIGHER_RANK = re.compile(r"(aceae|ales|mycetes|mycetidae|mycota|mycotina|mycetes)$")
# 第二個字若是這些，代表不是「種小名」（section／complex／hybrid／construct…）
_NOT_SPECIES_WORDS = {"section", "sect", "complex", "group", "hybrid", "construct", "species", "subgen", "x"}


def parse_query_id(query_id: str) -> dict[str, str]:
    query_id = query_id.split()[0]
    parts = query_id.split("|")
    return {
        "sample": parts[0] if len(parts) > 0 else "",
        "marker": parts[1] if len(parts) > 1 else "",
        "primers": parts[2] if len(parts) > 2 else "",
    }


def normalize_query_id(query: str) -> str:
    return query.split()[0]


def first_title(hit_def: str) -> str:
    return hit_def.split(" >")[0].strip()


def organism_guess(hit_def: str) -> tuple[str, str, str]:
    title = first_title(hit_def)
    title = re.sub(r"^[A-Z]{1,3}_?\d+(?:\.\d+)?\s+", "", title)
    title = re.sub(r"\s+", " ", title).strip()

    # cf.／aff. 要排在一般二名法之前，否則會被「Genus cf」吃掉而永遠用不到（🟢-1）
    patterns = [
        r"\b([A-Z][a-z][a-z-]+)\s+(cf\.)\s+([a-z][a-z-]+)\b",
        r"\b([A-Z][a-z][a-z-]+)\s+(aff\.)\s+([a-z][a-z-]+)\b",
        r"\b([A-Z][a-z][a-z-]+)\s+(sp\.)",
        r"\b([A-Z][a-z][a-z-]+)\s+([a-z][a-z-]+)\b",
    ]
    for pattern in patterns:
        match = re.search(pattern, title)
        if not match:
            continue
        genus = match.group(1)
        if genus in GENUS_STOPWORDS or _HIGHER_RANK.search(genus):
            continue
        groups = match.groups()
        if len(groups) == 2 and groups[1] == "sp":
            groups = (groups[0], "sp.")         # "Genus sp" 與 taxonomy 的 "Genus sp." 寫法統一
        if len(groups) == 2 and groups[1] in _NOT_SPECIES_WORDS:
            return genus, genus, ""            # section／complex／hybrid／construct：只到屬，不當種名
        if len(groups) == 2:
            species = f"{groups[0]} {groups[1]}"
        else:
            species = " ".join(groups)
        return species, genus, species

    genus_match = re.search(r"\b([A-Z][a-z][a-z-]+)\b", title)
    if genus_match and genus_match.group(1) not in GENUS_STOPWORDS and not _HIGHER_RANK.search(genus_match.group(1)):
        genus = genus_match.group(1)
        return genus, genus, ""
    return "", "", ""


def marker_from_record(record_id: str, description: str) -> str:
    parsed = parse_query_id(record_id)
    if parsed["marker"]:
        return parsed["marker"]
    match = re.search(r"\|([^|]+)\|", description)
    return match.group(1) if match else ""


def batch_records(records: list, by_marker: bool) -> dict[str, list]:
    if not by_marker:
        return {"all": records}
    grouped: dict[str, list] = defaultdict(list)
    for record in records:
        grouped[marker_from_record(record.id, record.description)].append(record)
    return dict(sorted(grouped.items()))


def query_seq_hashes(records: list) -> dict[str, str]:
    """{查詢 ID: 序列 SHA1}，用來確認 XML 真的是這些序列的結果（修正 🔴-6）。"""
    return {r.id: hashlib.sha1(str(r.seq).upper().encode()).hexdigest() for r in records}


def meta_path_for(xml_path: Path) -> Path:
    return xml_path.with_name(xml_path.name + ".meta.json")


def xml_query_lengths(xml_path: Path) -> dict[str, int]:
    """讀 XML 內每個查詢的 ID 與長度（舊 XML 沒有側車檔時的最低限度比對）。"""
    out: dict[str, int] = {}
    with xml_path.open() as handle:
        for rec in NCBIXML.parse(handle):
            out[normalize_query_id(rec.query)] = int(rec.query_length)
    return out


def verified_query_ids(xml_path: Path, records: list, params: dict | None = None) -> tuple[set[str], str]:
    """回傳 (XML 內「確定是這批序列」的查詢 ID 集合, 驗證方式)。

    有側車檔 -> 逐條比序列 SHA1（exact）；沒有（舊 XML／使用者上傳）-> 比 ID＋長度（legacy）。
    XML 不存在或壞掉 -> 空集合。
    """
    if not xml_path.exists() or xml_path.stat().st_size <= 1000:
        return set(), "missing"
    want = query_seq_hashes(records)
    meta_file = meta_path_for(xml_path)
    if meta_file.is_file():
        try:
            meta = json.loads(meta_file.read_text(encoding="utf-8"))
            have = meta.get("queries", {})
        except (OSError, ValueError):
            return set(), "bad_meta"
        ok = {qid for qid, sha in want.items() if have.get(qid) == sha}
        return ok, "exact"
    try:
        lengths = xml_query_lengths(xml_path)
    except Exception:  # noqa: BLE001
        return set(), "bad_xml"
    ok = {r.id for r in records if lengths.get(r.id) == len(r.seq)}
    return ok, "legacy"


def write_xml_meta(xml_path: Path, records: list, params: dict) -> None:
    meta = {
        "queries": query_seq_hashes(records),
        "database": params.get("database"),
        "hitlist_size": params.get("hitlist_size"),
        "entrez_query": params.get("entrez_query"),
        "written": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }
    atomic_write_text(meta_path_for(xml_path), json.dumps(meta, ensure_ascii=False, indent=2))


def run_blast_batch(
    batch_name: str,
    records: list,
    outdir: Path,
    database: str,
    hitlist_size: int,
    entrez_query: str,
    email: str | None,
    tool_name: str,
    offline: bool = False,
    xml_subdir: str = "raw_blast_xml",
    request_timeout_sec: float = DEFAULT_REQUEST_TIMEOUT_SEC,
    resend_timeout_sec: float = DEFAULT_RESEND_TIMEOUT_SEC,
    max_resends: int = DEFAULT_MAX_RESENDS,
    on_event=None,
) -> tuple[Path | None, set[str]]:
    """回傳 (XML 路徑或 None, 這份 XML 內可信的查詢 ID 集合)。

    xml_subdir：XML 存放子夾。全庫第一輪用 raw_blast_xml；type material 第二輪用 raw_blast_xml_type（🔴-5）。

    XML 只在「每一條查詢序列都對得上」時沿用；對不上：線上模式舊 XML 改名 .stale 後重送，
    離線模式只採用對得上的那幾條，其餘由呼叫端標 NOT_BLASTED（絕不把舊結果掛到新序列上）。
    """
    query_fasta = "".join(record.format("fasta") for record in records)
    raw_dir = outdir / xml_subdir
    raw_dir.mkdir(parents=True, exist_ok=True)
    xml_path = raw_dir / f"{batch_name}.xml"
    params = {"database": database, "hitlist_size": hitlist_size, "entrez_query": entrez_query}
    all_ids = {r.id for r in records}

    ok_ids, how = verified_query_ids(xml_path, records, params)
    if ok_ids == all_ids:
        if how == "exact":
            try:
                meta = json.loads(meta_path_for(xml_path).read_text(encoding="utf-8"))
                same = (meta.get("database"), meta.get("hitlist_size"), meta.get("entrez_query")) == (
                    database, hitlist_size, entrez_query)
            except (OSError, ValueError):
                same = False
            if not same and not offline:
                ok_ids = set()          # 搜尋條件換了（例如 hitlist 加大、改 type 篩選）-> 線上重送
                how = "params_changed"
        if how == "legacy" and not offline:
            # 沒有側車檔的舊 XML 無從得知當時的 hitlist／篩選條件（P5b 起 hitlist 拉大、加 type 篩選）-> 線上重送
            ok_ids = set()
        if ok_ids == all_ids:
            return xml_path, ok_ids

    if offline:
        if ok_ids:
            print(f"  警告：{batch_name}.xml 只有 {len(ok_ids)}/{len(all_ids)} 條查詢對得上目前序列，其餘標 NOT_BLASTED")
            return xml_path, ok_ids
        return None, set()

    if xml_path.exists():                       # 線上：舊 XML 與目前序列不符 -> 保留為 .stale，重送
        stale = xml_path.with_name(f"{xml_path.name}.stale-{time.strftime('%Y%m%d-%H%M%S')}")
        os.replace(xml_path, stale)
        meta_file = meta_path_for(xml_path)
        if meta_file.exists():
            os.replace(meta_file, stale.with_name(stale.name + ".meta.json"))
        print(f"  {batch_name}.xml 與目前序列不符（{how}），舊檔改名 {stale.name}，重送 NCBI")

    if not email:
        raise SystemExit(
            "線上 BLAST 需要聯絡信箱：請在 config.yaml 的 blast.email 填實驗室信箱"
            "（或設環境變數 SANGER_NCBI_EMAIL）。"
        )

    NCBIWWW.tool = tool_name
    if email:
        NCBIWWW.email = email

    def _submit() -> str:
        handle = NCBIWWW.qblast(
            "blastn",
            database,
            query_fasta,
            hitlist_size=hitlist_size,
            descriptions=hitlist_size,
            alignments=hitlist_size,
            expect=10,
            megablast=True,
            entrez_query=entrez_query,
            format_type="XML",
        )
        try:
            return handle.read()
        finally:
            handle.close()

    text = call_with_timeout_resend(
        _submit, request_timeout_sec, resend_timeout_sec, max_resends,
        label=f"{xml_subdir}/{batch_name}", on_event=on_event,
    )
    atomic_write_text(xml_path, text)
    write_xml_meta(xml_path, records, params)
    return xml_path, all_ids


def best_hsp(alignment):
    return sorted(
        alignment.hsps,
        key=lambda hsp: (hsp.expect, -hsp.bits, -(hsp.identities / hsp.align_length)),
    )[0]


_SEG_ACC = re.compile(r"^(?:(?:gi\|\d+\|)?[a-z]{2,4}\|)?([A-Z]{1,3}_?\d{4,9}(?:\.\d+)?)\|?\s*(.*)$")


def hit_records(alignment) -> list[tuple[str, str]]:
    """一個 BLAST hit 可能合併多筆「序列完全相同」的資料庫紀錄（Hit_def 用 ' >' 串接）。
    回傳 [(accession(含版本), 標題), ...]，第一筆是主要紀錄。🔴-8：合併的紀錄各自可能是不同物種，不能只看第一筆。
    """
    parts = [p for p in alignment.hit_def.split(" >") if p.strip()]
    acc = str(getattr(alignment, "accession", "") or "")
    version = re.search(re.escape(acc) + r"\.(\d+)", str(getattr(alignment, "hit_id", "") or "")) if acc else None
    records = [(f"{acc}.{version.group(1)}" if version else acc, parts[0].strip() if parts else "")]
    for part in parts[1:]:
        match = _SEG_ACC.match(part.strip())
        records.append((match.group(1), match.group(2).strip()) if match else ("", part.strip()))
    return records


def hit_row(
    record, alignment, rank: int, batch_name: str, database: str, blast_date: str, source: str = "all"
) -> dict[str, str]:
    hsp = best_hsp(alignment)
    pident = hsp.identities / hsp.align_length * 100 if hsp.align_length else 0.0
    q_start = min(hsp.query_start, hsp.query_end)
    q_end = max(hsp.query_start, hsp.query_end)
    span = q_end - q_start + 1
    qcov_query = span / record.query_length * 100 if record.query_length else 0.0
    # 🔴-7：覆蓋率分母改取 query 與參考序列中較短者（參考比我們短時不再被誤排除）
    hit_length = int(getattr(alignment, "length", 0) or 0)
    denom = min(record.query_length, hit_length) if hit_length and record.query_length else record.query_length
    qcov = min(100.0, span / denom * 100) if denom else 0.0
    organism, genus, species = organism_guess(alignment.hit_def)
    records = hit_records(alignment)
    merged = []
    for _acc, title in records[1:]:
        name = organism_guess(title)[2]
        if name and name != species and name not in merged:
            merged.append(name)
    query_id = normalize_query_id(record.query)
    parsed = parse_query_id(query_id)
    gaps = getattr(hsp, "gaps", 0) or 0
    mismatches = max(0, hsp.align_length - hsp.identities - gaps)
    strand = ""
    if getattr(hsp, "frame", None):
        strand = str(hsp.frame)

    return {
        "query_id": query_id,
        "sample": parsed["sample"],
        "marker": parsed["marker"],
        "primers": parsed["primers"],
        "query_length": str(record.query_length),
        "rank": str(rank),
        "accession": alignment.accession,
        "accession_version": records[0][0],
        "hit_id": alignment.hit_id,
        "organism_guess": organism,
        "genus_guess": genus,
        "species_guess": species,
        "merged_records": str(len(records) - 1),             # 同一筆 hit 合併的其他相同序列紀錄數
        "merged_species": "; ".join(merged),                  # 這些紀錄裡「不同於主要紀錄」的物種（🔴-8，同分 hit）
        "_records": records,                                  # 內部欄位（不寫進 CSV）
        "hit_title": first_title(alignment.hit_def),
        "pident": f"{pident:.3f}",
        "qcov_pct": f"{qcov:.3f}",
        "qcov_query_pct": f"{qcov_query:.3f}",
        "hit_length": str(hit_length),
        "align_len": str(hsp.align_length),
        "identities": str(hsp.identities),
        "mismatches": str(mismatches),
        "gaps": str(gaps),
        "evalue": f"{hsp.expect:.3g}",
        "bit_score": f"{hsp.bits:.1f}",
        "query_from": str(q_start),
        "query_to": str(q_end),
        "hit_from": str(hsp.sbjct_start),
        "hit_to": str(hsp.sbjct_end),
        "hit_strand_or_frame": strand,
        "batch": batch_name,
        "database": database,
        "ref_source": source,          # all＝全庫 nt；type＝type material 第二輪（🔴-5）
        "blast_date": blast_date,
    }


def _eligible_for_taxonomy(row: dict, min_pident: float, min_qcov: float, min_align_len: int) -> bool:
    """只替「可能進候選集合」的 hit 查 taxid（pident 夠高、覆蓋與對齊長度達標，或排名前 20），省 NCBI 請求。"""
    if int(row["rank"]) <= 20:
        return True
    return hit_is_high_quality(row, min_pident, min_qcov, min_align_len)


def apply_taxonomy(
    rows: list[dict],
    resolver,
    min_pident: float,
    min_qcov: float,
    min_align_len: int,
) -> dict[str, dict]:
    """accession -> taxid -> 現行名稱，改寫每列的 species_guess／genus_guess／organism_guess（🟡-1＋🟢-1）。

    原標題猜測另存 title_species_guess；查不到 taxid 的退回標題猜測（name_source=title）。
    回傳 {accession: 解析資訊}，供寫稽核表 taxonomy_resolution.csv。
    """
    need = []
    for row in rows:
        if _eligible_for_taxonomy(row, min_pident, min_qcov, min_align_len):
            need.extend(acc for acc, _title in row["_records"] if acc)
    info = resolver.resolve(need) if (resolver is not None and need) else {}
    for row in rows:
        records = row["_records"]
        row["title_species_guess"] = row["species_guess"]
        row["taxid"] = row["taxon_name"] = row["taxon_rank"] = ""
        row["name_source"] = "title"
        names: list[str] = []
        for index, (acc, title) in enumerate(records):
            inf = info.get(acc)
            if inf and inf["species"].endswith(" sp.") and organism_guess(title)[2] not in ("", inf["species"]):
                # taxid 只到「屬 sp.」（紀錄的 organism 欄沒指定種），但標題寫了具體種名：
                # 保守起見採標題的種名當候選（寧可多一個候選讓交集更嚴，也不要讓候選憑空消失）
                inf = None
            if inf:
                species, genus = inf["species"], inf["genus"]
                organism = species or genus          # 非種級（section／complex／環境樣本）只留屬，與舊行為一致
                source = "taxonomy"
            else:
                organism, genus, species = organism_guess(title)
                source = "title"
            if index == 0:
                row["organism_guess"], row["genus_guess"], row["species_guess"] = organism, genus, species
                row["name_source"] = source
                if info.get(acc):
                    t = info[acc]
                    row["taxid"], row["taxon_name"], row["taxon_rank"] = t["taxid"], t["name"], t["rank"]
            elif species and species != row["species_guess"] and species not in names:
                names.append(species)
        row["merged_species"] = "; ".join(n for n in names if n != row["species_guess"])
    return info


def write_taxonomy_resolution(path: Path, info: dict[str, dict]) -> None:
    """稽核表：每個 accession 對到哪個 taxid、現行名稱、階層（資料庫版本會變，留當時結果）。"""
    rows = [
        {
            "accession": acc,
            "taxid": (inf or {}).get("taxid", ""),
            "current_name": (inf or {}).get("name", ""),
            "rank": (inf or {}).get("rank", ""),
            "species": (inf or {}).get("species", ""),
            "genus": (inf or {}).get("genus", ""),
            "resolved": "yes" if inf else "no",
        }
        for acc, inf in sorted(info.items())
    ]
    if rows:
        write_csv(path, rows)


def read_blast_xml(
    xml_paths: list[Path],
    database: str,
    blast_date: str,
    allowed_ids: set[str] | None = None,
    queried: set[str] | None = None,
    source: str = "all",
) -> list[dict[str, str]]:
    """allowed_ids：只收這些查詢 ID 的結果（對不上目前序列的不收）；queried：回填「XML 內有跑過的查詢」。"""
    rows: list[dict[str, str]] = []
    for xml_path in xml_paths:
        batch_name = xml_path.stem
        with xml_path.open() as handle:
            for record in NCBIXML.parse(handle):
                qid = normalize_query_id(record.query)
                if allowed_ids is not None and qid not in allowed_ids:
                    continue
                if queried is not None:
                    queried.add(qid)
                for rank, alignment in enumerate(record.alignments, start=1):
                    if not alignment.hsps:
                        continue
                    rows.append(hit_row(record, alignment, rank, batch_name, database, blast_date, source))
    return rows


def not_blasted_summary(query_id: str) -> dict[str, str]:
    """沒有可信的 BLAST 結果可用（沒送、或 XML 與目前序列對不上）：與 NO_HITS（真的送了但沒比到）區分。"""
    row = summarize_query(query_id, [], 0.0, 0.0)
    row["review_flag"] = "NOT_BLASTED"
    return row


def hit_is_high_quality(
    row: dict[str, str], min_identity: float, min_qcov: float, min_align_len: int = 0
) -> bool:
    """高品質 hit：identity、覆蓋率達標，且對齊長度 >= min_align_len（🔴-7，擋只對到保守區的短 hit）。

    query 本身短於門檻時不套對齊長度門檻（覆蓋率已要求涵蓋 95% 以上）。
    """
    if float(row["pident"]) < min_identity or float(row["qcov_pct"]) < min_qcov:
        return False
    try:
        query_length = int(row.get("query_length") or 0)
        align_len = int(row.get("align_len") or 0)
    except ValueError:
        return True
    if min_align_len and query_length >= min_align_len and align_len < min_align_len:
        return False
    return True


def species_of_row(row: dict) -> list[str]:
    """一個 hit 列涵蓋的物種：主要紀錄＋同分合併紀錄（去重、保留順序）。"""
    names = [row["species_guess"]] if row.get("species_guess") else []
    for name in (row.get("merged_species") or "").split("; "):
        if name and name not in names:
            names.append(name)
    return names


def summarize_query(
    query_id: str,
    rows: list[dict[str, str]],
    min_identity: float,
    min_qcov: float,
    min_align_len: int = 0,
    hitlist_size: int = 0,
) -> dict[str, str]:
    parsed = parse_query_id(query_id)
    if not rows:
        return {
            "query_id": query_id,
            "sample": parsed["sample"],
            "marker": parsed["marker"],
            "primers": parsed["primers"],
            "hit_count": "0",
            "top_hit_organism": "",
            "top_hit_accession": "",
            "top_pident": "",
            "top_qcov_pct": "",
            "top_evalue": "",
            "candidate_genus": "",
            "candidate_species": "",
            "top_species_counts_top20": "",
            "top_genus_counts_top20": "",
            "high_quality_species_counts": "",
            "high_quality_species_all": "",
            "species_best_pident_all": "",
            "review_flag": "NO_HITS",
        }

    sorted_rows = sorted(rows, key=lambda row: int(row["rank"]))
    top = sorted_rows[0]
    top20 = sorted_rows[:20]
    high_quality = [row for row in sorted_rows if hit_is_high_quality(row, min_identity, min_qcov, min_align_len)]
    # 前 20 名的物種／屬分布也要算進「同一 hit 併入的相同序列紀錄」，否則併入的外屬紀錄（誤標）會被旗標漏掉
    species_counts: Counter = Counter()
    genus_counts: Counter = Counter()
    for row in top20:
        names = species_of_row(row)
        species_counts.update(names)
        genera = [row["genus_guess"]] if row["genus_guess"] else []
        genera += [n.split()[0] for n in names[1:] if n.split()[0] not in genera]
        genus_counts.update(genera)
    # 🔴-8：同一 hit 合併的「序列相同」紀錄也算候選（它們與主要紀錄同分，只是 BLAST 把它們收成一筆）
    hq_species_counts: Counter = Counter()
    for row in high_quality:
        for name in species_of_row(row):
            hq_species_counts[name] += 1
    # 🔴-4：每個物種的最高 identity（覆蓋率與對齊長度達標的 hit，不限 identity 門檻），供唯一物種「最高分／差距」判斷
    best_pident: dict[str, float] = {}
    for row in sorted_rows:
        if not hit_is_high_quality(row, 0.0, min_qcov, min_align_len):
            continue
        for name in species_of_row(row):
            best_pident[name] = max(best_pident.get(name, 0.0), float(row["pident"]))

    candidate_species = species_counts.most_common(1)[0][0] if species_counts else ""
    candidate_genus = genus_counts.most_common(1)[0][0] if genus_counts else ""
    review_flags = []
    if float(top["pident"]) < min_identity:
        review_flags.append("LOW_TOP_IDENTITY")
    if float(top["qcov_pct"]) < min_qcov:
        review_flags.append("LOW_TOP_QCOV")
    if len(species_counts) > 1:
        review_flags.append("MULTIPLE_SPECIES_IN_TOP20")
    if len(genus_counts) > 1:
        review_flags.append("MULTIPLE_GENERA_IN_TOP20")
    if not high_quality:
        review_flags.append("NO_HIGH_QUALITY_HITS")
    # 🔴-8：清單被塞滿且最後一筆仍達高品質 -> 候選集合可能在清單邊界被截斷，不可當完整集合
    if hitlist_size and len(sorted_rows) >= hitlist_size and hit_is_high_quality(
        sorted_rows[-1], min_identity, min_qcov, min_align_len
    ):
        review_flags.append("HITLIST_SATURATED")
    if not candidate_species:
        review_flags.append("NO_SPECIES_PARSED")

    return {
        "query_id": query_id,
        "sample": parsed["sample"],
        "marker": parsed["marker"],
        "primers": parsed["primers"],
        "hit_count": str(len(sorted_rows)),
        "top_hit_organism": top["organism_guess"],
        "top_hit_accession": top["accession"],
        "top_pident": top["pident"],
        "top_qcov_pct": top["qcov_pct"],
        "top_evalue": top["evalue"],
        "candidate_genus": candidate_genus,
        "candidate_species": candidate_species,
        "top_species_counts_top20": "; ".join(f"{name}:{count}" for name, count in species_counts.most_common(8)),
        "top_genus_counts_top20": "; ".join(f"{name}:{count}" for name, count in genus_counts.most_common(8)),
        # 前 8 名只給人看；交集判讀一律用下面完整不截斷的欄位（修正 🔴-1）
        "high_quality_species_counts": "; ".join(f"{name}:{count}" for name, count in hq_species_counts.most_common(8)),
        "high_quality_species_all": "; ".join(f"{name}:{count}" for name, count in hq_species_counts.most_common()),
        "species_best_pident_all": "; ".join(
            f"{name}:{score:.3f}" for name, score in sorted(best_pident.items(), key=lambda kv: -kv[1])
        ),
        "review_flag": "OK" if not review_flags else "|".join(review_flags),
    }


TYPE_OK, TYPE_NO_HITS, TYPE_NOT_RUN, TYPE_DISABLED = "TYPE_OK", "TYPE_NO_HITS", "TYPE_NOT_RUN", "TYPE_DISABLED"
# type 結果主判時，被取代成 type 版本（原全庫版另存 full_ 前綴欄位）的摘要欄位
PRIMARY_FIELDS = (
    "hit_count", "top_hit_organism", "top_hit_accession", "top_pident", "top_qcov_pct", "top_evalue",
    "candidate_genus", "candidate_species", "top_species_counts_top20", "top_genus_counts_top20",
    "high_quality_species_counts", "high_quality_species_all", "species_best_pident_all", "review_flag",
)
FULL_KEEP = (
    "hit_count", "top_hit_organism", "top_hit_accession", "top_pident", "candidate_species",
    "high_quality_species_all", "species_best_pident_all", "review_flag",
)


def combine_with_type(full: dict[str, str], tsum: dict[str, str] | None, state: str) -> dict[str, str]:
    """把全庫第一輪摘要與 type 第二輪摘要合成一列（🔴-5）。

    type 有高品質候選（TYPE_OK）-> 主判欄位換成 type 版本，全庫版另存 full_*，candidate_source=type。
    第二輪沒有高品質 hit（TYPE_NO_HITS）／沒跑（TYPE_NOT_RUN）-> 沿用全庫版，candidate_source=全庫，
    只降信心（review_flag 加 NO_TYPE_SUPPORT／TYPE_NOT_RUN），不判錯。TYPE_DISABLED＝設定關閉第二輪。
    """
    out = dict(full)
    for key in FULL_KEEP:
        out[f"full_{key}"] = full.get(key, "")
    out["type_status"] = state
    out["type_hit_count"] = tsum["hit_count"] if tsum else ""
    out["type_top_hit_organism"] = tsum["top_hit_organism"] if tsum else ""
    out["type_top_pident"] = tsum["top_pident"] if tsum else ""
    out["candidate_source"] = "全庫"
    if state == TYPE_OK and tsum:
        for key in PRIMARY_FIELDS:
            out[key] = tsum[key]
        out["candidate_source"] = "type"
        return out
    if full["review_flag"] == "NOT_BLASTED":
        out["candidate_source"] = ""
        return out
    extra = {TYPE_NO_HITS: "NO_TYPE_SUPPORT", TYPE_NOT_RUN: "TYPE_NOT_RUN"}.get(state)
    if extra:
        flags = [] if full["review_flag"] == "OK" else full["review_flag"].split("|")
        out["review_flag"] = "|".join(flags + [extra])
    return out


def type_state(qid: str, tsum: dict[str, str] | None, type_queried: set[str], type_on: bool) -> str:
    if not type_on:
        return TYPE_DISABLED
    if qid not in type_queried or tsum is None:
        return TYPE_NOT_RUN
    return TYPE_OK if tsum.get("high_quality_species_all") else TYPE_NO_HITS


def summarize_samples(query_summaries: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in query_summaries:
        grouped[row["sample"]].append(row)

    sample_rows = []
    for sample, rows in sorted(grouped.items()):
        genus_counts = Counter(row["candidate_genus"] for row in rows if row["candidate_genus"])
        species_counts = Counter(row["candidate_species"] for row in rows if row["candidate_species"])
        markers = ";".join(f"{row['marker']}:{row['candidate_species'] or row['candidate_genus'] or 'NA'}" for row in rows)
        flags = sorted({flag for row in rows for flag in row["review_flag"].split("|") if flag and flag != "OK"})
        marker_count = len(rows)
        genus_count = len(genus_counts)
        species_count = len(species_counts)
        if all(row["review_flag"] == "NOT_BLASTED" for row in rows):
            overall = "NOT_BLASTED"
        elif not genus_counts:
            overall = "NO_RELIABLE_MATCH"
        elif marker_count > 1 and genus_count > 1:
            overall = "CONFLICT_BETWEEN_MARKERS"
        elif species_count > 1:
            overall = "GENUS_OR_SPECIES_REVIEW"
        elif flags:
            overall = "REVIEW"
        else:
            overall = "CONSISTENT"

        sample_rows.append(
            {
                "sample": sample,
                "marker_count": str(marker_count),
                "candidate_genus": genus_counts.most_common(1)[0][0] if genus_counts else "",
                "candidate_species": species_counts.most_common(1)[0][0] if species_counts else "",
                "marker_calls": markers,
                "all_candidate_genera": "; ".join(f"{name}:{count}" for name, count in genus_counts.most_common()),
                "all_candidate_species": "; ".join(f"{name}:{count}" for name, count in species_counts.most_common()),
                "overall_review_status": overall,
                "review_flags": "|".join(flags) if flags else "OK",
            }
        )
    return sample_rows


def read_assembly_summary(path: Path | None) -> dict[tuple[str, str], dict[str, str]]:
    if not path or not path.exists():
        return {}
    with path.open(newline="", encoding="utf-8-sig") as handle:      # 容許 Excel 存過的 BOM
        return {
            (row["sample"], row["marker"]): row
            for row in csv.DictReader(handle)
        }


def suggested_action(row: dict[str, str], assembly_row: dict[str, str] | None) -> str:
    flags = row["review_flag"].split("|") if row["review_flag"] else []
    assembly_status = assembly_row.get("assembly_status", "") if assembly_row else ""
    if "NOT_BLASTED" in flags:
        return "NOT_BLASTED_RERUN_BLAST"
    if assembly_status and assembly_status != "PASS":
        return "CHECK_ASSEMBLY_OR_CHROMATOGRAM_FIRST"
    if "LOW_TOP_IDENTITY" in flags or "LOW_TOP_QCOV" in flags or "NO_HIGH_QUALITY_HITS" in flags:
        return "CHECK_TOP_HITS_AND_CONSIDER_ADDITIONAL_MARKER"

    if "MULTIPLE_GENERA_IN_TOP20" in flags:
        return "CHECK_TAXONOMY_SYNONYMS_OR_ADD_MARKER"
    if "MULTIPLE_SPECIES_IN_TOP20" in flags:
        return "SPECIES_NEEDS_MANUAL_TOP_HIT_REVIEW"
    if "HITLIST_SATURATED" in flags:
        return "HITLIST_SATURATED_CHECK_TYPE_SEQUENCES"
    if "NO_TYPE_SUPPORT" in flags or "TYPE_NOT_RUN" in flags:
        return "TYPE_SUPPORT_MISSING_REVIEW"       # 🔴-5：全庫有候選但沒有 type 序列支持，不給「可信」
    if row["review_flag"] == "OK":
        return "CANDIDATE_OK_IF_SEQUENCE_QC_ACCEPTED"
    return "MANUAL_REVIEW"


def taxonomy_decision_rows(
    query_summaries: list[dict[str, str]],
    assembly_by_key: dict[tuple[str, str], dict[str, str]],
) -> list[dict[str, str]]:
    rows = []
    for row in query_summaries:
        assembly = assembly_by_key.get((row["sample"], row["marker"]), {})
        rows.append(
            {
                "sample": row["sample"],
                "marker": row["marker"],
                "primers": row["primers"],
                "assembly_status": assembly.get("assembly_status", ""),
                "assembly_identity_pct": assembly.get("identity_pct", ""),
                "assembly_low_confidence_sites": assembly.get("low_confidence_sites", ""),
                "top_hit_organism": row["top_hit_organism"],
                "top_hit_accession": row["top_hit_accession"],
                "top_pident": row["top_pident"],
                "top_qcov_pct": row["top_qcov_pct"],
                "candidate_genus": row["candidate_genus"],
                "candidate_species": row["candidate_species"],
                "top_species_counts_top20": row["top_species_counts_top20"],
                "high_quality_species_counts": row["high_quality_species_counts"],
                "high_quality_species_all": row["high_quality_species_all"],
                "species_best_pident_all": row["species_best_pident_all"],
                "review_flag": row["review_flag"],
                "candidate_source": row.get("candidate_source", ""),
                "type_status": row.get("type_status", ""),
                "type_top_hit_organism": row.get("type_top_hit_organism", ""),
                "type_top_pident": row.get("type_top_pident", ""),
                "full_candidate_species": row.get("full_candidate_species", ""),
                "full_top_hit_organism": row.get("full_top_hit_organism", ""),
                "full_top_pident": row.get("full_top_pident", ""),
                "full_species_best_pident_all": row.get("full_species_best_pident_all", ""),
                "suggested_action": suggested_action(row, assembly),
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    """先寫暫存檔、完整寫完才改名，避免中途失敗留下半份表（修正 🟡-13）。"""
    import io  # noqa: PLC0415

    buf = io.StringIO(newline="")
    if rows:
        writer = csv.DictWriter(buf, fieldnames=[k for k in rows[0] if not k.startswith("_")], extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    atomic_write_text(path, buf.getvalue())


def make_resolver(blast_cfg: dict, offline_net: bool | None = None):
    """依設定建立 TaxonomyResolver。沒信箱（或 taxonomy_lookup=false）-> 只吃快取、不連線。"""
    from pipeline.taxonomy import TaxonomyResolver  # noqa: PLC0415

    email = blast_cfg.get("email") or None
    net_off = (not email) or (not blast_cfg.get("taxonomy_lookup", True)) if offline_net is None else offline_net
    return TaxonomyResolver(
        cache_path=Path(blast_cfg["taxonomy_cache_file"]) if blast_cfg.get("taxonomy_cache_file") else None,
        email=email,
        tool=blast_cfg.get("tool_name", "sanger-web-blast"),
        offline=net_off,
        max_age_days=int(blast_cfg.get("taxonomy_cache_days", 90)),
    )


def run_blast(
    fasta: Path,
    output_dir: Path,
    blast_cfg: dict,
    assembly_summary: Path | None = None,
    by_marker: bool = True,
    offline: bool = False,
    resolver=None,
    on_event=None,
) -> Path:
    """跑（或離線重解析）一個訂單的 BLAST 並寫出全部結果表。

    offline=True：只讀 raw_blast_xml 既有 XML，缺檔就報錯，絕不送 BLAST。
    學名正規化（accession->taxid->現行名稱）另由 blast.taxonomy_lookup 與信箱決定是否連 E-utilities；
    沒信箱或 taxonomy_lookup=false 時只吃本機快取，查不到的退回標題猜測。resolver 可由呼叫端注入（測試用）。
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    records = list(SeqIO.parse(str(fasta), "fasta"))
    if not records:
        raise SystemExit(f"No FASTA records found in {fasta}")

    database = blast_cfg["database"]
    min_identity = blast_cfg["min_identity"]
    min_qcov = blast_cfg["min_qcov"]
    min_align_len = int(blast_cfg.get("min_align_len", 0))
    hitlist = int(blast_cfg["hitlist_size"])
    type_on = bool(blast_cfg.get("type_round", True))
    type_query = blast_cfg.get("type_entrez_query", "sequence from type[filter]")
    blast_date = datetime.now().isoformat(timespec="seconds")
    xml_paths: list[Path] = []
    type_xml_paths: list[Path] = []
    trusted_ids: set[str] = set()
    type_trusted: set[str] = set()
    type_errors: list[str] = []
    type_timeouts: list[str] = []
    timeouts = dict(
        request_timeout_sec=float(blast_cfg.get("request_timeout_sec", DEFAULT_REQUEST_TIMEOUT_SEC)),
        resend_timeout_sec=float(blast_cfg.get("resend_timeout_sec", DEFAULT_RESEND_TIMEOUT_SEC)),
        max_resends=int(blast_cfg.get("max_resends", DEFAULT_MAX_RESENDS)),
        on_event=on_event,
    )
    batches = batch_records(records, by_marker=by_marker)
    for index, (batch_name, batch_records_list) in enumerate(batches.items(), start=1):
        safe_batch = re.sub(r"[^A-Za-z0-9_.-]+", "_", batch_name or f"batch_{index}")
        print(f"Submitting/parsing BLAST batch {safe_batch}: {len(batch_records_list)} sequence(s)")
        args = (
            safe_batch, batch_records_list, output_dir, database, blast_cfg["hitlist_size"],
        )
        tail = (blast_cfg.get("email") or None, blast_cfg["tool_name"])
        xml_path, ok_ids = run_blast_batch(
            args[0], args[1], args[2], args[3], args[4], blast_cfg["entrez_query"], *tail, offline=offline,
            **timeouts,
        )
        if xml_path is not None:
            xml_paths.append(xml_path)
            trusted_ids |= ok_ids
        if type_on:
            # 🔴-5：type material 第二輪。失敗只降信心（該 marker 標 TYPE_NOT_RUN），不讓整單中斷
            if not offline:
                time.sleep(blast_cfg["delay_between_markers_sec"])
            try:
                txml, tok = run_blast_batch(
                    args[0], args[1], args[2], args[3], args[4], type_query, *tail,
                    offline=offline, xml_subdir="raw_blast_xml_type", **timeouts,
                )
            except SystemExit:
                raise
            except BlastTimeout as exc:
                txml, tok = None, set()
                type_errors.append(f"{safe_batch}: TIMEOUT_ABANDONED {exc}")
                type_timeouts.append(safe_batch)
                print(f"  警告：{safe_batch} type 第二輪逾時放棄（{exc}），該 marker 只降信心、沿用全庫結果")
            except Exception as exc:  # noqa: BLE001
                txml, tok = None, set()
                type_errors.append(f"{safe_batch}: {exc}")
                print(f"  警告：{safe_batch} type 第二輪失敗（{exc}），該 marker 只降信心、沿用全庫結果")
            if txml is not None:
                type_xml_paths.append(txml)
                type_trusted |= tok
        if index < len(batches) and not offline:
            time.sleep(blast_cfg["delay_between_markers_sec"])

    if offline and not trusted_ids:
        raise SystemExit(f"離線模式：找不到與目前序列對得上的既有 BLAST XML（{output_dir / 'raw_blast_xml'}）")
    queried: set[str] = set()
    hit_rows = read_blast_xml(xml_paths, database, blast_date, allowed_ids=trusted_ids, queried=queried)
    hit_rows = sorted(hit_rows, key=lambda row: (row["sample"], row["marker"], int(row["rank"])))
    type_queried: set[str] = set()
    type_hit_rows: list[dict[str, str]] = []
    if type_on:
        type_hit_rows = read_blast_xml(
            type_xml_paths, database, blast_date, allowed_ids=type_trusted, queried=type_queried, source="type"
        )
        type_hit_rows = sorted(type_hit_rows, key=lambda row: (row["sample"], row["marker"], int(row["rank"])))
    # 🟡-1＋🟢-1：學名改走 accession -> taxid -> 現行名稱（全庫與 type 兩批 hit 一次查）
    if resolver is None and blast_cfg.get("taxonomy_lookup", True):
        resolver = make_resolver(blast_cfg)
    tax_info: dict[str, dict] = {}
    if resolver is not None:
        tax_info = apply_taxonomy(hit_rows + type_hit_rows, resolver, float(blast_cfg.get("taxonomy_min_pident", 95.0)),
                                  min_qcov, min_align_len)
        write_taxonomy_resolution(output_dir / "taxonomy_resolution.csv", tax_info)
        for err in resolver.errors[:3]:
            print(f"  警告：taxonomy 查詢 {err}（該批退回標題猜測）")
    else:
        for row in hit_rows + type_hit_rows:
            row.update(title_species_guess=row["species_guess"], taxid="", taxon_name="", taxon_rank="", name_source="title")
    write_csv(output_dir / "blast_top100_hits.csv", hit_rows)
    if type_on:
        write_csv(output_dir / "blast_type_hits.csv", type_hit_rows)

    grouped_hits: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in hit_rows:
        grouped_hits[row["query_id"]].append(row)
    grouped_type: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in type_hit_rows:
        grouped_type[row["query_id"]].append(row)

    query_summaries = []
    for record in records:
        if record.id in queried:
            full = summarize_query(record.id, grouped_hits.get(record.id, []), min_identity, min_qcov, min_align_len, hitlist)
        else:
            full = not_blasted_summary(record.id)
        tsum = (
            summarize_query(record.id, grouped_type.get(record.id, []), min_identity, min_qcov, min_align_len, hitlist)
            if record.id in type_queried else None
        )
        state = type_state(record.id, tsum, type_queried, type_on)
        query_summaries.append(combine_with_type(full, tsum, state))
    query_summaries = sorted(query_summaries, key=lambda row: (row["sample"], row["marker"]))
    write_csv(output_dir / "blast_query_summary.csv", query_summaries)
    write_csv(output_dir / "species_decision_table.csv", summarize_samples(query_summaries))
    assembly_by_key = read_assembly_summary(assembly_summary)
    write_csv(
        output_dir / "taxonomy_decision_table.csv",
        taxonomy_decision_rows(query_summaries, assembly_by_key),
    )
    write_csv(
        output_dir / "blast_needs_review.csv",
        [row for row in query_summaries if row["review_flag"] != "OK"],
    )

    metadata = {
        "fasta": str(fasta),
        "database": database,
        "hitlist_size": str(blast_cfg["hitlist_size"]),
        "entrez_query": blast_cfg["entrez_query"],
        "type_round": str(type_on),
        "type_entrez_query": type_query if type_on else "",
        "type_round_errors": "; ".join(type_errors),
        "type_round_timeouts": "; ".join(type_timeouts),
        "taxonomy_accessions": str(len(tax_info)),
        "taxonomy_unresolved": str(sum(1 for v in tax_info.values() if not v)),
        "taxonomy_errors": "; ".join(resolver.errors[:3]) if resolver is not None else "",
        "min_identity": str(min_identity),
        "min_qcov": str(min_qcov),
        "min_align_len": str(min_align_len),
        "qcov_basis": "min(query_length, hit_length)",
        "blast_date": blast_date,
        "xml_files": ";".join(str(path) for path in xml_paths),
        "type_xml_files": ";".join(str(path) for path in type_xml_paths),
    }
    with (output_dir / "blast_run_metadata.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(metadata.keys()))
        writer.writeheader()
        writer.writerow(metadata)

    return output_dir


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta", type=Path)
    parser.add_argument("-o", "--output-dir", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=None)
    # 以下參數預設取自 config.yaml 的 blast 區塊；命令列給了就覆蓋
    parser.add_argument("--database")
    parser.add_argument("--hitlist-size", type=int)
    parser.add_argument("--entrez-query")
    parser.add_argument("--email")
    parser.add_argument("--tool-name")
    parser.add_argument("--min-identity", type=float)
    parser.add_argument("--min-qcov", type=float)
    parser.add_argument("--assembly-summary", type=Path)
    parser.add_argument("--no-marker-batches", action="store_true")
    parser.add_argument(
        "--reuse-xml",
        action="store_true",
        help="離線模式：只重解析 raw_blast_xml 內既有 XML，不連線 NCBI。",
    )
    args = parser.parse_args()

    blast_cfg = blast_settings(load_config(args.config))
    for key in ("database", "hitlist_size", "entrez_query", "email", "tool_name", "min_identity", "min_qcov"):
        value = getattr(args, key)
        if value is not None:
            blast_cfg[key] = value
    out = run_blast(
        args.fasta,
        args.output_dir,
        blast_cfg,
        assembly_summary=args.assembly_summary,
        by_marker=not args.no_marker_batches,
        offline=args.reuse_xml,
    )
    print(out)


if __name__ == "__main__":
    main()
