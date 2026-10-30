#!/usr/bin/env python3
"""Step 6：同樣本多 marker 交集判讀（最終回報口徑）。規則照 SKILL.md Step 6。

1. 每個 marker 取 taxonomy_decision_table.csv 的 high_quality_species_all（完整不截斷）作候選集合
2. 只有明確到 species 的二名法學名納入唯一物種判斷
3. `X sp.`、含 clone / gene / genomic 的名稱只放備註
4. 同一樣本 >= 2 組 marker -> 取候選集合交集：1 個 = 唯一物種；0 = 無共同；>1 = 仍有多個
5. 只有 1 組 marker -> 「只有一組 marker，不能做兩組交集確認」

額外（Phase 2 工單）：BLAST 結果只到屬（沒有任何二名法候選）時，查 data/gene_lookup.json，
輸出該屬「下次建議 marker」與「建議引子」。

輸出：{OUTPUT}/YYYY-MM-DD_new_orders_marker_intersection_summary.csv（UTF-8 with BOM）

用法：python3 -m pipeline.consensus_call --orders S2604280595 ... [--config ...] [--output-dir DIR]
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

from pipeline.config import gene_lookup_path, load_config, output_dir

HEADER = [
    "訂單", "樣本", "Marker 數", "是否能比對到唯一物種", "唯一物種", "判讀狀態",
    "各 marker 高品質候選", "共同候選交集", "備註", "下次建議 marker", "建議引子",
    "判讀屬", "待確認原因", "合併訂單", "判讀來源",
]


STATUS_UNIQUE = "兩組 marker 交集為唯一物種"
STATUS_NONE = "兩組 marker 無共同高品質物種"
STATUS_MULTI = "兩組 marker 交集仍有多個物種"
STATUS_SINGLE = "只有一組 marker，不能做兩組交集確認"
# 🔴-2：交集唯一，但組裝非 PASS 或有阻斷旗標 -> 黃燈，不給綠燈
STATUS_UNIQUE_PENDING = "唯一物種（待確認）"
# 🔴-4：交集只有 1 個，但「最高分／與次佳種差距」未過 -> 黃燈，不給綠燈
STATUS_INTERSECT_PENDING = "交集唯一（待覆核）"
# 🔴-6：該樣本所有 marker 都沒有可信的 BLAST 結果（與「真的比了但沒比到」區分）
STATUS_NOT_BLASTED = "尚未比對（NOT_BLASTED）"
# 🔴-3：該屬 ITS 無種級解析力，且只有 ITS -> 只報屬／section，不給種名
STATUS_GENUS_ONLY = "僅到屬／section，種名待第二 marker"
NO_CANDIDATE = "(no high-quality species candidate)"

_BINOMIAL = re.compile(r"^[A-Z][a-z]+(?:-[a-z]+)? [a-z][a-z-]*$")
_NOT_EPITHET = {"sp", "spp", "cf", "aff", "gen", "nov"}


def parse_counts(text: str) -> list[str]:
    """'A b:5; C d:2' -> ['A b', 'C d']（保留原順序）。"""
    names = []
    for part in (text or "").split(";"):
        part = part.strip()
        if not part:
            continue
        name = part.rsplit(":", 1)[0].strip() if re.search(r":\d+$", part) else part
        if name and name not in names:
            names.append(name)
    return names


def is_species_binomial(name: str, exclude_keywords: list[str]) -> bool:
    """明確到 species 的二名法學名：Genus epithet，且不是 sp. / cf. / clone / gene 等。"""
    name = name.strip()
    lowered = name.lower()
    for kw in exclude_keywords:
        if re.search(r"\b" + re.escape(kw) + r"\b", lowered):
            return False
    if not _BINOMIAL.match(name):
        return False
    epithet = name.split(" ", 1)[1]
    return epithet.rstrip(".") not in _NOT_EPITHET


def load_lookup(path: Path | None) -> dict:
    if not path or not Path(path).is_file():
        return {}
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def _primer_text(primers: dict) -> str:
    parts = []
    for key in ("fwd", "rev"):
        primer = primers.get(key)
        if primer and primer.get("name") and primer.get("seq"):
            parts.append(f"{primer['name']} {primer['seq']}")
    return " / ".join(parts)


def next_marker_suggestion(genus: str, done_markers: set[str], lookup: dict) -> tuple[str, str]:
    """回傳 (下次建議 marker, 建議引子)。優先序最高且尚未做過的 marker；同 marker 的備用引子一併列出。"""
    entry = lookup.get(genus)
    if not entry:
        return ("gene_lookup 無此屬資料，需人工查證", "")
    todo = sorted(
        (m for m in entry.get("markers", []) if m.get("marker") not in done_markers),
        key=lambda m: m.get("priority", 99),
    )
    if not todo:
        return ("gene_lookup 列出的 marker 皆已定序，需人工判斷", "")
    first_marker = todo[0]["marker"]
    chosen = [m for m in todo if m["marker"] == first_marker]
    label = "；".join(
        f"{m['marker']}（{m['gene']}）" if m.get("gene") and m["gene"] != m["marker"] else m["marker"]
        for m in chosen
    )
    # 引子欄位是 null（知識庫尚未查證）時不可輸出 "None None"，改標「引子待查，不可下單」
    primers = "；".join(
        f"{m.get('gene') or m['marker']}: {_primer_text(m['primers']) or '引子待查，不可下單'}"
        for m in chosen
        if m.get("primers")
    )
    if entry.get("its_sufficient") is True:
        label = f"ITS 即足夠；若仍要加跑：{label}"
    return (label, primers)


DEFAULT_TOP_TOLERANCE_PCT = 0.2
DEFAULT_MIN_GAP_PCT = 0.8


def parse_scores(text: str) -> dict[str, float]:
    """'A b:99.812; C d:98.1' -> {'A b': 99.812, 'C d': 98.1}"""
    out: dict[str, float] = {}
    for part in (text or "").split(";"):
        part = part.strip()
        if ":" not in part:
            continue
        name, _, val = part.rpartition(":")
        try:
            out[name.strip()] = max(out.get(name.strip(), 0.0), float(val))
        except ValueError:
            continue
    return out


def are_synonyms(a: str, b: str, synonym_groups=()) -> bool:
    """最小同義判斷：同種小名但屬名不同（例 Hypocrea lixii / Trichoderma lixii），或在 config 的同義群組內。

    taxid／現行名稱正規化已在 BLAST 層做（pipeline/taxonomy.py）；這裡只補 taxid 查不到、退回標題猜測時的舊屬名。
    """
    if a == b:
        return True
    ga, _, ea = a.partition(" ")
    gb, _, eb = b.partition(" ")
    if ea and ea == eb and ga != gb:
        return True
    return any(a in grp and b in grp for grp in synonym_groups)


def rank_check(
    by_marker: dict[str, list[dict[str, str]]],
    species: str,
    exclude_keywords: list[str],
    tolerance: float = DEFAULT_TOP_TOLERANCE_PCT,
    min_gap: float = DEFAULT_MIN_GAP_PCT,
    synonym_groups=(),
) -> list[str]:
    """🔴-4 唯一物種條件：該種在每個 marker 為最高分（差 <= tolerance），
    且至少一個 marker 與次佳非同義種差 >= min_gap（交集的鑑別力來自某個 marker）。
    回傳未通過原因；空清單＝通過。"""
    reasons: list[str] = []
    gaps: list[float] = []
    for marker in sorted(by_marker):
        scores: dict[str, float] = {}
        for r in by_marker[marker]:
            for name, val in parse_scores(r.get("species_best_pident_all", "")).items():
                if is_species_binomial(name, exclude_keywords):
                    scores[name] = max(scores.get(name, 0.0), val)
        if species not in scores:
            reasons.append(f"{marker} 無 {species} 的排名資料，無法驗證最高分與差距")
            continue
        mine = scores[species]
        others = {n: v for n, v in scores.items() if not are_synonyms(n, species, synonym_groups)}
        best_other_name = max(others, key=others.get) if others else ""
        best_other = others.get(best_other_name, 0.0)
        if best_other_name and best_other > mine + tolerance:
            reasons.append(f"{marker} 最高分是 {best_other_name} {best_other:.2f}%，{species} 只有 {mine:.2f}%")
        gaps.append(mine - best_other if best_other_name else float("inf"))
    if not reasons and gaps and max(gaps) < min_gap:
        reasons.append(f"各 marker 與次佳種差距均 < {min_gap}%（最大 {max(gaps):.2f}%）")
    return reasons


DEFAULT_BLOCKING_FLAGS = ("LOW_TOP_QCOV", "MULTIPLE_GENERA_IN_TOP20", "NO_SPECIES_PARSED", "HITLIST_SATURATED")


def type_support(by_marker: dict[str, list[dict[str, str]]]) -> tuple[list[str], str]:
    """🔴-5：type material 第二輪是否支持各 marker 的候選。

    回傳 (降信心原因, 判讀來源文字)。type 結果主判、全庫為輔：
    - TYPE_OK：候選取自 type 序列（來源 type）
    - TYPE_NO_HITS／TYPE_NOT_RUN：沿用全庫候選，只降信心（唯一物種降為「交集唯一（待覆核）」），不判錯
    - 欄位不存在（舊格式輸出）或 TYPE_DISABLED：不評論
    """
    reasons: list[str] = []
    parts: list[str] = []
    for marker in sorted(by_marker):
        states = {r.get("type_status") for r in by_marker[marker] if r.get("type_status") is not None}
        if not states:
            parts.append(f"{marker}:全庫(無 type 欄位)")
            continue
        if "TYPE_NOT_RUN" in states:
            reasons.append(f"{marker} 未做 type 第二輪")
            parts.append(f"{marker}:全庫(未做 type)")
        elif "TYPE_NO_HITS" in states:
            reasons.append(f"{marker} 無 type 序列支持")
            parts.append(f"{marker}:全庫(無 type 支持)")
        elif "TYPE_DISABLED" in states:
            parts.append(f"{marker}:全庫(type 已停用)")
        else:
            parts.append(f"{marker}:type")
    return reasons, "；".join(parts)


def full_db_disagreement(by_marker, species, exclude_keywords, tolerance, synonym_groups=()) -> list[str]:
    """type 主判時，全庫（輔助）若有更高分的另一個種，提醒（可能是誤標或 type 庫未收），不改判讀。"""
    notes = []
    for marker in sorted(by_marker):
        for r in by_marker[marker]:
            if r.get("candidate_source") != "type":
                continue
            full = {n: v for n, v in parse_scores(r.get("full_species_best_pident_all", "")).items()
                    if is_species_binomial(n, exclude_keywords)}
            mine = parse_scores(r.get("species_best_pident_all", "")).get(species)
            others = {n: v for n, v in full.items() if not are_synonyms(n, species, synonym_groups)}
            if mine is not None and others:
                top = max(others, key=others.get)
                if others[top] > mine + tolerance:
                    notes.append(f"{marker} 全庫另有更高分 {top} {others[top]:.2f}%（非 type 紀錄，可能為誤標，請人工看）")
    return notes


def gate_reasons(by_marker: dict[str, list[dict[str, str]]], blocking_flags) -> list[str]:
    """🔴-2 綠燈閘門：回傳「哪個 marker 為什麼不能亮綠燈」的清單；空清單＝全部通過。"""
    reasons: list[str] = []
    for marker in sorted(by_marker):
        for row in by_marker[marker]:
            status = (row.get("assembly_status") or "").strip()
            if status != "PASS":
                reasons.append(f"{marker} 組裝 {status or '狀態未知'}")
            flags = [f for f in (row.get("review_flag") or "").split("|") if f and f != "OK"]
            hit = [f for f in flags if f in blocking_flags]
            if hit:
                reasons.append(f"{marker} 旗標 {'|'.join(hit)}")
    return list(dict.fromkeys(reasons))


def call_sample(
    order: str,
    sample: str,
    marker_rows: list[dict[str, str]],
    exclude_keywords: list[str],
    lookup: dict,
    blocking_flags=DEFAULT_BLOCKING_FLAGS,
    unrun: dict[str, str] | None = None,
    rank_cfg: dict | None = None,
    merged_orders: list[str] | None = None,
) -> dict[str, str]:
    """對單一（訂單, 樣本）做交集判讀。marker_rows 為該樣本在 taxonomy_decision_table.csv 的列。

    unrun：這個樣本有組裝結果、卻沒進 BLAST 判讀的 marker -> 原因（只進備註，避免看起來像沒送）。
    """
    by_marker: dict[str, list[dict[str, str]]] = defaultdict(list)
    unrun = dict(unrun or {})
    for row in marker_rows:
        if "NOT_BLASTED" in (row.get("review_flag") or "").split("|"):
            unrun[row["marker"]] = "尚未送 BLAST（NOT_BLASTED）"       # 🔴-6：不算候選、不算矛盾
            continue
        by_marker[row["marker"]].append(row)
    markers = sorted(by_marker)
    if not markers:
        return {
            "訂單": order, "樣本": sample, "Marker 數": "0", "是否能比對到唯一物種": "否", "唯一物種": "",
            "判讀狀態": STATUS_NOT_BLASTED, "各 marker 高品質候選": "", "共同候選交集": "",
            "備註": "；".join(f"{m} 未進判讀：{why}" for m, why in sorted(unrun.items())),
            "下次建議 marker": "重跑 BLAST（不是重新定序）", "建議引子": "", "判讀屬": "", "待確認原因": "",
            "合併訂單": "、".join(merged_orders) if merged_orders and len(merged_orders) > 1 else "",
            "判讀來源": "",
        }

    candidates: dict[str, list[str]] = {}
    notes: list[str] = []
    genus_votes: list[str] = []
    for marker in markers:
        species: list[str] = []
        nonspecies: list[str] = []
        for row in by_marker[marker]:
            # 🔴-1：交集用完整不截斷的 high_quality_species_all；舊格式輸出沒有此欄時退回顯示欄
            full = row.get("high_quality_species_all")
            if full is None:
                full = row.get("high_quality_species_counts", "")
                if len(parse_counts(full)) >= 8:
                    notes.append(f"{marker} 候選欄為舊格式（可能被截斷在 8 種），請重跑 BLAST 解析")
            for name in parse_counts(full):
                if is_species_binomial(name, exclude_keywords):
                    if name not in species:
                        species.append(name)
                elif name not in nonspecies:
                    nonspecies.append(name)
            if row.get("candidate_genus"):
                genus_votes.append(row["candidate_genus"])
            flag = row.get("review_flag", "")
            if flag and flag != "OK":
                notes.append(f"{marker} flag: {flag}")
        candidates[marker] = species
        if nonspecies:
            notes.insert(0, f"{marker} genus-level hits: " + " | ".join(nonspecies))
            genus_votes.extend(n.split()[0] for n in nonspecies)

    per_marker_text = " || ".join(
        f"{m}: {' | '.join(candidates[m]) or NO_CANDIDATE}" for m in markers
    )
    sets = [set(candidates[m]) for m in markers]
    inter = sorted(set.intersection(*sets)) if sets else []

    unique = ""
    if len(markers) >= 2:
        if len(inter) == 1:
            status, unique = STATUS_UNIQUE, inter[0]
        elif not inter:
            status = STATUS_NONE
        else:
            status = STATUS_MULTI
    else:
        status = STATUS_SINGLE
        inter = sorted(sets[0]) if sets else []

    if merged_orders and len(merged_orders) > 1:
        notes.append("跨訂單合併判讀：" + "、".join(merged_orders))
    for marker, why in sorted((unrun or {}).items()):
        notes.append(f"{marker} 未進判讀：{why}")

    # 判讀屬：各列 candidate_genus 與候選學名的屬，多數決
    all_genus = list(genus_votes) + [n.split()[0] for m in markers for n in candidates[m]]
    genus_call = max(sorted(set(all_genus)), key=all_genus.count) if all_genus else ""
    # 🟡-1：現行屬名在知識庫沒有鍵、但舊屬名／異名有時，用 config 的 genus_aliases 對到知識庫的鍵
    aliases = (rank_cfg or {}).get("genus_aliases") or {}
    lookup_genus = genus_call if genus_call in lookup else aliases.get(genus_call, genus_call)
    entry = lookup.get(lookup_genus) or {}

    pending_reason = ""
    rank_cfg = rank_cfg or {}
    type_reasons, source_text = type_support(by_marker)
    if len(markers) >= 2:
        if status == STATUS_UNIQUE:
            reasons = gate_reasons(by_marker, blocking_flags)
            rank_reasons = rank_check(
                by_marker, unique, exclude_keywords,
                rank_cfg.get("top_tolerance_pct", DEFAULT_TOP_TOLERANCE_PCT),
                rank_cfg.get("min_gap_pct", DEFAULT_MIN_GAP_PCT),
                rank_cfg.get("synonym_groups", ()),
            ) + type_reasons
            notes.extend(full_db_disagreement(
                by_marker, unique, exclude_keywords,
                rank_cfg.get("top_tolerance_pct", DEFAULT_TOP_TOLERANCE_PCT), rank_cfg.get("synonym_groups", ()),
            ))
            if reasons:
                status = STATUS_UNIQUE_PENDING
                pending_reason = "；".join(reasons + rank_reasons)
                notes.append(f"唯一物種待確認：{pending_reason}，請先看色譜圖")
            elif rank_reasons:
                status = STATUS_INTERSECT_PENDING
                pending_reason = "；".join(rank_reasons)
                notes.append(f"交集唯一但待覆核：{pending_reason}")
    elif markers == ["ITS"] and entry.get("its_sufficient") is False:
        # 🔴-3：ITS 在這個屬分不開種，單看 ITS 不給種名
        status = STATUS_GENUS_ONLY
        pending_reason = f"{genus_call} 的 ITS 無種級解析力（gene_lookup its_sufficient=false）"
        notes.append(f"{pending_reason}；ITS 候選僅供參考：{' | '.join(inter) or NO_CANDIDATE}")
        inter = []

    # 🟡-3：只要沒有「確定的唯一物種」就給下次建議 marker（排除這一株已做過的 marker），
    # 不再只限「完全沒有種名候選」。已有唯一物種但待確認（組裝／旗標問題）時，要做的是看色譜圖，不加跑 marker。
    next_marker = next_primers = ""
    if not any(candidates.values()) and genus_votes:
        notes.append(f"BLAST 只到屬：{max(sorted(set(genus_votes)), key=genus_votes.count)}")
    if status not in (STATUS_UNIQUE, STATUS_UNIQUE_PENDING) and genus_call:
        done = set(markers) | {m for m, why in (unrun or {}).items() if "NOT_BLASTED" in why}
        next_marker, next_primers = next_marker_suggestion(lookup_genus, done, lookup)

    return {
        "訂單": order,
        "樣本": sample,
        "Marker 數": str(len(markers)),
        "是否能比對到唯一物種": (
            "待確認" if unique and status == STATUS_UNIQUE_PENDING else "待覆核" if unique and status == STATUS_INTERSECT_PENDING else ("是" if unique else "否")
        ),
        "唯一物種": unique,
        "判讀狀態": status,
        "各 marker 高品質候選": per_marker_text,
        "共同候選交集": " | ".join(inter),
        "備註": "；".join(notes),
        "下次建議 marker": next_marker,
        "建議引子": next_primers,
        "判讀屬": genus_call,
        "待確認原因": pending_reason,
        "合併訂單": "、".join(merged_orders) if merged_orders and len(merged_orders) > 1 else "",
        "判讀來源": source_text,
    }


def read_taxonomy(path: Path) -> dict[str, list[dict[str, str]]]:
    by_sample: dict[str, list[dict[str, str]]] = defaultdict(list)
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            by_sample[row["sample"]].append(row)
    return by_sample


def build_rows(out: Path, orders: list[str], cfg: dict) -> list[dict[str, str]]:
    """每個樣本一列。同一株（樣本編號相同）在不同訂單補做的 marker 會合併判讀（修正 🟡-14）：
    樣本編號永遠同一株、不重用（Rose 2026-10-30 確認），所以以樣本編號為鍵，
    連同輸出夾內其他訂單（歷史單）的同編號樣本一起交集。「訂單」欄放本次處理的最前面一單，「合併訂單」列全部來源。
    """
    exclude = cfg["consensus"]["exclude_keywords"]
    blocking = cfg["consensus"].get("blocking_flags", DEFAULT_BLOCKING_FLAGS)
    rank_cfg = cfg["consensus"]
    merge = cfg["consensus"].get("merge_across_orders", True)
    lookup = load_lookup(gene_lookup_path(cfg))

    tables: dict[str, dict[str, list[dict[str, str]]]] = {}
    for order in orders:
        path = out / f"{order}_ncbi_blast" / "taxonomy_decision_table.csv"
        if path.is_file():
            tables[order] = read_taxonomy(path)
    wanted = {s for t in tables.values() for s in t}
    history: dict[str, dict[str, list[dict[str, str]]]] = {}
    if merge and wanted:
        for path in sorted(out.glob("*_ncbi_blast/taxonomy_decision_table.csv")):
            other = path.parent.name.removesuffix("_ncbi_blast")
            if other in tables:
                continue
            table = {s: r for s, r in read_taxonomy(path).items() if s in wanted}
            if table:
                history[other] = table

    unrun_cache: dict[str, dict[str, dict[str, str]]] = {}

    def unrun_of(order: str) -> dict[str, dict[str, str]]:
        if order not in unrun_cache:
            unrun_cache[order] = unrun_markers(out, order)
        return unrun_cache[order]

    keyed: list[tuple[str, str]] = sorted(
        {(min(o for o, t in tables.items() if s in t), s) for s in wanted}
    )
    rows: list[dict[str, str]] = []
    for primary, sample in keyed:
        sources = {o: t[sample] for o, t in tables.items() if sample in t}
        sources.update({o: t[sample] for o, t in history.items() if sample in t})
        src_orders = sorted(sources)
        marker_rows = [r for o in src_orders for r in sources[o]]
        have = {r["marker"] for r in marker_rows}
        unrun: dict[str, str] = {}
        for o in src_orders:
            for m, why in unrun_of(o).get(sample, {}).items():
                if m not in have:
                    unrun[m] = why
        rows.append(
            call_sample(primary, sample, marker_rows, exclude, lookup, blocking, unrun, rank_cfg, src_orders)
        )
    return rows


def unrun_markers(out: Path, order: str) -> dict[str, dict[str, str]]:
    """{樣本: {marker: 原因}}：組裝有紀錄、但沒進 BLAST 判讀表的 marker（例如組裝 FAIL）。"""
    path = out / f"{order}_combined_assembly" / "assembly_summary.csv"
    result: dict[str, dict[str, str]] = {}
    if not path.is_file():
        return result
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for r in csv.DictReader(handle):
            result.setdefault(r["sample"], {})[r["marker"]] = f"組裝 {r.get('assembly_status', '') or '無結果'}"
    return result


def write_rows(rows: list[dict[str, str]], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=HEADER)
        writer.writeheader()
        writer.writerows(rows)
    return path


def report(rows: list[dict[str, str]]) -> str:
    unique = [r for r in rows if r["判讀狀態"] == STATUS_UNIQUE]
    notunique = [r for r in rows if r["判讀狀態"] in (STATUS_NONE, STATUS_MULTI)]
    single = [r for r in rows if r["判讀狀態"] in (STATUS_SINGLE, STATUS_GENUS_ONLY)]
    pending = [r for r in rows if r["判讀狀態"] in (STATUS_UNIQUE_PENDING, STATUS_INTERSECT_PENDING)]
    lines = [f"共 {len(rows)} 個樣本"]
    lines.append(f"能以兩組 marker 交集確認唯一物種：{len(unique)}")
    lines += [f"  - {r['訂單']} {r['樣本']}: {r['唯一物種']}" for r in unique]
    lines.append(f"交集唯一但待確認／待覆核（組裝、旗標或排名差距未過）：{len(pending)}")
    lines += [f"  - {r['訂單']} {r['樣本']}: {r['唯一物種']}（{r['待確認原因']}）" for r in pending]
    lines.append(f"兩組 marker 交集後仍不唯一：{len(notunique)}")
    lines += [f"  - {r['訂單']} {r['樣本']}: {r['判讀狀態']}" for r in notunique]
    lines.append(f"只有一組 marker、不能做雙 marker 確認：{len(single)}")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orders", nargs="+", required=True)
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None, help="結果寫到這裡；預設 = config 的輸出夾")
    parser.add_argument("--input-dir", type=Path, default=None, help="讀 BLAST 結果的夾；預設 = 輸出夾（唯讀來源可與 --output-dir 分開）")
    parser.add_argument("--date", default=None)
    args = parser.parse_args()
    cfg = load_config(args.config)
    out = args.output_dir or output_dir(cfg)
    date = args.date or datetime.now(ZoneInfo("Asia/Taipei")).strftime("%Y-%m-%d")
    rows = build_rows(args.input_dir or out, args.orders, cfg)
    path = write_rows(rows, out / f"{date}_new_orders_marker_intersection_summary.csv")
    print(report(rows))
    print(f"輸出：{path}")


if __name__ == "__main__":
    main()
