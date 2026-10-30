"""輸入雜湊與階段失效機制（修正 🔴-6）。

問題：續跑與快取原本只看「資料夾／檔案存不存在」，改名、重傳、分批到貨後舊結果會掛到新序列上。
做法：每個訂單存一份「輸入清單」{order}_input_manifest.json，記錄
  - 各子批次 ab1 的 SHA1（QC 依這個）
  - 該訂單全部 ab1 的 SHA1 + qc/assembly 設定 + gene_lookup 雜湊（Assembly 依這個）
  - Assembly 共識 FASTA 的 SHA1 + BLAST 相關設定（BLAST 解析表依這個）
上游指紋一變，下游就視為失效：舊輸出夾改名搬進 _stale/（不刪除），重新產生。
BLAST XML 另有逐條查詢序列雜湊（見 sanger_ncbi_blast.py），不符就重送。
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import tempfile
import time
from contextlib import contextmanager
from pathlib import Path

MANIFEST_VERSION = 1
# 影響 BLAST 解析表的設定（改了就要重新解析；XML 是否重送由 XML 側車檔另外判斷）
BLAST_PARSE_KEYS = (
    "database", "hitlist_size", "entrez_query", "min_identity", "min_qcov", "min_align_len",
    "type_round", "type_entrez_query",
)
# 解析邏輯版本：解析規則變動（例如 qcov 分母改法）時加 1，讓舊解析表失效
PARSE_LOGIC_VERSION = 6        # 3＝P5b：type 第二輪、合併 defline 物種、taxid 正規化

_file_hash_cache: dict[tuple[str, int, int], str] = {}


def sha1_file(path: Path) -> str:
    st = path.stat()
    key = (str(path), st.st_size, st.st_mtime_ns)
    if key not in _file_hash_cache:
        h = hashlib.sha1()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1 << 20), b""):
                h.update(chunk)
        _file_hash_cache[key] = h.hexdigest()
    return _file_hash_cache[key]


def sha1_text(text: str) -> str:
    return hashlib.sha1(text.encode("utf-8")).hexdigest()


def sha1_json(obj) -> str:
    return sha1_text(json.dumps(obj, sort_keys=True, ensure_ascii=False, default=str))


def ab1_hashes(folder: Path) -> dict[str, str]:
    """{檔名: sha1}；收檔規則與 QC／Assembly 的 glob("*.ab1") 一致。"""
    return {p.name: sha1_file(p) for p in sorted(folder.glob("*.ab1")) if p.is_file()}


def lookup_hashes() -> dict[str, str]:
    """引子→marker 映射實際讀的 gene_lookup（markers.py 寫死讀的那份）。"""
    from pipeline.markers import LOOKUP_PATH  # noqa: PLC0415

    return {"markers_lookup": sha1_file(LOOKUP_PATH) if LOOKUP_PATH.is_file() else ""}


def qc_fingerprint(files: dict[str, str], cfg: dict) -> str:
    return sha1_json({"v": MANIFEST_VERSION, "files": files, "qc": cfg["qc"], **lookup_hashes()})


def assembly_fingerprint(files_by_sub: dict[str, dict[str, str]], cfg: dict) -> str:
    return sha1_json(
        {"v": MANIFEST_VERSION, "files": files_by_sub, "qc": cfg["qc"], "assembly": cfg["assembly"], **lookup_hashes()}
    )


def blast_fingerprint(assembly_dir: Path, cfg: dict) -> str | None:
    fasta = assembly_dir / "paired_consensus_for_blast.fasta"
    if not fasta.is_file():
        return None
    parse_cfg = {k: cfg["blast"].get(k) for k in BLAST_PARSE_KEYS}
    return sha1_json({"fasta": sha1_file(fasta), "blast": parse_cfg, "logic": PARSE_LOGIC_VERSION})


# ---------------------------------------------------------------- 輸入清單檔
def manifest_path(out: Path, order: str) -> Path:
    return out / f"{order}_input_manifest.json"


def load_manifest(out: Path, order: str) -> dict:
    path = manifest_path(out, order)
    if path.is_file():
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            return {}
    return {}


def save_manifest(out: Path, order: str, manifest: dict) -> None:
    manifest["version"] = MANIFEST_VERSION
    path = manifest_path(out, order)
    path.parent.mkdir(parents=True, exist_ok=True)
    atomic_write_text(path, json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True))


def atomic_write_text(path: Path, text: str) -> None:
    fd, tmp = tempfile.mkstemp(dir=str(path.parent), prefix=path.name + ".", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            handle.write(text)
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


# ---------------------------------------------------------------- 失效判定
def current_inputs(raw: Path, order_info: dict, cfg: dict) -> dict:
    """算出一個訂單目前的輸入指紋：{"files": {sub: {name: sha}}, "qc": {sub: fp}, "assembly": fp}。"""
    files = {sub: ab1_hashes(raw / sub) for sub in order_info["subfolders"] if (raw / sub).is_dir()}
    return {
        "files": files,
        "qc": {sub: qc_fingerprint(f, cfg) for sub, f in files.items()},
        "assembly": assembly_fingerprint(files, cfg),
    }


def judge(output_exists: bool, recorded: str | None, current: str | None) -> str:
    """回傳 'missing'（要跑）/ 'fresh'（可沿用）/ 'stale'（上游變了，要重做）/ 'legacy'（有輸出但沒輸入紀錄）。"""
    if not output_exists:
        return "missing"
    if recorded is None:
        return "legacy"
    return "fresh" if recorded == current else "stale"


def stage_report(raw: Path, out: Path, order: str, order_info: dict, cfg: dict) -> dict:
    """唯讀：各階段狀態。scan 與 run_all 共用。"""
    cur = current_inputs(raw, order_info, cfg)
    man = load_manifest(out, order)
    qc = {
        sub: judge((out / f"{sub}_processed").is_dir(), (man.get("qc") or {}).get(sub, {}).get("fp"), cur["qc"][sub])
        for sub in cur["qc"]
    }
    asm_dir = out / f"{order}_combined_assembly"
    asm = judge(asm_dir.is_dir(), (man.get("assembly") or {}).get("fp"), cur["assembly"])
    if asm_dir.is_dir() and any(st == "missing" for st in qc.values()):
        asm = "stale"   # 有子批次還沒做過 QC（例如晚到的第二批）：既有 Assembly 可能漏了它，重做（離線、不連 NCBI）
    blast_dir = out / f"{order}_ncbi_blast"
    cur_blast = blast_fingerprint(asm_dir, cfg)
    blast = judge((blast_dir / "taxonomy_decision_table.csv").is_file(), (man.get("blast") or {}).get("fp"), cur_blast)
    if asm in ("missing", "stale"):
        blast = "missing" if blast == "missing" else "stale"
    return {"current": cur, "manifest": man, "qc": qc, "assembly": asm, "blast": blast, "blast_fp": cur_blast}


# ---------------------------------------------------------------- 輸出夾：暫存→成功才改名
def retire_to_stale(out: Path, target: Path) -> Path | None:
    """舊輸出夾搬進 out/_stale/（不刪除）。"""
    if not target.exists():
        return None
    stale_root = out / "_stale"
    stale_root.mkdir(parents=True, exist_ok=True)
    dest = stale_root / f"{target.name}_{time.strftime('%Y%m%d-%H%M%S')}"
    n = 1
    while dest.exists():
        dest = stale_root / f"{target.name}_{time.strftime('%Y%m%d-%H%M%S')}_{n}"
        n += 1
    shutil.move(str(target), str(dest))
    return dest


@contextmanager
def atomic_dir(target: Path):
    """with atomic_dir(target) as tmp: 寫進 tmp；沒丟例外才改名成 target，失敗就清掉暫存夾（不留半成品空夾）。"""
    tmp = target.with_name(f".{target.name}.tmp-{os.getpid()}")
    if tmp.exists():
        shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    try:
        yield tmp
    except BaseException:
        shutil.rmtree(tmp, ignore_errors=True)
        raise
    if target.exists():
        raise FileExistsError(f"{target} 已存在；先用 retire_to_stale 搬走舊版")
    os.replace(tmp, target)
