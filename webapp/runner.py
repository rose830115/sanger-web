"""網頁介面背後的邏輯（不依賴 streamlit，方便單元測試）。

- 上傳 zip 的驗證與解壓（檔名格式、zip slip、大小上限）
- 組出 pipeline 要的 rawdata／output 工作夾，呼叫 run_all.run_pipeline
- 讀回結果並整理成「每樣本一列」的結果表
- 打包結果 zip
- 安全閘：輸出絕不可落在 NAS 根目錄底下
"""

from __future__ import annotations

import csv
import hashlib
import io
import re
import shutil
import tempfile
import time
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

from pipeline import consensus_call
from pipeline.config import PROJECT_ROOT, load_config
from pipeline.markers import UNKNOWN_MARKER, marker_for_primer
from run_all import run_pipeline

MAX_UNCOMPRESSED_BYTES = 400 * 1024 * 1024   # 解壓後上限（Community Cloud 資源有限）
MAX_FILES = 2000
FILENAME_EXAMPLE = "S260428059517_63aSID322a_ITS4_01.ab1"
WORK_BASE = Path(tempfile.gettempdir()) / "sanger_web_jobs"
MINUTES_PER_MARKER_ONLINE = 4             # 線上 BLAST 每個 marker 批次的粗估分鐘數
NO_SUGGESTION = "—"


class UploadError(Exception):
    """給使用者看的錯誤訊息（繁體中文）。"""


# ---------------------------------------------------------------- 安全閘
def assert_not_under_nas(path: Path, cfg: dict) -> None:
    """輸出夾若在 NAS 根目錄底下就直接擋下（NAS 只准讀）。

    nas_root 留空（公開 repo 的預設）時不可把空字串解析成專案根目錄，
    否則所有本機輸出都會被誤擋；此時仍保留 /Volumes/ 的通用防線。
    """
    from pipeline.config import _resolve  # noqa: PLC0415

    target = Path(path).expanduser().resolve()
    if str(target).startswith("/Volumes/"):
        raise UploadError(f"輸出位置不可在 NAS（只准讀取）：{target}")
    nas_root = cfg["paths"].get("nas_root") or ""
    if nas_root:
        nas = _resolve(cfg, nas_root).resolve()
        if target == nas or nas in target.parents:
            raise UploadError(f"輸出位置不可在 NAS（只准讀取）：{target}")


# ---------------------------------------------------------------- 上傳驗證
def check_ab1_name(stem: str) -> str | None:
    """回傳錯誤原因；合格回傳 None。格式：流水號_樣本_引子_孔位。"""
    parts = stem.split("_")
    if len(parts) < 4:
        return f"底線分隔不足 4 段（現在 {len(parts)} 段）"
    if any(not p.strip() for p in parts[:4]):
        return "有空白的欄位（出現連續底線或開頭結尾底線）"
    return None


def _safe_name(name: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "-", name).strip("-.")
    return cleaned or "upload"


@dataclass
class UploadReport:
    subfolders: dict[str, list[str]] = field(default_factory=dict)   # 子批次 -> ab1 檔名
    unknown_primers: dict[str, int] = field(default_factory=dict)    # 引子 -> 檔數
    ignored: int = 0

    @property
    def n_files(self) -> int:
        return sum(len(v) for v in self.subfolders.values())

    @property
    def orders(self) -> list[str]:
        return sorted({s.split("_", 1)[0] for s in self.subfolders})


def extract_uploads(files: list[tuple[str, bytes]], raw_dir: Path) -> UploadReport:
    """把一個或多個 zip 解到 raw_dir/<子批次>/*.ab1；不合格就丟 UploadError（一次列出全部問題）。"""
    report = UploadReport()
    bad: list[str] = []
    total = 0
    count = 0
    pending: list[tuple[str, str, bytes]] = []   # (子批次, 檔名, 內容)

    for zip_name, data in files:
        try:
            zf = zipfile.ZipFile(io.BytesIO(data))
        except zipfile.BadZipFile as exc:
            raise UploadError(f"「{zip_name}」不是有效的 zip 檔，請重新下載後再上傳。") from exc
        stem_zip = _safe_name(Path(zip_name).stem)
        for info in zf.infolist():
            if info.is_dir():
                continue
            parts = [p for p in Path(info.filename.replace("\\", "/")).parts if p not in ("/", "")]
            if (
                not parts
                or ".." in parts
                or Path(info.filename).is_absolute()
                or parts[0] == "__MACOSX"
                or parts[-1].startswith(".")
            ):
                report.ignored += 1
                continue
            fname = parts[-1]
            if not fname.lower().endswith(".ab1"):
                report.ignored += 1
                continue
            count += 1
            total += info.file_size
            if count > MAX_FILES or total > MAX_UNCOMPRESSED_BYTES:
                raise UploadError("上傳內容太大（檔案數或解壓後容量超過上限），請分批上傳。")
            reason = check_ab1_name(Path(fname).stem)
            if reason:
                bad.append(f"{zip_name} → {fname}：{reason}")
                continue
            sub = _safe_name(parts[-2]) if len(parts) >= 2 else stem_zip
            pending.append((sub, fname, zf.read(info)))

    if bad:
        shown = "\n".join(f"- {b}" for b in bad[:15])
        more = f"\n…還有 {len(bad) - 15} 個" if len(bad) > 15 else ""
        raise UploadError(
            "以下檔名不符合格式，請確認是定序公司給的原始檔名（沒有被改過）。\n"
            f"正確格式：流水號_樣本名_引子_孔位.ab1，例如 {FILENAME_EXAMPLE}\n{shown}{more}"
        )
    if not pending:
        raise UploadError("zip 裡找不到任何 .ab1 檔。請確認上傳的是定序公司的結果壓縮檔。")

    for sub, fname, content in pending:
        target = raw_dir / sub
        target.mkdir(parents=True, exist_ok=True)
        (target / fname).write_bytes(content)
        report.subfolders.setdefault(sub, []).append(fname)
        primer = fname[:-4].split("_")[2]
        if marker_for_primer(primer) == UNKNOWN_MARKER:
            report.unknown_primers[primer] = report.unknown_primers.get(primer, 0) + 1
    return report


def job_id_for(files: list[tuple[str, bytes]]) -> str:
    """同樣的檔案內容 -> 同一個工作夾，重新整理後可斷點續跑。"""
    h = hashlib.sha1()
    for name, data in sorted(files, key=lambda x: x[0]):
        h.update(name.encode())
        h.update(hashlib.sha1(data).digest())
    return h.hexdigest()[:12]


def clear_job(job_id: str) -> bool:
    """「清除重跑」：刪掉這個工作夾（只限 WORK_BASE 底下、且不在 NAS）。回傳有沒有刪到東西。"""
    base = WORK_BASE.resolve()
    target = (WORK_BASE / job_id).resolve()
    if target.parent != base or not re.fullmatch(r"[0-9a-f]{6,40}", job_id):
        raise UploadError("工作夾位置不合法，拒絕清除。")
    if not target.is_dir():
        return False
    shutil.rmtree(target, ignore_errors=True)
    return True


def cleanup_old_jobs(max_age_hours: float = 24) -> None:
    if not WORK_BASE.is_dir():
        return
    cutoff = time.time() - max_age_hours * 3600
    for d in WORK_BASE.iterdir():
        try:
            if d.is_dir() and d.stat().st_mtime < cutoff:
                shutil.rmtree(d, ignore_errors=True)
        except OSError:
            pass


# ---------------------------------------------------------------- 預估時間
def estimate_minutes(report: UploadReport, offline: bool) -> int:
    if offline:
        return 1
    markers = {
        marker_for_primer(f[:-4].split("_")[2])
        for names in report.subfolders.values()
        for f in names
    }
    return max(1, len(markers) * len(report.orders) * MINUTES_PER_MARKER_ONLINE)


# ---------------------------------------------------------------- 執行
STEP_LABELS = {
    "qc": "檢查定序品質（QC）",
    "assembly": "正反向合併（Assembly）",
    "blast": "NCBI 比對（BLAST）",
    "summarize": "整理結果摘要",
    "consensus": "多 marker 交集判讀",
    "cross": "整理菌種總表",
}


def run_job(
    cfg: dict,
    raw_dir: Path,
    out_dir: Path,
    orders: list[str],
    offline: bool,
    on_progress: Callable[[float, str], None] | None = None,
    date: str | None = None,
) -> dict:
    """跑完整流程（QC -> Assembly -> BLAST -> 摘要 -> 交集 -> 總表）。回傳 run_pipeline 的結果。"""
    assert_not_under_nas(out_dir, cfg)
    out_dir.mkdir(parents=True, exist_ok=True)
    from pipeline.scan_batches import scan  # noqa: PLC0415

    info = scan(raw_dir, out_dir, cfg)["orders"]
    n_sub = sum(len(info[o]["subfolders"]) for o in orders if o in info)
    total = max(1, n_sub + 2 * len(orders) + 3)
    done = {"n": 0}

    def progress(step: str, detail: str) -> None:
        if on_progress:
            on_progress(min(done["n"] / total, 0.99), f"{STEP_LABELS.get(step, step)} {detail}".strip())
        done["n"] += 1

    result = run_pipeline(
        cfg,
        orders=orders,
        raw=raw_dir,
        out=out_dir,
        offline=offline,
        date=date,
        progress=progress,
    )
    if on_progress:
        on_progress(1.0, "完成")
    return result


def seed_blast_xml(out_dir: Path, orders: list[str], xml_files: list[tuple[str, bytes]]) -> int:
    """離線重解析用：把使用者提供的 BLAST 原始 XML 放進各訂單的 raw_blast_xml。"""
    n = 0
    for order in orders:
        target = out_dir / f"{order}_ncbi_blast" / "raw_blast_xml"
        target.mkdir(parents=True, exist_ok=True)
        for name, data in xml_files:
            (target / _safe_name(Path(name).name)).write_bytes(data)
            n += 1
    return n


def copy_blast_xml_from(src_dir: Path, out_dir: Path, orders: list[str]) -> int:
    """本機模式：從既有結果夾（唯讀）複製 raw_blast_xml；來源只讀不寫。"""
    n = 0
    for order in orders:
        src = Path(src_dir) / f"{order}_ncbi_blast" / "raw_blast_xml"
        if not src.is_dir():
            continue
        target = out_dir / f"{order}_ncbi_blast" / "raw_blast_xml"
        target.mkdir(parents=True, exist_ok=True)
        for x in src.glob("*.xml"):
            shutil.copy2(x, target / x.name)
            n += 1
    return n


# ---------------------------------------------------------------- 結果整理
def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def _find_one(out: Path, pattern: str) -> Path | None:
    hits = sorted(out.glob(pattern))
    return hits[-1] if hits else None


def format_primer_pair(primers: str) -> str:
    """'rpb2: 5F GAYG... / 7CR CCCA...' -> '5F / 7CR'（只取引子名；序列留在詳細欄）。"""
    out = []
    for chunk in primers.split("；"):
        body = chunk.split(":", 1)[-1]
        names = [seg.strip().split(" ")[0] for seg in body.split(" / ") if seg.strip()]
        if names:
            out.append(" / ".join(names))
    return "；備用 ".join(out)


def build_result_rows(out: Path, orders: list[str]) -> dict:
    """讀輸出夾，回傳 {"rows": [...], "warnings": [...], "counts": {...}}。每樣本一列。"""
    consensus_path = _find_one(out, "*_new_orders_marker_intersection_summary.csv")
    consensus = _read_csv(consensus_path) if consensus_path else []
    rows: list[dict[str, str]] = []
    seen: set[tuple[str, str]] = set()
    warnings: list[str] = []

    markers_by_sample: dict[str, list[str]] = {}
    unknown: dict[str, set[str]] = {}
    merged_orders = {m for r in consensus for m in (r.get("合併訂單") or "").split("、") if m}
    for order in list(orders) + sorted(merged_orders - set(orders)):     # 跨訂單合併的歷史單也要讀 marker
        for r in _read_csv(out / f"{order}_ncbi_blast" / "taxonomy_decision_table.csv"):
            markers_by_sample.setdefault(r["sample"], []).append(r["marker"])
        asm = _read_csv(out / f"{order}_combined_assembly" / "assembly_summary.csv")
        for r in asm:
            if r["marker"] == UNKNOWN_MARKER:
                unknown.setdefault(order, set()).add(r["sample"])

    for r in consensus:
        key = (r["訂單"], r["樣本"])
        seen.add(key)
        for merged in (r.get("合併訂單") or "").split("、"):
            if merged:
                seen.add((merged, r["樣本"]))
        if r["判讀狀態"] == consensus_call.STATUS_UNIQUE:
            species, level = r["唯一物種"], "ok"
        elif r["判讀狀態"] == consensus_call.STATUS_UNIQUE_PENDING:
            species, level = f"{r['唯一物種']}（待確認）", "warn"
        elif r["判讀狀態"] == consensus_call.STATUS_INTERSECT_PENDING:
            species, level = f"{r['唯一物種']}（交集唯一，待覆核）", "warn"
        elif r["判讀狀態"] == consensus_call.STATUS_NOT_BLASTED:
            species, level = "（尚未比對：BLAST 沒跑完或結果與目前序列對不上）", "warn"
        elif r["判讀狀態"] == consensus_call.STATUS_GENUS_ONLY:
            species, level = f"{r['判讀屬']}（僅到屬／section，種名待第二 marker）", "warn"
        elif r["共同候選交集"]:
            species, level = r["共同候選交集"], "warn"
        else:
            species, level = r["各 marker 高品質候選"] or "（無高品質候選）", "warn"
        marker = r["下次建議 marker"].strip()
        primers = r["建議引子"].strip()
        suggestion = NO_SUGGESTION
        pair = ""
        if marker:
            suggestion = marker
            pair = format_primer_pair(primers) if primers else ""
        mk = "、".join(dict.fromkeys(markers_by_sample.get(key[1], [])))
        note = ""
        if key[1] in unknown.get(key[0], set()):
            note = "引子名未辨識，marker 歸為「其他」，請人工確認"
            level = "warn"
        rows.append(
            {
                "訂單": r["訂單"],
                "樣本": r["樣本"],
                "Marker": mk or r["Marker 數"],
                "判讀狀態": r["判讀狀態"],
                "物種／候選": species,
                "下次建議 Marker": suggestion,
                "建議引子對": pair,
                "引子序列": primers,
                "備註": "；".join(x for x in [r["備註"], note] if x),
                "等級": level,
            }
        )

    for order in orders:
        for r in _read_csv(out / f"{order}_combined_assembly" / "assembly_summary.csv"):
            if (order, r["sample"]) in seen:
                continue
            seen.add((order, r["sample"]))
            note = "引子名未辨識，marker 歸為「其他」" if r["marker"] == UNKNOWN_MARKER else ""
            rows.append(
                {
                    "訂單": order,
                    "樣本": r["sample"],
                    "Marker": r["marker"],
                    "判讀狀態": f"未能比對（Assembly：{r['assembly_status']}）",
                    "物種／候選": "（無結果）",
                    "下次建議 Marker": "重新定序或人工檢查色譜圖",
                    "建議引子對": "",
                    "引子序列": "",
                    "備註": note,
                    "等級": "fail",
                }
            )

    for order, samples in unknown.items():
        warnings.append(
            f"訂單 {order}：{len(samples)} 個樣本的引子名程式不認得，marker 歸為「其他」（資料都有保留，請人工確認）。"
        )
    counts = {
        "samples": len(rows),
        "ok": sum(r["等級"] == "ok" for r in rows),
        "warn": sum(r["等級"] == "warn" for r in rows),
        "fail": sum(r["等級"] == "fail" for r in rows),
    }
    return {"rows": rows, "warnings": warnings, "counts": counts}


RESULT_COLUMNS = [
    "訂單", "樣本", "Marker", "判讀狀態", "物種／候選", "下次建議 Marker", "建議引子對", "引子序列", "備註",
]


def rows_to_csv_bytes(rows: list[dict[str, str]]) -> bytes:
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=RESULT_COLUMNS, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)
    return buf.getvalue().encode("utf-8-sig")


def build_result_zip(out: Path, rows: list[dict[str, str]]) -> bytes:
    """打包：結果總表 + 輸出夾內全部 CSV／FASTA／txt／BLAST XML。"""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("結果總表.csv", rows_to_csv_bytes(rows))
        for p in sorted(out.rglob("*")):
            if p.is_file() and p.suffix.lower() in {".csv", ".fasta", ".txt", ".xml"}:
                zf.write(p, p.relative_to(out).as_posix())
    return buf.getvalue()


def new_job_dir(job_id: str) -> Path:
    d = WORK_BASE / job_id
    d.mkdir(parents=True, exist_ok=True)
    return d


def app_settings(cfg: dict | None = None) -> dict:
    cfg = cfg or load_config()
    base = {"enable_local_mode": False, "local_output_dir": "local_output"}
    base.update(cfg.get("app") or {})
    return base


LOCAL_OUTPUT_DEFAULT = PROJECT_ROOT / "local_output"
