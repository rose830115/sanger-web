#!/usr/bin/env python3
"""一鍵串全流程（SKILL.md 七步），可從任一步續跑。

步驟（--from-step 可選）：
  scan      掃描 rawdata 找未處理批次
  qc        逐子批次 QC
  assembly  逐訂單 Assembly
  blast     逐訂單 BLAST（線上；--offline 則只重解析既有 XML）
  summarize 結果摘要表（Step 5）
  consensus 多 marker 交集判讀（Step 6）
  cross     跨批次菌種候選彙整（Step 7）

範例：
  python3 run_all.py                          # 偵測並處理所有未完成批次
  python3 run_all.py --from-step blast        # 從 BLAST 續跑（QC／Assembly 之前已完成）
  python3 run_all.py --from-step consensus --orders S2604280595 S2605070135
  python3 run_all.py --dry-run                # 只看會處理什麼

每步都是「缺什麼補什麼」：輸出有輸入雜湊紀錄（{訂單}_input_manifest.json），上游（ab1／設定）一變，下游自動失效重做，
舊輸出搬進 _stale/ 不刪除；--force 強制重做 QC／Assembly；BLAST XML 逐條驗證序列雜湊，對不上才重送。
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path
from typing import Callable

from pipeline import consensus_call, cross_batch_summary, sanger_pair_assemble, sanger_qc, summarize
from pipeline.config import blast_settings, gene_lookup_path, load_config, output_dir, rawdata_dir
from pipeline.sanger_ncbi_blast import run_blast
from pipeline.scan_batches import format_report, scan
from pipeline.stage_guard import (
    atomic_dir,
    blast_fingerprint,
    load_manifest,
    retire_to_stale,
    save_manifest,
    stage_report,
)

STEPS = ["scan", "qc", "assembly", "blast", "summarize", "consensus", "cross"]


def log(msg: str) -> None:
    print(msg, flush=True)


def run_pipeline(
    cfg: dict,
    from_step: str = "scan",
    to_step: str = "cross",
    orders: list[str] | None = None,
    raw: Path | None = None,
    out: Path | None = None,
    force: bool = False,
    offline: bool = False,
    dry_run: bool = False,
    date: str | None = None,
    progress: Callable[[str, str], None] | None = None,
) -> dict:
    """執行流程，回傳 {"failures": [...], "orders": [...]}。

    progress(step, detail)：每個工作單位開始前呼叫（給網頁進度條用），預設不做事。
    """
    notify = progress or (lambda step, detail: None)
    raw = raw or rawdata_dir(cfg)
    out = out or output_dir(cfg)
    first, last = STEPS.index(from_step), STEPS.index(to_step)
    active = STEPS[first : last + 1]
    failures: list[str] = []

    result = scan(raw, out, cfg)
    log(format_report(result))
    todo = list(orders) if orders else list(result["pending_orders"])
    unknown = [o for o in todo if o not in result["orders"]]
    if unknown:
        raise SystemExit(f"rawdata 找不到訂單：{', '.join(unknown)}")
    if not todo and any(s in active for s in STEPS[1:]):
        log("沒有待處理訂單；要重跑摘要／判讀請用 --orders 指定訂單。")
        return {"failures": failures, "orders": todo}
    if dry_run:
        log(f"[dry-run] 將處理訂單：{', '.join(todo)}；步驟：{' -> '.join(active)}")
        return {"failures": failures, "orders": todo}

    blast_cfg = blast_settings(cfg)
    if "blast" in active and not offline and not blast_cfg.get("email"):
        raise SystemExit(
            "線上 BLAST 需要聯絡信箱：請在 config.yaml 的 blast.email 填實驗室信箱，"
            "或設環境變數 SANGER_NCBI_EMAIL；只想離線重解析請加 --offline。"
        )

    # 每個訂單先算一次「目前輸入指紋」與各階段狀態（修正 🔴-6）；上游變動 -> 下游失效
    reports = {o: stage_report(raw, out, o, result["orders"][o], cfg) for o in todo}

    def remember(order: str, key: str, value) -> None:
        man = load_manifest(out, order)
        if key == "qc":
            man.setdefault("qc", {}).update(value)
        else:
            man[key] = value
        man["files"] = reports[order]["current"]["files"]
        save_manifest(out, order, man)

    if "qc" in active:
        for order in todo:
            for sub in result["orders"][order]["subfolders"]:
                target = out / f"{sub}_processed"
                state = reports[order]["qc"].get(sub, "missing")
                fp = reports[order]["current"]["qc"].get(sub)
                if state == "legacy" and not force:
                    log(f"[qc] {sub} 沿用既有輸出（舊輸出沒有輸入紀錄，無法驗證；補記目前輸入）")
                    remember(order, "qc", {sub: {"fp": fp}})
                    continue
                if state == "fresh" and not force:
                    continue
                log(f"[qc] {sub}" + ("（輸入已變，重做）" if state == "stale" else ""))
                notify("qc", sub)
                try:
                    if target.exists():
                        retire_to_stale(out, target)
                    with atomic_dir(target) as tmp:
                        sanger_qc.process(raw / sub, tmp, cfg["qc"])
                    remember(order, "qc", {sub: {"fp": fp}})
                except (SystemExit, Exception) as exc:  # noqa: BLE001
                    failures.append(f"qc {sub}: {exc}")
                    log(f"  失敗：{exc}")

    if "assembly" in active:
        for order in todo:
            target = out / f"{order}_combined_assembly"
            state = reports[order]["assembly"]
            fp = reports[order]["current"]["assembly"]
            if state == "legacy" and not force:
                log(f"[assembly] {order} 沿用既有輸出（舊輸出沒有輸入紀錄，無法驗證；補記目前輸入）")
                remember(order, "assembly", {"fp": fp})
                continue
            if state == "fresh" and not force:
                continue
            log(f"[assembly] {order}" + ("（輸入已變，重做）" if state == "stale" else ""))
            notify("assembly", order)
            try:
                dirs = [raw / sub for sub in result["orders"][order]["subfolders"]]
                if target.exists():
                    retire_to_stale(out, target)
                with atomic_dir(target) as tmp:
                    sanger_pair_assemble.process(dirs, tmp, cfg["qc"], cfg["assembly"])
                remember(order, "assembly", {"fp": fp})
            except (SystemExit, Exception) as exc:  # noqa: BLE001
                failures.append(f"assembly {order}: {exc}")
                log(f"  失敗：{exc}")

    if "blast" in active:
        pending_blast = []
        for o in todo:
            asm_dir = out / f"{o}_combined_assembly"
            fp_now = blast_fingerprint(asm_dir, cfg)          # Assembly 剛重做過就會是新的指紋
            table = (out / f"{o}_ncbi_blast" / "taxonomy_decision_table.csv").is_file()
            recorded = (load_manifest(out, o).get("blast") or {}).get("fp")
            if not table:
                pending_blast.append(o)
            elif recorded is None:
                log(f"[blast] {o} 沿用既有解析表（舊輸出沒有輸入紀錄；補記目前輸入）")
                if fp_now:
                    remember(o, "blast", {"fp": fp_now})
            elif recorded != fp_now:
                log(f"[blast] {o} 上游或設定已變，重新解析（XML 逐條驗證，對不上的才重送）")
                pending_blast.append(o)
        for index, order in enumerate(pending_blast):
            asm_dir = out / f"{order}_combined_assembly"
            fasta = asm_dir / "paired_consensus_for_blast.fasta"
            if not fasta.is_file():
                failures.append(f"blast {order}: 缺 {fasta.name}（Assembly 未完成？）")
                log(f"[blast] {order} 略過：缺 {fasta}")
                continue
            log(f"[blast] {order}")
            notify("blast", order)
            try:
                run_blast(
                    fasta,
                    out / f"{order}_ncbi_blast",
                    blast_cfg,
                    assembly_summary=asm_dir / "assembly_summary.csv",
                    offline=offline,
                    on_event=lambda msg, _o=order: notify("blast", f"{_o}　{msg}"),
                )
                remember(order, "blast", {"fp": blast_fingerprint(asm_dir, cfg)})
            except (SystemExit, Exception) as exc:  # noqa: BLE001
                failures.append(f"blast {order}: {exc}")
                log(f"  失敗：{exc}")
            if not offline and index < len(pending_blast) - 1:
                time.sleep(blast_cfg["delay_between_orders_sec"])

    done_orders = [
        o for o in todo if (out / f"{o}_ncbi_blast" / "taxonomy_decision_table.csv").is_file()
    ]
    stamp = date or summarize.today_taipei()
    if "summarize" in active:
        notify("summarize", "")
        rows = summarize.summarize_orders(out, done_orders, consensus_call.load_lookup(gene_lookup_path(cfg)))
        path = summarize.write_summary(rows, out / f"{stamp}_new_orders_taxonomy_summary.csv")
        log(f"[summarize] {len(rows)} 列 -> {path}")
    if "consensus" in active:
        notify("consensus", "")
        rows = consensus_call.build_rows(out, done_orders, cfg)
        path = consensus_call.write_rows(
            rows, out / f"{stamp}_new_orders_marker_intersection_summary.csv"
        )
        log(f"[consensus] -> {path}")
        log(consensus_call.report(rows))
    if "cross" in active:
        notify("cross", "")
        log("[cross]")
        cross_batch_summary.run(cfg, out)

    if failures:
        log("\n以下步驟失敗，需處理後從該步續跑：")
        for item in failures:
            log(f"- {item}")
    return {"failures": failures, "orders": todo}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--from-step", choices=STEPS, default="scan")
    parser.add_argument("--to-step", choices=STEPS, default="cross")
    parser.add_argument("--orders", nargs="+", help="只處理這些訂單；省略 = 自動偵測未完成的")
    parser.add_argument("--rawdata-dir", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--force", action="store_true", help="QC／Assembly 即使已有輸出也重做")
    parser.add_argument("--offline", action="store_true", help="BLAST 只重解析既有 XML，不連線 NCBI")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--date", default=None, help="摘要檔名日期 YYYY-MM-DD，預設今天（台北）")
    args = parser.parse_args()
    if STEPS.index(args.from_step) > STEPS.index(args.to_step):
        parser.error("--from-step 不能晚於 --to-step")
    result = run_pipeline(
        load_config(args.config),
        args.from_step,
        args.to_step,
        args.orders,
        args.rawdata_dir,
        args.output_dir,
        args.force,
        args.offline,
        args.dry_run,
        args.date,
    )
    sys.exit(1 if result["failures"] else 0)


if __name__ == "__main__":
    main()
