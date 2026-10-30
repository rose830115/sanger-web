#!/usr/bin/env python3
"""Step 1：比對 rawdata 與輸出夾，找出尚未處理完的批次（規則見 SKILL.md Step 1）。

- rawdata 子資料夾格式 {OrderID}_{BatchID}；同 OrderID 的子批次合併處理
- 沒有 {subfolder}_processed/                              -> 需要 QC
- 沒有 {OrderID}_combined_assembly/                        -> 需要 Assembly
- 沒有 {OrderID}_ncbi_blast/taxonomy_decision_table.csv    -> 需要 BLAST

用法：python3 -m pipeline.scan_batches [--config config.yaml] [--json]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from pipeline.config import load_config, output_dir, rawdata_dir
from pipeline.stage_guard import stage_report


def order_id_of(subfolder: str) -> str:
    """`_` 前的部分為訂單號。"""
    return subfolder.split("_", 1)[0]


def scan(raw: Path, out: Path, cfg: dict | None = None) -> dict:
    """回傳 {"orders": {OrderID: {...}}, "pending_*": [...]}。

    給 cfg：除了「輸出存不存在」，還用輸入清單雜湊判斷上游是否已變動（修正 🔴-6），
    變動的階段列入 stale_*，視為待處理。不給 cfg：維持只看存在與否（舊行為）。
    """
    if not raw.is_dir():
        raise SystemExit(f"找不到 rawdata 目錄：{raw}（NAS 是否已掛載？）")
    subfolders = sorted(
        p.name for p in raw.iterdir() if p.is_dir() and not p.name.startswith(".")
    )
    orders: dict[str, dict] = {}
    for sub in subfolders:
        info = orders.setdefault(
            order_id_of(sub),
            {"subfolders": [], "need_qc": [], "need_assembly": False, "need_blast": False},
        )
        info["subfolders"].append(sub)
        if not (out / f"{sub}_processed").is_dir():
            info["need_qc"].append(sub)
    for order_id, info in orders.items():
        info["need_assembly"] = not (out / f"{order_id}_combined_assembly").is_dir()
        info["need_blast"] = not (
            out / f"{order_id}_ncbi_blast" / "taxonomy_decision_table.csv"
        ).is_file()
        info["stale_qc"], info["stale_assembly"], info["stale_blast"] = [], False, False
        if cfg is not None:
            rep = stage_report(raw, out, order_id, info, cfg)
            info["stale_qc"] = [s for s, st in rep["qc"].items() if st == "stale"]
            info["stale_assembly"] = rep["assembly"] == "stale"
            info["stale_blast"] = rep["blast"] == "stale"
        info["pending"] = bool(
            info["need_qc"] or info["need_assembly"] or info["need_blast"]
            or info["stale_qc"] or info["stale_assembly"] or info["stale_blast"]
        )

    pending_orders = [o for o, i in orders.items() if i["pending"]]
    return {
        "orders": orders,
        "pending_orders": pending_orders,
        "pending_subfolders": [s for o in pending_orders for s in orders[o]["subfolders"]],
    }


def scan_from_config(cfg: dict) -> dict:
    return scan(rawdata_dir(cfg), output_dir(cfg), cfg)


def format_report(result: dict) -> str:
    lines = [
        f"偵測到 {len(result['pending_subfolders'])} 個新子批次，"
        f"屬於 {len(result['pending_orders'])} 個訂單，準備開始處理。"
    ]
    for order_id in result["pending_orders"]:
        info = result["orders"][order_id]
        todo = []
        if info["need_qc"]:
            todo.append(f"QC x{len(info['need_qc'])}")
        if info.get("stale_qc"):
            todo.append(f"QC 重做 x{len(info['stale_qc'])}（輸入已變）")
        if info["need_assembly"]:
            todo.append("Assembly")
        elif info.get("stale_assembly"):
            todo.append("Assembly 重做（輸入已變）")
        if info["need_blast"]:
            todo.append("BLAST")
        elif info.get("stale_blast"):
            todo.append("BLAST 重新解析（上游已變）")
        lines.append(f"- {order_id}（{', '.join(info['subfolders'])}）：{' / '.join(todo)}")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--json", action="store_true", help="以 JSON 輸出")
    args = parser.parse_args()
    result = scan_from_config(load_config(args.config))
    if args.json:
        print(json.dumps(result, ensure_ascii=False, indent=2))
    else:
        print(format_report(result))


if __name__ == "__main__":
    main()
