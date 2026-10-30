#!/usr/bin/env python3
"""Step 5：把各訂單的 taxonomy_decision_table.csv 彙整成結果摘要表。

輸出：{OUTPUT}/YYYY-MM-DD_new_orders_taxonomy_summary.csv（UTF-8 with BOM）
欄位沿用既有成品：訂單,樣本,Marker,候選菌種,Top hit,Identity (%),Query cover (%),Review flag,建議動作

用法：python3 -m pipeline.summarize --orders S2604280595 S2605070135 [--config ...] [--output-dir DIR]
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

from pipeline.config import gene_lookup_path, load_config, output_dir
from pipeline.consensus_call import load_lookup

HEADER = [
    "訂單", "樣本", "Marker", "候選菌種", "Top hit",
    "Identity (%)", "Query cover (%)", "Review flag", "建議動作",
]

# SKILL.md Step 5 的 suggested_action 說明
ACTION_MEANING = {
    "CANDIDATE_OK_IF_SEQUENCE_QC_ACCEPTED": "可信，確認 QC 後接受",
    "SPECIES_NEEDS_MANUAL_TOP_HIT_REVIEW": "需手動確認 top hit",
    "CHECK_TAXONOMY_SYNONYMS_OR_ADD_MARKER": "可能有分類學同義詞，建議加跑其他 marker",
    "CHECK_TOP_HITS_AND_CONSIDER_ADDITIONAL_MARKER": "identity 偏低，建議加跑 marker",
    "CHECK_ASSEMBLY_OR_CHROMATOGRAM_FIRST": "Assembly 有問題，先查色譜圖",
    "HITLIST_SATURATED_CHECK_TYPE_SEQUENCES": "比對清單被高品質 hit 塞滿，候選可能不完整，請看 type 序列結果",
    "TYPE_SUPPORT_MISSING_REVIEW": "全庫有候選但沒有 type 序列支持，需人工覆核（降信心，不代表錯）",
    "ITS_GENUS_LEVEL_ONLY": "此屬 ITS 無種級解析力，只能報到屬／section，種名待第二 marker",
}


def today_taipei() -> str:
    return datetime.now(ZoneInfo("Asia/Taipei")).strftime("%Y-%m-%d")


def order_of_dir(name: str) -> str:
    return name.removesuffix("_ncbi_blast")


def downgrade_action(r: dict[str, str], lookup: dict | None) -> str:
    """🔴-3：ITS 在 its_sufficient=false 的屬，不可輸出「可信，確認 QC 後接受」。"""
    action = r["suggested_action"]
    if (
        lookup
        and action == "CANDIDATE_OK_IF_SEQUENCE_QC_ACCEPTED"
        and r["marker"] == "ITS"
        and (lookup.get(r.get("candidate_genus", "")) or {}).get("its_sufficient") is False
    ):
        return "ITS_GENUS_LEVEL_ONLY"
    return action


def summarize_orders(out: Path, orders: list[str], lookup: dict | None = None) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for order in orders:
        path = out / f"{order}_ncbi_blast" / "taxonomy_decision_table.csv"
        if not path.is_file():
            continue
        with path.open(newline="", encoding="utf-8-sig") as handle:
            for r in csv.DictReader(handle):
                rows.append(
                    {
                        "訂單": order,
                        "樣本": r["sample"],
                        "Marker": r["marker"],
                        "候選菌種": r["candidate_species"],
                        "Top hit": r["top_hit_organism"],
                        "Identity (%)": r["top_pident"],
                        "Query cover (%)": r["top_qcov_pct"],
                        "Review flag": r["review_flag"],
                        "建議動作": downgrade_action(r, lookup),
                    }
                )
    return rows


def write_summary(rows: list[dict[str, str]], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=HEADER)
        writer.writeheader()
        writer.writerows(rows)
    return path


def discover_orders(out: Path) -> list[str]:
    return sorted(
        order_of_dir(p.name)
        for p in out.glob("*_ncbi_blast")
        if (p / "taxonomy_decision_table.csv").is_file()
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orders", nargs="+", help="訂單號；省略 = 輸出夾內全部已有 BLAST 結果的訂單")
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None, help="結果寫到這裡；預設 = config 的輸出夾")
    parser.add_argument("--input-dir", type=Path, default=None, help="讀 BLAST 結果的夾；預設 = 輸出夾")
    parser.add_argument("--date", default=None, help="檔名日期 YYYY-MM-DD，預設今天（台北）")
    args = parser.parse_args()
    cfg = load_config(args.config)
    out = args.output_dir or output_dir(cfg)
    src = args.input_dir or out
    orders = args.orders or discover_orders(src)
    rows = summarize_orders(src, orders, load_lookup(gene_lookup_path(cfg)))
    path = write_summary(rows, out / f"{args.date or today_taipei()}_new_orders_taxonomy_summary.csv")
    print(f"完成，{len(orders)} 個訂單、{len(rows)} 列，輸出至 {path}")
    for row in rows:
        meaning = ACTION_MEANING.get(row["建議動作"], "")
        print(f"{row['訂單']} | {row['樣本']} | {row['Marker']} | {row['候選菌種']} | {meaning or row['建議動作']}")


if __name__ == "__main__":
    main()
