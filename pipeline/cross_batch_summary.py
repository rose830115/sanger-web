#!/usr/bin/env python3
"""Step 7：跨批次菌種候選彙整（備查，不作最終唯一物種判讀）。

以 SKILL.md Step 7 內嵌程式為基底，路徑／門檻改吃 config.yaml。
掃描全輸出夾的 */*_ncbi_blast/blast_top100_hits.csv；條件 mismatch <= 10；
名稱含 clone / gene / genes / genomic（全詞）-> 「屬名 sp.」；sp / sp. 統一成 sp.。
輸出 {OUTPUT}/blast_species_summary.csv（UTF-8 with BOM）。

用法：python3 -m pipeline.cross_batch_summary [--config ...] [--output-dir DIR]
"""

from __future__ import annotations

import argparse
import csv
import glob
import re
from collections import OrderedDict
from pathlib import Path

from pipeline.config import load_config, output_dir


def normalize_species(name: str | None, keywords: list[str]) -> str | None:
    if not name:
        return None
    name_lower = name.lower()
    for kw in keywords:
        if re.search(r"\b" + kw + r"\b", name_lower):
            genus = name.strip().split()[0]
            return f"{genus} sp."
    return re.sub(r"\bsp\.?\s*$", "sp.", name.strip(), flags=re.IGNORECASE).strip()


def collect(out: Path, keywords: list[str], mismatch_limit: int) -> "OrderedDict[str, list[str]]":
    files = sorted(glob.glob(f"{glob.escape(str(out))}/*_ncbi_blast/blast_top100_hits.csv"))
    results: "OrderedDict[str, list[str]]" = OrderedDict()
    for fpath in files:
        with open(fpath, newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                try:
                    if int(row["mismatches"]) > mismatch_limit:
                        continue
                except (ValueError, KeyError):
                    continue
                sample = row.get("sample", "").strip()
                species = normalize_species(row.get("organism_guess", "").strip(), keywords)
                if not sample or not species:
                    continue
                if sample not in results:
                    results[sample] = []
                if species.lower() not in [s.lower() for s in results[sample]]:
                    results[sample].append(species)
    return results


def write_summary(results: "OrderedDict[str, list[str]]", path: Path, mismatch_limit: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow(["樣本名稱", f"可能菌種（mismatch ≤ {mismatch_limit}）", "菌種數量"])
        for sample, sp_list in results.items():
            writer.writerow([sample, " | ".join(sp_list), len(sp_list)])


def run(cfg: dict, out: Path | None = None, src: Path | None = None) -> Path:
    sc = cfg["species_summary"]
    out = out or output_dir(cfg)
    results = collect(src or out, sc["genus_sp_keywords"], sc["mismatch_limit"])
    path = out / sc["output_name"]
    write_summary(results, path, sc["mismatch_limit"])
    print(f"完成，共 {len(results)} 個樣本，輸出至 {path}")
    return path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None, help="結果寫到這裡；預設 = config 的輸出夾")
    parser.add_argument("--input-dir", type=Path, default=None, help="讀 blast_top100_hits.csv 的夾；預設 = 輸出夾")
    args = parser.parse_args()
    run(load_config(args.config), args.output_dir, args.input_dir)


if __name__ == "__main__":
    main()
