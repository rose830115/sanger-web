#!/usr/bin/env python3
"""Batch QC and trimming for Sanger AB1 files.

門檻來自 config.yaml 的 qc 區塊（原版寫死 250 / 50 / 85 / 25 / 40）。
用法：python3 -m pipeline.sanger_qc INPUT_DIR -o OUTPUT_DIR [--config config.yaml]
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from Bio import SeqIO

from pipeline.config import load_config
from pipeline.markers import UNKNOWN_MARKER, marker_for_primer


def parse_name(stem: str) -> dict[str, str]:
    parts = stem.split("_")
    primer = parts[2] if len(parts) > 2 else ""
    marker = marker_for_primer(primer)

    return {
        "run_id": parts[0] if parts else stem,
        "sample": parts[1] if len(parts) > 1 else "",
        "primer": primer,
        "marker": marker,
        "well": parts[3] if len(parts) > 3 else "",
    }


def status_for(
    trimmed_len: int, mean_q: float, q20_pct: float, qc: dict | None = None
) -> str:
    qc = qc or load_config()["qc"]
    if (
        trimmed_len < qc["min_trimmed_len"]
        or q20_pct < qc["fail_q20_pct"]
        or mean_q < qc["fail_mean_q"]
    ):
        return "FAIL_RESEQUENCE_OR_MANUAL_REVIEW"
    if q20_pct < qc["warn_q20_pct"] or mean_q < qc["warn_mean_q"]:
        return "WARN_MANUAL_REVIEW"
    return "PASS"


def write_pairing_report(rows: list[dict[str, str]], out_path: Path) -> None:
    groups: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in rows:
        groups.setdefault((row["sample"], row["marker"]), []).append(row)

    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sample", "marker", "read_count", "primers", "pairing_status"],
        )
        writer.writeheader()
        for (sample, marker), group in sorted(groups.items()):
            primers = ",".join(sorted({item["primer"] for item in group}))
            read_count = len(group)
            pairing_status = "POTENTIAL_PAIR_OR_MULTI_READ" if read_count >= 2 else "UNPAIRED"
            writer.writerow(
                {
                    "sample": sample,
                    "marker": marker,
                    "read_count": read_count,
                    "primers": primers,
                    "pairing_status": pairing_status,
                }
            )


def process(indir: Path, outdir: Path, qc: dict | None = None) -> None:
    qc = qc or load_config()["qc"]
    outdir.mkdir(parents=True, exist_ok=True)
    ab1_files = sorted(indir.glob("*.ab1"))
    if not ab1_files:
        raise SystemExit(f"No .ab1 files found in {indir}")

    rows: list[dict[str, str]] = []
    raw_records = []
    trimmed_records = []

    for ab1 in ab1_files:
        try:
            raw = SeqIO.read(str(ab1), "abi")
            trimmed = SeqIO.read(str(ab1), "abi-trim")
        except Exception:  # noqa: BLE001  # 🟡-13：壞檔只標這一條，不讓整個子批次 QC 崩潰
            meta = parse_name(ab1.stem)
            rows.append(
                {
                    "file": ab1.name,
                    **meta,
                    "raw_len": "0",
                    "trimmed_len": "0",
                    "removed_bp": "0",
                    "mean_q_raw": "0.0",
                    "q20_raw_pct": "0.0",
                    "status": "FAIL_READ_UNUSABLE",
                }
            )
            continue
        qualities = raw.letter_annotations.get("phred_quality", [])
        mean_q = sum(qualities) / len(qualities) if qualities else 0.0
        q20_pct = sum(q >= 20 for q in qualities) / len(qualities) * 100 if qualities else 0.0
        meta = parse_name(ab1.stem)
        status = status_for(len(trimmed.seq), mean_q, q20_pct, qc)
        record_id = (
            f"{meta['sample']}|{meta['marker']}|{meta['primer']}|{ab1.stem}"
        )

        raw.id = raw.name = record_id
        raw.description = "raw ABI basecall"
        trimmed.id = trimmed.name = record_id
        trimmed.description = f"trimmed ABI basecall status={status}"
        raw_records.append(raw)
        trimmed_records.append(trimmed)

        rows.append(
            {
                "file": ab1.name,
                **meta,
                "raw_len": str(len(raw.seq)),
                "trimmed_len": str(len(trimmed.seq)),
                "removed_bp": str(len(raw.seq) - len(trimmed.seq)),
                "mean_q_raw": f"{mean_q:.1f}",
                "q20_raw_pct": f"{q20_pct:.1f}",
                "status": status,
            }
        )

    if not rows:
        raise SystemExit(f"No readable .ab1 files in {indir}")
    SeqIO.write(raw_records, str(outdir / "raw_ab1_sequences.fasta"), "fasta")
    SeqIO.write(trimmed_records, str(outdir / "trimmed_sequences.fasta"), "fasta")

    with (outdir / "qc_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    with (outdir / "needs_review.txt").open("w") as handle:
        flagged = [row for row in rows if row["status"] != "PASS"]
        unknown = [row for row in rows if row["marker"] == UNKNOWN_MARKER]
        if not flagged and not unknown:
            handle.write("No sequences flagged for review.\n")
        else:
            if flagged:
                handle.write("Sequences flagged for manual chromatogram review or resequencing:\n")
            for row in flagged:
                handle.write(
                    f"- {row['file']}: {row['status']}; "
                    f"trimmed_len={row['trimmed_len']}; "
                    f"meanQ={row['mean_q_raw']}; Q20%={row['q20_raw_pct']}\n"
                )
            if unknown:
                handle.write("WARN unrecognized primer name (marker set to 其他, data kept):\n")
                for row in unknown:
                    handle.write(f"- {row['file']}: primer={row['primer']}\n")

    write_pairing_report(rows, outdir / "pairing_report.csv")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir", type=Path)
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        help="Defaults to INPUT_DIR_processed next to the input directory.",
    )
    parser.add_argument("--config", type=Path, default=None)
    args = parser.parse_args()
    cfg = load_config(args.config)
    outdir = args.output_dir or args.input_dir.with_name(f"{args.input_dir.name}_processed")
    process(args.input_dir, outdir, cfg["qc"])
    print(outdir)


if __name__ == "__main__":
    main()
