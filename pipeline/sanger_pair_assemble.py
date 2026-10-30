#!/usr/bin/env python3
"""Pair forward/reverse Sanger AB1 reads and build simple consensus sequences.

門檻來自 config.yaml 的 qc / assembly 區塊。
用法：python3 -m pipeline.sanger_pair_assemble DIR1 [DIR2 ...] -o OUTPUT_DIR [--config config.yaml]
"""

from __future__ import annotations

import argparse
import csv
import warnings
from pathlib import Path

from Bio import BiopythonDeprecationWarning
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pipeline.config import load_config
from pipeline.markers import is_forward
from pipeline.sanger_qc import parse_name, status_for


warnings.simplefilter("ignore", BiopythonDeprecationWarning)
from Bio import pairwise2  # noqa: E402


def read_trimmed_ab1(path: Path, qc: dict | None = None) -> dict[str, object]:
    raw = SeqIO.read(str(path), "abi")
    trimmed = SeqIO.read(str(path), "abi-trim")
    qualities = raw.letter_annotations.get("phred_quality", [])
    mean_q = sum(qualities) / len(qualities) if qualities else 0.0
    q20_pct = sum(q >= 20 for q in qualities) / len(qualities) * 100 if qualities else 0.0
    meta = parse_name(path.stem)
    status = status_for(len(trimmed.seq), mean_q, q20_pct, qc)
    return {
        "path": path,
        "file": path.name,
        **meta,
        "seq": str(trimmed.seq).upper(),
        "qual": trimmed.letter_annotations.get("phred_quality", []),
        "trimmed_len": len(trimmed.seq),
        "mean_q": mean_q,
        "q20_pct": q20_pct,
        "qc_status": status,
    }


def reverse_read(read: dict[str, object]) -> dict[str, object]:
    copy = dict(read)
    copy["seq"] = str(Seq(str(read["seq"])).reverse_complement())
    copy["qual"] = list(reversed(read["qual"]))  # type: ignore[arg-type]
    return copy


def align_pair(seq_a: str, seq_b: str):
    return pairwise2.align.globalms(
        seq_a,
        seq_b,
        2,
        -1,
        -5,
        -0.5,
        one_alignment_only=True,
    )[0]


def alignment_stats(alignment) -> dict[str, float]:
    aligned_a, aligned_b = alignment.seqA, alignment.seqB
    overlap = sum(a != "-" and b != "-" for a, b in zip(aligned_a, aligned_b))
    matches = sum(a == b and a != "-" for a, b in zip(aligned_a, aligned_b))
    mismatches = sum(
        a != "-"
        and b != "-"
        and a != b
        for a, b in zip(aligned_a, aligned_b)
    )
    identity = matches / overlap * 100 if overlap else 0.0
    return {
        "score": float(alignment.score),
        "overlap_bp": float(overlap),
        "matches": float(matches),
        "mismatches": float(mismatches),
        "identity_pct": identity,
    }


def choose_orientation(read_a: dict[str, object], read_b: dict[str, object]):
    seq_a = str(read_a["seq"])
    direct_alignment = align_pair(seq_a, str(read_b["seq"]))
    reverse = reverse_read(read_b)
    reverse_alignment = align_pair(seq_a, str(reverse["seq"]))
    direct_stats = alignment_stats(direct_alignment)
    reverse_stats = alignment_stats(reverse_alignment)
    if reverse_stats["score"] >= direct_stats["score"]:
        return reverse, reverse_alignment, "reverse_complement", reverse_stats
    return read_b, direct_alignment, "direct", direct_stats


def base_quality(qualities: list[int], index: int) -> int:
    if 0 <= index < len(qualities):
        return int(qualities[index])
    return 0


def trim_n_edges(seq: str) -> str:
    return seq.strip("N")


def consensus_quality(
    base_a: str,
    q_a: int,
    base_b: str,
    q_b: int,
    consensus_base: str,
    decision: str,
) -> int:
    if consensus_base == "N":
        return 0
    if decision == "two_read_agreement":
        return min(60, q_a + q_b)
    if decision == "single_read_a":
        return q_a
    if decision == "single_read_b":
        return q_b
    if decision == "conflict_use_read_a":
        return max(0, q_a - q_b)
    if decision == "conflict_use_read_b":
        return max(0, q_b - q_a)
    return 0


def trim_consensus_rows(rows: list[dict[str, object]]) -> tuple[str, list[dict[str, object]]]:
    start = 0
    end = len(rows)
    while start < end and rows[start]["consensus_base"] == "N":
        start += 1
    while end > start and rows[end - 1]["consensus_base"] == "N":
        end -= 1

    trimmed_rows = rows[start:end]
    for index, row in enumerate(trimmed_rows, start=1):
        row["consensus_pos"] = index
    return "".join(str(row["consensus_base"]) for row in trimmed_rows), trimmed_rows


def consensus_from_alignment(
    sample: str,
    marker: str,
    read_a: dict[str, object],
    read_b: dict[str, object],
    alignment,
    asm: dict | None = None,
) -> tuple[str, int, int, list[dict[str, object]]]:
    asm = asm or load_config()["assembly"]
    qa = list(read_a["qual"])  # type: ignore[arg-type]
    qb = list(read_b["qual"])  # type: ignore[arg-type]
    pos_a = 0
    pos_b = 0
    base_rows: list[dict[str, object]] = []
    ambiguous = 0
    conflicts = 0

    for alignment_pos, (base_a, base_b) in enumerate(zip(alignment.seqA, alignment.seqB), start=1):
        read_a_pos = pos_a + 1 if base_a != "-" else ""
        read_b_pos = pos_b + 1 if base_b != "-" else ""
        q_a = base_quality(qa, pos_a) if base_a != "-" else 0
        q_b = base_quality(qb, pos_b) if base_b != "-" else 0
        decision = "gap"

        if base_a != "-" and base_b != "-":
            if base_a == base_b:
                consensus_base = base_a
                decision = "two_read_agreement"
            else:
                conflicts += 1
                if abs(q_a - q_b) >= asm["conflict_q_diff"]:
                    consensus_base = base_a if q_a > q_b else base_b
                    decision = "conflict_use_read_a" if q_a > q_b else "conflict_use_read_b"
                else:
                    consensus_base = "N"
                    decision = "conflict_to_N"
                    ambiguous += 1
        elif base_a != "-":
            if q_a >= asm["single_read_min_q"]:
                consensus_base = base_a
                decision = "single_read_a"
            else:
                consensus_base = "N"
                decision = "single_read_a_low_q_to_N"
                ambiguous += 1
        elif base_b != "-":
            if q_b >= asm["single_read_min_q"]:
                consensus_base = base_b
                decision = "single_read_b"
            else:
                consensus_base = "N"
                decision = "single_read_b_low_q_to_N"
                ambiguous += 1
        else:
            consensus_base = "N"

        consensus_q = consensus_quality(base_a, q_a, base_b, q_b, consensus_base, decision)
        needs_review = (
            consensus_base == "N"
            or consensus_q < 20
            or decision.startswith("conflict")
        )
        base_rows.append(
            {
                "sample": sample,
                "marker": marker,
                "alignment_pos": alignment_pos,
                "consensus_pos": "",
                "consensus_base": consensus_base,
                "consensus_q_est": consensus_q,
                "read_a_primer": read_a["primer"],
                "read_a_pos": read_a_pos,
                "read_a_base": base_a,
                "read_a_q": q_a,
                "read_b_primer": read_b["primer"],
                "read_b_pos": read_b_pos,
                "read_b_base": base_b,
                "read_b_q": q_b,
                "decision": decision,
                "needs_review": "yes" if needs_review else "no",
            }
        )

        if base_a != "-":
            pos_a += 1
        if base_b != "-":
            pos_b += 1

    consensus, trimmed_rows = trim_consensus_rows(base_rows)
    return consensus, conflicts, ambiguous, trimmed_rows


def write_alignment_checkpoint(
    out_path: Path,
    read_a: dict[str, object],
    read_b: dict[str, object],
    alignment,
    base_rows: list[dict[str, object]],
    stats: dict[str, float],
    status: str,
    width: int = 80,
) -> None:
    by_alignment_pos = {int(row["alignment_pos"]): row for row in base_rows}
    consensus_chars = []
    review_chars = []
    for index, _ in enumerate(alignment.seqA, start=1):
        row = by_alignment_pos.get(index)
        consensus_chars.append(str(row["consensus_base"]) if row else " ")
        review_chars.append("!" if row and row["needs_review"] == "yes" else " ")

    match_line = "".join(
        "|" if a == b and a != "-" else " "
        for a, b in zip(alignment.seqA, alignment.seqB)
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as handle:
        handle.write(f"status: {status}\n")
        handle.write(f"read_a: {read_a['file']} ({read_a['primer']})\n")
        handle.write(f"read_b: {read_b['file']} ({read_b['primer']})\n")
        handle.write(
            f"identity: {stats['identity_pct']:.2f}% | overlap: {int(stats['overlap_bp'])} bp | "
            f"mismatches: {int(stats['mismatches'])}\n"
        )
        handle.write("!: conflict, N, or estimated consensus Q < 20\n\n")
        for start in range(0, len(alignment.seqA), width):
            end = start + width
            handle.write(f"read_a     {alignment.seqA[start:end]}\n")
            handle.write(f"match      {match_line[start:end]}\n")
            handle.write(f"read_b     {alignment.seqB[start:end]}\n")
            handle.write(f"consensus  {''.join(consensus_chars[start:end])}\n")
            handle.write(f"review     {''.join(review_chars[start:end])}\n\n")


def assembly_status(
    reads: list[dict[str, object]],
    stats: dict[str, float],
    ambiguous: int,
    asm: dict | None = None,
) -> str:
    asm = asm or load_config()["assembly"]
    qc_statuses = {str(read["qc_status"]) for read in reads}
    if any(status.startswith("FAIL") for status in qc_statuses):
        return "FAIL_RESEQUENCE_OR_MANUAL_REVIEW"
    if stats["overlap_bp"] < asm["min_overlap_bp"] or stats["identity_pct"] < asm["fail_identity_pct"]:
        return "FAIL_ASSEMBLY_CHECK"
    if any(status.startswith("WARN") for status in qc_statuses):
        return "WARN_MANUAL_REVIEW"
    if stats["identity_pct"] < asm["warn_identity_pct"] or ambiguous > asm["max_ambiguous_n"]:
        return "WARN_CONSENSUS_CHECK"
    return "PASS"


def primer_order(primer: str) -> tuple[int, str]:
    return (0 if is_forward(primer) else 1, primer)


ROW_FIELDS = [
    "sample", "marker", "read_count", "files", "orientation", "overlap_bp", "identity_pct", "mismatches",
    "conflicts", "ambiguous_N", "low_confidence_sites", "consensus_len", "read_qc_statuses", "assembly_status",
]


def _fail_row(sample: str, marker: str, reads: list[dict[str, object]], status: str) -> dict[str, str]:
    row = {k: "" for k in ROW_FIELDS}
    row.update(
        {
            "sample": sample,
            "marker": marker,
            "read_count": str(len(reads)),
            "files": ";".join(str(read["file"]) for read in reads),
            "read_qc_statuses": ";".join(str(read["qc_status"]) for read in reads),
            "assembly_status": status,
        }
    )
    return row


def _assemble_pair(sample: str, marker: str, reads: list[dict[str, object]], outdir: Path, asm: dict) -> dict:
    reads = sorted(reads, key=lambda read: primer_order(str(read["primer"])))
    read_a, read_b = reads
    oriented_b, alignment, orientation, stats = choose_orientation(read_a, read_b)
    consensus, conflicts, ambiguous, base_rows = consensus_from_alignment(
        sample, marker, read_a, oriented_b, alignment, asm,
    )
    pair_low_confidence = [row for row in base_rows if row["needs_review"] == "yes"]
    status = assembly_status([read_a, read_b], stats, ambiguous, asm)
    primers = f"{read_a['primer']}+{read_b['primer']}"
    record_id = f"{sample}|{marker}|{primers}"
    safe_id = f"{sample}_{marker}_{read_a['primer']}_{read_b['primer']}".replace("/", "-")
    write_alignment_checkpoint(
        outdir / "alignment_checkpoints" / f"{safe_id}.txt",
        read_a, oriented_b, alignment, base_rows, stats, status,
    )
    record = SeqRecord(
        Seq(consensus),
        id=record_id,
        description=(
            f"status={status} orientation={orientation} "
            f"identity={stats['identity_pct']:.2f}% overlap={int(stats['overlap_bp'])}"
        ),
    )
    row = {
        "sample": sample,
        "marker": marker,
        "read_count": "2",
        "files": f"{read_a['file']};{read_b['file']}",
        "orientation": orientation,
        "overlap_bp": str(int(stats["overlap_bp"])),
        "identity_pct": f"{stats['identity_pct']:.2f}",
        "mismatches": str(int(stats["mismatches"])),
        "conflicts": str(conflicts),
        "ambiguous_N": str(ambiguous),
        "low_confidence_sites": str(len(pair_low_confidence)),
        "consensus_len": str(len(consensus)),
        "read_qc_statuses": f"{read_a['qc_status']};{read_b['qc_status']}",
        "assembly_status": status,
    }
    review_line = ""
    if status != "PASS":
        review_line = (
            f"- {sample} {marker}: {status}; identity={stats['identity_pct']:.2f}%; "
            f"overlap={int(stats['overlap_bp'])} bp; files={read_a['file']}, {read_b['file']}"
        )
    return {
        "base_rows": base_rows, "low_confidence": pair_low_confidence, "record": record,
        "status": status, "row": row, "review_line": review_line,
    }


def process(
    input_dirs: list[Path],
    outdir: Path,
    qc: dict | None = None,
    asm: dict | None = None,
) -> None:
    cfg = None if (qc and asm) else load_config()
    qc = qc or cfg["qc"]
    asm = asm or cfg["assembly"]
    outdir.mkdir(parents=True, exist_ok=True)
    reads_by_pair: dict[tuple[str, str], list[dict[str, object]]] = {}
    bad_by_pair: dict[tuple[str, str], list[str]] = {}   # 🟡-13：讀不了的 ab1 不炸整單
    for input_dir in input_dirs:
        for path in sorted(input_dir.glob("*.ab1")):
            try:
                read = read_trimmed_ab1(path, qc)
            except Exception as exc:  # noqa: BLE001
                meta = parse_name(path.stem)
                bad_by_pair.setdefault((meta["sample"], meta["marker"]), []).append(f"{path.name}（{exc}）")
                continue
            reads_by_pair.setdefault((str(read["sample"]), str(read["marker"])), []).append(read)

    rows: list[dict[str, str]] = []
    records: list[SeqRecord] = []
    blast_records: list[SeqRecord] = []
    pass_records: list[SeqRecord] = []
    review_lines: list[str] = []
    all_base_rows: list[dict[str, object]] = []
    low_confidence_rows: list[dict[str, object]] = []

    for key in bad_by_pair:
        reads_by_pair.setdefault(key, [])

    for (sample, marker), reads in sorted(reads_by_pair.items()):
        bad = bad_by_pair.get((sample, marker), [])
        if bad:
            rows.append(_fail_row(sample, marker, reads, "FAIL_READ_UNUSABLE"))
            review_lines.append(f"- {sample} {marker}: FAIL_READ_UNUSABLE; 讀不了的檔案：{'; '.join(bad)}")
            continue
        if len(reads) != 2:
            rows.append(_fail_row(sample, marker, reads, "FAIL_PAIRING_CHECK"))
            review_lines.append(f"- {sample} {marker}: expected 2 reads, found {len(reads)}")
            continue

        try:
            pair = _assemble_pair(sample, marker, reads, outdir, asm)
        except Exception as exc:  # noqa: BLE001  # 🟡-13：一對壞 read（例如修剪後 0 bp）只讓這一對失敗
            rows.append(_fail_row(sample, marker, reads, "FAIL_READ_UNUSABLE"))
            review_lines.append(f"- {sample} {marker}: FAIL_READ_UNUSABLE; {type(exc).__name__}: {exc}")
            continue
        all_base_rows.extend(pair["base_rows"])
        low_confidence_rows.extend(pair["low_confidence"])
        records.append(pair["record"])
        if not pair["status"].startswith("FAIL"):
            blast_records.append(pair["record"])
        if pair["status"] == "PASS":
            pass_records.append(pair["record"])
        rows.append(pair["row"])
        if pair["review_line"]:
            review_lines.append(pair["review_line"])

    SeqIO.write(records, str(outdir / "paired_consensus.fasta"), "fasta")
    SeqIO.write(blast_records, str(outdir / "paired_consensus_for_blast.fasta"), "fasta")
    SeqIO.write(pass_records, str(outdir / "paired_consensus_pass_only.fasta"), "fasta")
    with (outdir / "assembly_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ROW_FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    base_fieldnames = [
        "sample",
        "marker",
        "alignment_pos",
        "consensus_pos",
        "consensus_base",
        "consensus_q_est",
        "read_a_primer",
        "read_a_pos",
        "read_a_base",
        "read_a_q",
        "read_b_primer",
        "read_b_pos",
        "read_b_base",
        "read_b_q",
        "decision",
        "needs_review",
    ]
    with (outdir / "consensus_base_support.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=base_fieldnames)
        writer.writeheader()
        writer.writerows(all_base_rows)

    with (outdir / "low_confidence_sites.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=base_fieldnames)
        writer.writeheader()
        writer.writerows(low_confidence_rows)

    with (outdir / "assembly_needs_review.txt").open("w") as handle:
        if review_lines:
            handle.write("Assemblies flagged for review:\n")
            handle.write("\n".join(review_lines))
            handle.write("\n")
        else:
            handle.write("No assemblies flagged for review.\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dirs", type=Path, nargs="+")
    parser.add_argument("-o", "--output-dir", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=None)
    args = parser.parse_args()
    cfg = load_config(args.config)
    process(args.input_dirs, args.output_dir, cfg["qc"], cfg["assembly"])
    print(args.output_dir)


if __name__ == "__main__":
    main()
