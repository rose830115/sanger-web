"""QC 門檻、scan_batches、跨批次正規化的基本測試。"""

import csv

from pipeline.config import load_config
from pipeline.cross_batch_summary import collect, normalize_species
from pipeline.sanger_qc import parse_name, status_for
from pipeline.scan_batches import order_id_of, scan

QC = load_config()["qc"]
KW = ["clone", "gene", "genes", "genomic"]


def test_qc_thresholds_from_config():
    assert status_for(300, 45, 90, QC) == "PASS"
    assert status_for(300, 35, 90, QC) == "WARN_MANUAL_REVIEW"
    assert status_for(300, 45, 70, QC) == "WARN_MANUAL_REVIEW"
    assert status_for(249, 45, 90, QC) == "FAIL_RESEQUENCE_OR_MANUAL_REVIEW"
    assert status_for(300, 45, 49, QC) == "FAIL_RESEQUENCE_OR_MANUAL_REVIEW"
    assert status_for(300, 24, 90, QC) == "FAIL_RESEQUENCE_OR_MANUAL_REVIEW"


def test_marker_from_primer_name():
    assert parse_name("S1_sampleA_ITS5_01")["marker"] == "ITS"
    assert parse_name("S1_sampleA_Bt2b_01")["marker"] == "Btub"
    assert parse_name("S1_sampleA_EF-1_01")["marker"] == "TEF1"
    assert parse_name("S1_sampleA_CF4_01")["marker"] == "CF"


def test_config_has_no_email_or_personal_path():
    from pipeline.config import DEFAULT_CONFIG_PATH

    assert load_config()["blast"]["email"] == ""
    text = DEFAULT_CONFIG_PATH.read_text(encoding="utf-8")
    personal_markers = ["mac" + "-4", "g" + "mail", "@", "/Users/"]
    assert not [m for m in personal_markers if m in text]


def test_scan_batches(tmp_path):
    raw, out = tmp_path / "raw", tmp_path / "out"
    for name in ["S1_A", "S1_B", "S2_C", "S3_D"]:
        (raw / name).mkdir(parents=True)
    out.mkdir()
    # S1：QC 都做了、Assembly 做了、BLAST 完成 -> 不待處理
    (out / "S1_A_processed").mkdir()
    (out / "S1_B_processed").mkdir()
    (out / "S1_combined_assembly").mkdir()
    (out / "S1_ncbi_blast").mkdir()
    (out / "S1_ncbi_blast" / "taxonomy_decision_table.csv").write_text("x")
    # S2：只做 QC
    (out / "S2_C_processed").mkdir()
    res = scan(raw, out)
    assert order_id_of("S1_A") == "S1"
    assert res["pending_orders"] == ["S2", "S3"]
    assert res["orders"]["S2"]["need_assembly"] and res["orders"]["S2"]["need_blast"]
    assert res["orders"]["S2"]["need_qc"] == []
    assert res["orders"]["S3"]["need_qc"] == ["S3_D"]
    assert res["pending_subfolders"] == ["S2_C", "S3_D"]


def test_normalize_species():
    assert normalize_species("Fusarium clone 12", KW) == "Fusarium sp."
    assert normalize_species("Aspergillus genomic", KW) == "Aspergillus sp."
    assert normalize_species("Aspergillus sp", KW) == "Aspergillus sp."
    assert normalize_species("Aspergillus flavus", KW) == "Aspergillus flavus"
    assert normalize_species("", KW) is None


def test_cross_batch_collect(tmp_path):
    d = tmp_path / "S1_ncbi_blast"
    d.mkdir()
    with (d / "blast_top100_hits.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["sample", "organism_guess", "mismatches"])
        w.writeheader()
        w.writerows([
            {"sample": "A", "organism_guess": "Aspergillus flavus", "mismatches": "2"},
            {"sample": "A", "organism_guess": "Aspergillus flavus", "mismatches": "3"},   # 重複
            {"sample": "A", "organism_guess": "Aspergillus oryzae", "mismatches": "11"},  # 超過門檻
            {"sample": "A", "organism_guess": "Penicillium clone", "mismatches": "1"},
        ])
    assert collect(tmp_path, KW, 10) == {"A": ["Aspergillus flavus", "Penicillium sp."]}


def test_summarize_its_in_non_sufficient_genus_not_candidate_ok(tmp_path):
    """🔴-3：Penicillium 只有 ITS 全同種名，不可輸出 CANDIDATE_OK。"""
    from pipeline.consensus_call import load_lookup
    from pipeline.config import gene_lookup_path
    from pipeline.summarize import summarize_orders

    d = tmp_path / "S1_ncbi_blast"
    d.mkdir()
    cols = ["sample", "marker", "candidate_species", "top_hit_organism", "top_pident",
            "top_qcov_pct", "review_flag", "suggested_action", "candidate_genus"]
    with (d / "taxonomy_decision_table.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for genus, sp in [("Penicillium", "Penicillium chrysogenum"), ("Lichtheimia", "Lichtheimia corymbifera")]:
            w.writerow({"sample": genus, "marker": "ITS", "candidate_species": sp, "top_hit_organism": sp,
                        "top_pident": "100", "top_qcov_pct": "100", "review_flag": "OK",
                        "suggested_action": "CANDIDATE_OK_IF_SEQUENCE_QC_ACCEPTED", "candidate_genus": genus})
    lookup = load_lookup(gene_lookup_path(load_config()))
    got = {r["樣本"]: r["建議動作"] for r in summarize_orders(tmp_path, ["S1"], lookup)}
    assert got["Penicillium"] == "ITS_GENUS_LEVEL_ONLY"
    assert got["Lichtheimia"] == "CANDIDATE_OK_IF_SEQUENCE_QC_ACCEPTED"


# ---------------------------------------------------------------- 🔴-7 覆蓋率分母與最短對齊長度
from types import SimpleNamespace as NS  # noqa: E402

from pipeline.sanger_ncbi_blast import hit_is_high_quality, hit_row, summarize_query  # noqa: E402


def _fake_hit(query_len, hit_len, span, ident=None):
    hsp = NS(query_start=1, query_end=span, identities=ident or span, align_length=span, gaps=0,
             expect=0.0, bits=1000.0, sbjct_start=1, sbjct_end=span, frame=None)
    aln = NS(hsps=[hsp], hit_def="X1 Penicillium copticola strain Y", accession="X1", hit_id="gb|X1|", length=hit_len)
    rec = NS(query="S|CF|p", query_length=query_len)
    return hit_row(rec, aln, 1, "CF", "nt", "2026-01-01")


def test_qcov_uses_shorter_of_query_and_reference():
    # 共識 770 bp，參考只有 690 bp，全長吻合：舊算法 88.4%（被排除），新算法 100%
    row = _fake_hit(770, 690, 680)
    assert float(row["qcov_query_pct"]) < 95
    assert float(row["qcov_pct"]) >= 95
    assert hit_is_high_quality(row, 99.0, 95.0, 300)


def test_qcov_query_shorter_than_reference_unchanged():
    row = _fake_hit(500, 2000, 400)         # 長參考：分母仍是 query
    assert abs(float(row["qcov_pct"]) - 80.0) < 0.01


def test_min_align_len_blocks_short_conserved_hits():
    row = _fake_hit(700, 150, 150)           # 短參考完整對上 150 bp：覆蓋率 100% 但太短
    assert float(row["qcov_pct"]) == 100.0
    assert not hit_is_high_quality(row, 99.0, 95.0, 300)
    assert hit_is_high_quality(row, 99.0, 95.0, 0)


def test_min_align_len_not_applied_when_query_itself_short():
    row = _fake_hit(280, 600, 270)
    assert hit_is_high_quality(row, 99.0, 95.0, 300)


def test_summarize_query_applies_min_align_len():
    short = _fake_hit(700, 150, 150)
    longer = _fake_hit(700, 700, 690)
    longer["species_guess"] = "Penicillium rubens"
    out = summarize_query("S|CF|p", [{**short, "rank": "1"}, {**longer, "rank": "2"}], 99.0, 95.0, 300)
    assert "Penicillium rubens" in out["high_quality_species_all"]
    assert "copticola" not in out["high_quality_species_all"]


def test_config_has_min_align_len():
    assert load_config()["blast"]["min_align_len"] == 300


def test_summarize_query_writes_species_best_pident_sorted():
    a = _fake_hit(700, 700, 690, ident=689)       # 99.855%
    a["species_guess"] = "Penicillium rubens"
    b = _fake_hit(700, 700, 690, ident=683)       # 98.986%
    b["species_guess"] = "Penicillium chrysogenum"
    c = _fake_hit(700, 700, 690, ident=680)       # 較低的同種 hit 不覆蓋較高分
    c["species_guess"] = "Penicillium rubens"
    out = summarize_query("S|CF|p", [{**a, "rank": "1"}, {**b, "rank": "2"}, {**c, "rank": "3"}], 99.0, 95.0, 300)
    parts = out["species_best_pident_all"].split("; ")
    assert parts[0].startswith("Penicillium rubens:99.8") and parts[1].startswith("Penicillium chrysogenum:98.9")
