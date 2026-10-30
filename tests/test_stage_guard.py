"""🔴-6／🟡-13：輸入雜湊、階段失效、XML 逐條驗證、NOT_BLASTED、壞 read 不炸整單、暫存改名。"""

import csv
from pathlib import Path

import pytest
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import run_all
from pipeline import consensus_call as cc
from pipeline import sanger_ncbi_blast as nb
from pipeline import sanger_pair_assemble as spa
from pipeline import stage_guard as sg
from pipeline.config import gene_lookup_path, load_config
from pipeline.scan_batches import scan

CFG = load_config()


# ------------------------------------------------------------ 測試用假資料
def make_raw(raw: Path, subs: dict[str, dict[str, bytes]]):
    for sub, files in subs.items():
        (raw / sub).mkdir(parents=True, exist_ok=True)
        for name, data in files.items():
            (raw / sub / name).write_bytes(data)


def xml_for(queries: dict[str, int]) -> str:
    """最小可解析的 NCBI BLAST XML：每條查詢一個 Iteration，附一筆 100% 的 hit。"""
    its = ""
    for i, (qid, qlen) in enumerate(queries.items(), start=1):
        its += f"""
<Iteration><Iteration_iter-num>{i}</Iteration_iter-num><Iteration_query-ID>Query_{i}</Iteration_query-ID>
<Iteration_query-def>{qid}</Iteration_query-def><Iteration_query-len>{qlen}</Iteration_query-len>
<Iteration_hits><Hit><Hit_num>1</Hit_num><Hit_id>gb|AB1|</Hit_id><Hit_def>Aspergillus niger strain X</Hit_def>
<Hit_accession>AB1</Hit_accession><Hit_len>{qlen}</Hit_len><Hit_hsps><Hsp><Hsp_num>1</Hsp_num>
<Hsp_bit-score>1000</Hsp_bit-score><Hsp_score>500</Hsp_score><Hsp_evalue>0</Hsp_evalue>
<Hsp_query-from>1</Hsp_query-from><Hsp_query-to>{qlen}</Hsp_query-to><Hsp_hit-from>1</Hsp_hit-from><Hsp_hit-to>{qlen}</Hsp_hit-to>
<Hsp_query-frame>1</Hsp_query-frame><Hsp_hit-frame>1</Hsp_hit-frame><Hsp_identity>{qlen}</Hsp_identity>
<Hsp_positive>{qlen}</Hsp_positive><Hsp_gaps>0</Hsp_gaps><Hsp_align-len>{qlen}</Hsp_align-len>
<Hsp_qseq>{'A' * 5}</Hsp_qseq><Hsp_hseq>{'A' * 5}</Hsp_hseq><Hsp_midline>|||||</Hsp_midline></Hsp></Hit_hsps></Hit></Iteration_hits>
<Iteration_stat><Statistics><Statistics_db-num>1</Statistics_db-num><Statistics_db-len>1</Statistics_db-len>
<Statistics_hsp-len>0</Statistics_hsp-len><Statistics_eff-space>0</Statistics_eff-space><Statistics_kappa>0.1</Statistics_kappa>
<Statistics_lambda>1.0</Statistics_lambda><Statistics_entropy>0.1</Statistics_entropy></Statistics></Iteration_stat></Iteration>"""
    return f"""<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput><BlastOutput_program>blastn</BlastOutput_program><BlastOutput_version>BLASTN 2.16.0+</BlastOutput_version>
<BlastOutput_reference>x</BlastOutput_reference><BlastOutput_db>nt</BlastOutput_db><BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
<BlastOutput_query-def>x</BlastOutput_query-def><BlastOutput_query-len>1</BlastOutput_query-len>
<BlastOutput_param><Parameters><Parameters_expect>10</Parameters_expect></Parameters></BlastOutput_param>
<BlastOutput_iterations>{its}
</BlastOutput_iterations></BlastOutput>"""


def recs(spec: dict[str, str]):
    return [SeqRecord(Seq(seq), id=qid, description="") for qid, seq in spec.items()]


def write_fasta(path: Path, spec: dict[str, str]):
    path.write_text("".join(f">{k}\n{v}\n" for k, v in spec.items()))


SEQ = "ACGT" * 100          # 400 bp


# ------------------------------------------------------------ 雜湊與失效判定
def test_judge_states():
    assert sg.judge(False, None, "a") == "missing"
    assert sg.judge(True, None, "a") == "legacy"
    assert sg.judge(True, "a", "a") == "fresh"
    assert sg.judge(True, "a", "b") == "stale"


def test_scan_detects_changed_ab1_and_late_subbatch(tmp_path):
    raw, out = tmp_path / "raw", tmp_path / "out"
    make_raw(raw, {"S1_A": {"S1x_s1_ITS4_01.ab1": b"one", "S1x_s1_ITS5_02.ab1": b"two"}})
    info = scan(raw, out, CFG)["orders"]["S1"]
    assert info["need_qc"] == ["S1_A"]

    # 模擬完整跑完並記錄指紋
    (out / "S1_A_processed").mkdir(parents=True)
    (out / "S1_combined_assembly").mkdir()
    (out / "S1_ncbi_blast").mkdir()
    (out / "S1_ncbi_blast" / "taxonomy_decision_table.csv").write_text("x")
    rep = sg.stage_report(raw, out, "S1", info, CFG)
    sg.save_manifest(out, "S1", {"qc": {"S1_A": {"fp": rep["current"]["qc"]["S1_A"]}},
                                 "assembly": {"fp": rep["current"]["assembly"]},
                                 "blast": {"fp": "whatever"}})
    res = scan(raw, out, CFG)
    assert res["pending_orders"] == [] or res["orders"]["S1"]["stale_blast"]   # blast 指紋用假值 -> 只有 blast 失效
    assert not res["orders"]["S1"]["stale_qc"] and not res["orders"]["S1"]["stale_assembly"]

    # 1) 孔位標錯，內容被換掉（檔名不變、內容變）-> QC 與 Assembly 都失效
    (raw / "S1_A" / "S1x_s1_ITS4_01.ab1").write_bytes(b"changed")
    res = scan(raw, out, CFG)
    assert res["orders"]["S1"]["stale_qc"] == ["S1_A"]
    assert res["orders"]["S1"]["stale_assembly"] is True
    assert "S1" in res["pending_orders"]


def test_scan_late_second_subbatch_marks_assembly_stale(tmp_path):
    """③ 同訂單第二子批次晚到：QC 做了、Assembly／BLAST 以前會被跳過，現在要重做。"""
    raw, out = tmp_path / "raw", tmp_path / "out"
    make_raw(raw, {"S1_A": {"S1x_s1_ITS4_01.ab1": b"one"}})
    (out / "S1_A_processed").mkdir(parents=True)
    (out / "S1_combined_assembly").mkdir()
    (out / "S1_ncbi_blast").mkdir()
    (out / "S1_ncbi_blast" / "taxonomy_decision_table.csv").write_text("x")
    assert scan(raw, out)["pending_orders"] == []          # 舊行為（不給 cfg）看不出來
    make_raw(raw, {"S1_B": {"S1y_s2_ITS4_01.ab1": b"two"}})
    res = scan(raw, out, CFG)
    assert res["pending_orders"] == ["S1"]
    assert res["orders"]["S1"]["need_qc"] == ["S1_B"]
    assert res["orders"]["S1"]["stale_assembly"] is True


def test_blast_fingerprint_changes_with_fasta_and_config(tmp_path):
    d = tmp_path / "asm"
    d.mkdir()
    assert sg.blast_fingerprint(d, CFG) is None
    write_fasta(d / "paired_consensus_for_blast.fasta", {"a|ITS|p": SEQ})
    f1 = sg.blast_fingerprint(d, CFG)
    write_fasta(d / "paired_consensus_for_blast.fasta", {"a|ITS|p": SEQ[::-1]})
    f2 = sg.blast_fingerprint(d, CFG)
    cfg2 = {**CFG, "blast": {**CFG["blast"], "min_align_len": 123}}
    f3 = sg.blast_fingerprint(d, cfg2)
    assert len({f1, f2, f3}) == 3


def test_atomic_dir_renames_only_on_success(tmp_path):
    target = tmp_path / "S1_combined_assembly"
    with pytest.raises(RuntimeError):
        with sg.atomic_dir(target) as tmp:
            (tmp / "partial.txt").write_text("x")
            raise RuntimeError("boom")
    assert not target.exists()
    assert not [p for p in tmp_path.iterdir()]          # 暫存夾也清掉，不留空夾
    with sg.atomic_dir(target) as tmp:
        (tmp / "ok.txt").write_text("x")
    assert (target / "ok.txt").is_file()


def test_retire_to_stale_moves_not_deletes(tmp_path):
    d = tmp_path / "S1_combined_assembly"
    d.mkdir()
    (d / "f.txt").write_text("keep")
    dest = sg.retire_to_stale(tmp_path, d)
    assert not d.exists() and (dest / "f.txt").read_text() == "keep"


# ------------------------------------------------------------ run_all：上游變 -> 下游失效（重現報告的錯配）
@pytest.fixture
def fake_stages(monkeypatch):
    calls = {"qc": [], "assembly": [], "blast": []}

    def fake_qc(indir, outdir, qc=None):
        calls["qc"].append(indir.name)
        (outdir / "qc_summary.csv").write_text("x")

    def fake_asm(dirs, outdir, qc=None, asm=None):
        calls["assembly"].append(sorted(d.name for d in dirs))
        write_fasta(outdir / "paired_consensus_for_blast.fasta", {"s1|ITS|p": SEQ + "A" * sum(len(f.read_bytes()) for d in dirs for f in d.glob("*.ab1"))})
        (outdir / "assembly_summary.csv").write_text("sample,marker,assembly_status\ns1,ITS,PASS\n")

    def fake_blast(fasta, outdir, blast_cfg, assembly_summary=None, by_marker=True, offline=False, on_event=None):
        calls["blast"].append(fasta.parent.name)
        outdir.mkdir(parents=True, exist_ok=True)
        (outdir / "taxonomy_decision_table.csv").write_text("sample,marker\n")

    monkeypatch.setattr(run_all.sanger_qc, "process", fake_qc)
    monkeypatch.setattr(run_all.sanger_pair_assemble, "process", fake_asm)
    monkeypatch.setattr(run_all, "run_blast", fake_blast)
    return calls


def test_run_all_skips_fresh_and_redoes_downstream_on_change(tmp_path, fake_stages):
    raw, out = tmp_path / "raw", tmp_path / "out"
    make_raw(raw, {"S1_A": {"S1x_s1_ITS4_01.ab1": b"one", "S1x_s1_ITS5_02.ab1": b"two"}})
    run = lambda **kw: run_all.run_pipeline(CFG, "qc", "blast", ["S1"], raw, out, offline=True, **kw)  # noqa: E731

    run()
    assert fake_stages == {"qc": ["S1_A"], "assembly": [["S1_A"]], "blast": ["S1_combined_assembly"]}
    assert sg.manifest_path(out, "S1").is_file()
    assert "files" in sg.load_manifest(out, "S1")           # 訂單輸入清單有逐檔雜湊

    run()                                                    # 沒變 -> 全部沿用
    assert len(fake_stages["qc"]) == 1 and len(fake_stages["assembly"]) == 1 and len(fake_stages["blast"]) == 1

    (raw / "S1_A" / "S1x_s1_ITS4_01.ab1").write_bytes(b"renamed-or-replaced")   # 上游變
    run()
    assert len(fake_stages["qc"]) == 2 and len(fake_stages["assembly"]) == 2 and len(fake_stages["blast"]) == 2
    assert list((out / "_stale").glob("S1_combined_assembly_*"))      # 舊版保留不刪

    make_raw(raw, {"S1_B": {"S1y_s1_ITS4_01.ab1": b"late"}})          # 第二子批次晚到
    run()
    assert fake_stages["assembly"][-1] == ["S1_A", "S1_B"]
    assert len(fake_stages["blast"]) == 3


def test_run_all_force_does_not_trust_old_blast_when_fasta_changes(tmp_path, fake_stages):
    raw, out = tmp_path / "raw", tmp_path / "out"
    make_raw(raw, {"S1_A": {"S1x_s1_ITS4_01.ab1": b"one"}})
    run_all.run_pipeline(CFG, "qc", "blast", ["S1"], raw, out, offline=True)
    make_raw(raw, {"S1_B": {"S1y_s9_ITS4_01.ab1": b"late"}})
    run_all.run_pipeline(CFG, "qc", "blast", ["S1"], raw, out, offline=True, force=True)
    assert len(fake_stages["blast"]) == 2                   # --force 之後 BLAST 解析表也跟著重做


# ------------------------------------------------------------ BLAST XML 逐條驗證、NOT_BLASTED
def blast_cfg(**kw):
    c = dict(CFG["blast"])
    c.update(kw)
    return c


def test_offline_stale_query_marked_not_blasted_not_no_hits(tmp_path):
    out = tmp_path / "blast"
    raw_xml = out / "raw_blast_xml"
    raw_xml.mkdir(parents=True)
    old = {"s1|ITS|p": SEQ, "s2|ITS|p": SEQ}
    (raw_xml / "ITS.xml").write_text(xml_for({k: len(v) for k, v in old.items()}))   # 舊 XML（沒側車檔）
    new = {"s1|ITS|p": SEQ, "s2|ITS|p": SEQ, "s3|ITS|p": SEQ}        # 新增 s3 沒送過 BLAST
    fasta = tmp_path / "c.fasta"
    write_fasta(fasta, new)
    nb.run_blast(fasta, out, blast_cfg(), offline=True)
    rows = {r["query_id"]: r for r in csv.DictReader((out / "blast_query_summary.csv").open())}
    assert rows["s3|ITS|p"]["review_flag"] == "NOT_BLASTED"
    assert rows["s1|ITS|p"]["review_flag"] != "NOT_BLASTED"
    tax = {r["sample"]: r for r in csv.DictReader((out / "taxonomy_decision_table.csv").open())}
    assert tax["s3"]["suggested_action"] == "NOT_BLASTED_RERUN_BLAST"


def test_sidecar_hash_rejects_same_id_same_length_but_changed_sequence(tmp_path):
    out = tmp_path / "blast"
    raw_xml = out / "raw_blast_xml"
    raw_xml.mkdir(parents=True)
    old = {"s1|ITS|p": SEQ}
    xml = raw_xml / "ITS.xml"
    xml.write_text(xml_for({"s1|ITS|p": len(SEQ)}))
    nb.write_xml_meta(xml, recs(old), {"database": "nt", "hitlist_size": 100, "entrez_query": "(none)"})
    changed = {"s1|ITS|p": "TGCA" * 100}                        # 同 ID、同長度、序列不同（孔位標錯改名後重組裝）
    ok, how = nb.verified_query_ids(xml, recs(changed))
    assert how == "exact" and ok == set()
    ok2, _ = nb.verified_query_ids(xml, recs(old))
    assert ok2 == {"s1|ITS|p"}


def test_online_stale_xml_is_renamed_and_resent(tmp_path, monkeypatch):
    out = tmp_path / "blast"
    raw_xml = out / "raw_blast_xml"
    raw_xml.mkdir(parents=True)
    xml = raw_xml / "ITS.xml"
    xml.write_text(xml_for({"s1|ITS|p": len(SEQ)}))
    nb.write_xml_meta(xml, recs({"s1|ITS|p": SEQ}), {"database": "nt", "hitlist_size": 100, "entrez_query": "(none)"})
    sent = []

    class FakeHandle:
        def read(self):
            return xml_for({"s1|ITS|p": 400})

        def close(self):
            pass

    monkeypatch.setattr(nb.NCBIWWW, "qblast", lambda *a, **k: sent.append(a) or FakeHandle())
    new = recs({"s1|ITS|p": "TGCA" * 100})
    path, ids = nb.run_blast_batch("ITS", new, out, "nt", 100, "(none)", "lab@example.org", "t", offline=False)
    assert sent and ids == {"s1|ITS|p"}
    assert list(raw_xml.glob("ITS.xml.stale-*"))                   # 舊 XML 保留為 .stale
    ok, how = nb.verified_query_ids(path, new)
    assert how == "exact" and ok == {"s1|ITS|p"}                   # 新 XML 有側車檔，序列對得上


def test_matching_xml_is_reused_without_sending(tmp_path, monkeypatch):
    out = tmp_path / "blast"
    raw_xml = out / "raw_blast_xml"
    raw_xml.mkdir(parents=True)
    xml = raw_xml / "ITS.xml"
    xml.write_text(xml_for({"s1|ITS|p": len(SEQ)}))
    nb.write_xml_meta(xml, recs({"s1|ITS|p": SEQ}), {"database": "nt", "hitlist_size": 100, "entrez_query": "(none)"})
    monkeypatch.setattr(nb.NCBIWWW, "qblast", lambda *a, **k: pytest.fail("不該重送"))
    path, ids = nb.run_blast_batch("ITS", recs({"s1|ITS|p": SEQ}), out, "nt", 100, "(none)", "lab@example.org", "t")
    assert path == xml and ids == {"s1|ITS|p"}


def test_offline_with_no_matching_xml_still_raises(tmp_path):
    out = tmp_path / "blast"
    fasta = tmp_path / "c.fasta"
    write_fasta(fasta, {"s1|ITS|p": SEQ})
    with pytest.raises(SystemExit):
        nb.run_blast(fasta, out, blast_cfg(), offline=True)


def test_consensus_not_blasted_marker_is_excluded_not_counted_as_conflict():
    lookup = cc.load_lookup(gene_lookup_path(CFG))
    rows = [
        {"sample": "X", "marker": "ITS", "high_quality_species_all": "Aspergillus niger:3",
         "species_best_pident_all": "Aspergillus niger:100", "candidate_genus": "Aspergillus",
         "review_flag": "OK", "assembly_status": "PASS"},
        {"sample": "X", "marker": "Btub", "high_quality_species_all": "", "species_best_pident_all": "",
         "candidate_genus": "", "review_flag": "NOT_BLASTED", "assembly_status": "PASS"},
    ]
    r = cc.call_sample("S1", "X", rows, ["clone"], lookup)
    assert r["Marker 數"] == "1"
    assert "Btub 未進判讀：尚未送 BLAST" in r["備註"]
    assert r["判讀狀態"] != cc.STATUS_NONE

    only = cc.call_sample("S1", "X", rows[1:], ["clone"], lookup)
    assert only["判讀狀態"] == cc.STATUS_NOT_BLASTED
    assert "重跑 BLAST" in only["下次建議 marker"]


# ------------------------------------------------------------ 壞 read 不炸整單
def test_corrupt_ab1_marks_pair_unusable_instead_of_crashing(tmp_path):
    d = tmp_path / "raw"
    d.mkdir()
    (d / "S1_s1_ITS4_01.ab1").write_bytes(b"this is not an ab1 file")
    (d / "S1_s1_ITS5_02.ab1").write_bytes(b"neither is this")
    out = tmp_path / "asm"
    spa.process([d], out, CFG["qc"], CFG["assembly"])
    rows = list(csv.DictReader((out / "assembly_summary.csv").open()))
    assert [r["assembly_status"] for r in rows] == ["FAIL_READ_UNUSABLE"]
    assert "FAIL_READ_UNUSABLE" in (out / "assembly_needs_review.txt").read_text()
    assert (out / "paired_consensus_for_blast.fasta").read_text() == ""


def test_one_bad_pair_does_not_kill_other_pairs(tmp_path, monkeypatch):
    def fake_read(path, qc=None):
        meta = spa.parse_name(path.stem)
        return {"path": path, "file": path.name, **meta, "seq": "ACGT", "qual": [40] * 4, "trimmed_len": 4,
                "mean_q": 40.0, "q20_pct": 100.0, "qc_status": "PASS"}

    def fake_pair(sample, marker, reads, outdir, asm):
        if sample == "bad":
            raise IndexError("trimmed to 0 bp")
        rec = SeqRecord(Seq("ACGT" * 80), id=f"{sample}|{marker}|p", description="")
        row = {k: "" for k in spa.ROW_FIELDS}
        row.update({"sample": sample, "marker": marker, "assembly_status": "PASS", "read_count": "2"})
        return {"base_rows": [], "low_confidence": [], "record": rec, "status": "PASS", "row": row, "review_line": ""}

    monkeypatch.setattr(spa, "read_trimmed_ab1", fake_read)
    monkeypatch.setattr(spa, "_assemble_pair", fake_pair)
    d = tmp_path / "raw"
    d.mkdir()
    for s in ("good", "bad"):
        for p in ("ITS4", "ITS5"):
            (d / f"S1_{s}_{p}_01.ab1").write_bytes(b"x")
    out = tmp_path / "asm"
    spa.process([d], out, CFG["qc"], CFG["assembly"])
    st = {r["sample"]: r["assembly_status"] for r in csv.DictReader((out / "assembly_summary.csv").open())}
    assert st == {"bad": "FAIL_READ_UNUSABLE", "good": "PASS"}
    assert ">good|ITS|p" in (out / "paired_consensus_for_blast.fasta").read_text()


def test_qc_corrupt_file_marked_not_crash(tmp_path):
    from pipeline import sanger_qc

    d = tmp_path / "raw"
    d.mkdir()
    (d / "S1_s1_ITS4_01.ab1").write_bytes(b"junk")
    out = tmp_path / "qc"
    sanger_qc.process(d, out, CFG["qc"])
    rows = list(csv.DictReader((out / "qc_summary.csv").open()))
    assert rows[0]["status"] == "FAIL_READ_UNUSABLE"


# ------------------------------------------------------------ 網頁：清除重跑
def test_clear_job_only_inside_work_base(tmp_path, monkeypatch):
    from webapp import runner

    monkeypatch.setattr(runner, "WORK_BASE", tmp_path / "jobs")
    d = runner.new_job_dir("abc123def456")
    (d / "x").write_text("1")
    assert runner.clear_job("abc123def456") is True
    assert not d.exists()
    assert runner.clear_job("abc123def456") is False
    with pytest.raises(runner.UploadError):
        runner.clear_job("../../etc")
