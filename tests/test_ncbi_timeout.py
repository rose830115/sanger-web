"""NCBI BLAST 逾時自動重發：mock qblast 傳輸層（用 Event 卡住假 qblast 模擬佇列壅塞），逾時秒數壓到 0.1 秒內。不連線。"""

import threading

import pytest

from pipeline import sanger_ncbi_blast as nb
from tests.test_p5b_blast import QID, SEQ, TYPE_Q, FakeHandle, cfg, make_xml, read, write_fasta

GOOD = make_xml({QID: [{"acc": "Z1", "defn": "Aspergillus niger ITS"}]})
FAST = dict(request_timeout_sec=0.1, resend_timeout_sec=0.2, max_resends=1)


@pytest.fixture
def stuck():
    """放行用 Event：測試結束放掉所有被卡住的背景 thread。"""
    ev = threading.Event()
    yield ev
    ev.set()


def setup(tmp_path, monkeypatch, fake):
    monkeypatch.setattr(nb.NCBIWWW, "qblast", fake)
    monkeypatch.setattr(nb.time, "sleep", lambda s: None)
    fasta = tmp_path / "c.fasta"
    write_fasta(fasta, {QID: SEQ})
    return fasta, tmp_path / "blast"


def test_timeout_then_resend_succeeds(tmp_path, monkeypatch, stuck, caplog):
    calls = {"full": 0}

    def fake(program, db, fasta, **kw):
        if kw["entrez_query"] == "(none)":
            calls["full"] += 1
            if calls["full"] == 1:          # 第一次送出被佇列卡住
                stuck.wait(5)
        return FakeHandle(GOOD)

    fasta, out = setup(tmp_path, monkeypatch, fake)
    events = []
    with caplog.at_level("WARNING", logger="sanger.blast"):
        nb.run_blast(fasta, out, cfg(email="lab@example.org", **FAST), on_event=events.append)
    assert calls["full"] == 2
    assert any("BLAST_RESENT" in r.message for r in caplog.records)
    assert events and "佇列壅塞" in events[0]
    row = read(out / "taxonomy_decision_table.csv")[0]
    assert row["type_status"] == "TYPE_OK"


def test_full_round_double_timeout_fails_batch_and_resumes(tmp_path, monkeypatch, stuck):
    sent = []
    mode = {"hang": True}

    def fake(program, db, fasta, **kw):
        sent.append(kw["entrez_query"])
        if mode["hang"]:
            stuck.wait(5)
        return FakeHandle(GOOD)

    fasta, out = setup(tmp_path, monkeypatch, fake)
    c = cfg(email="lab@example.org", **FAST)
    with pytest.raises(nb.BlastTimeout, match="TIMEOUT_ABANDONED"):
        nb.run_blast(fasta, out, c)
    assert sent == ["(none)", "(none)"]                     # 原送 + 重發一次，type 輪沒跑
    assert not (out / "raw_blast_xml" / "ITS.xml").exists()  # 失敗批次不留半套 XML
    # 斷點續跑：網路恢復後重跑同單，可正常完成
    mode["hang"] = False
    nb.run_blast(fasta, out, c)
    assert read(out / "taxonomy_decision_table.csv")[0]["type_status"] == "TYPE_OK"


def test_type_round_double_timeout_only_downgrades(tmp_path, monkeypatch, stuck):
    def fake(program, db, fasta, **kw):
        if kw["entrez_query"] == TYPE_Q:
            stuck.wait(5)
        return FakeHandle(GOOD)

    fasta, out = setup(tmp_path, monkeypatch, fake)
    nb.run_blast(fasta, out, cfg(email="lab@example.org", **FAST))
    row = read(out / "taxonomy_decision_table.csv")[0]
    assert row["type_status"] == "TYPE_NOT_RUN" and row["high_quality_species_all"] == "Aspergillus niger:1"
    meta = read(out / "blast_run_metadata.csv")[0]
    assert "TIMEOUT_ABANDONED" in meta["type_round_errors"] and meta["type_round_timeouts"] == "ITS"


def test_network_failure_not_flagged_as_timeout(tmp_path, monkeypatch):
    def fake(program, db, fasta, **kw):
        if kw["entrez_query"] == TYPE_Q:
            raise RuntimeError("NCBI busy")
        return FakeHandle(GOOD)

    fasta, out = setup(tmp_path, monkeypatch, fake)
    nb.run_blast(fasta, out, cfg(email="lab@example.org", **FAST))
    meta = read(out / "blast_run_metadata.csv")[0]
    assert "NCBI busy" in meta["type_round_errors"] and "TIMEOUT" not in meta["type_round_errors"]
    assert meta["type_round_timeouts"] == ""


def test_zero_timeout_disables_watchdog():
    assert nb.call_with_timeout_resend(lambda: 7, 0, 0) == 7
