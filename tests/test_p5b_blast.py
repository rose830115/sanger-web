"""P5b BLAST 層：🔴-5 type 第二輪、🔴-8 hitlist／同分合併、🟡-1＋🟢-1 taxid 正規化。

全部離線：用手寫的最小 NCBI XML＋mock（qblast、E-utilities），不連線。
"""

import csv
import json
from datetime import date, timedelta
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pipeline import consensus_call as cc
from pipeline import sanger_ncbi_blast as nb
from pipeline import taxonomy as tx
from pipeline.config import gene_lookup_path, load_config

CFG = load_config()
SEQ = "ACGT" * 100          # 400 bp


# ------------------------------------------------------------ 假資料
def hit_xml(hits: list[dict], qlen: int) -> str:
    out = ""
    for i, h in enumerate(hits, start=1):
        alen = h.get("alen", qlen)
        ident = h.get("ident", alen)
        out += f"""<Hit><Hit_num>{i}</Hit_num><Hit_id>{h.get('hit_id', 'gb|' + h['acc'] + '.1|')}</Hit_id>
<Hit_def>{h['defn']}</Hit_def><Hit_accession>{h['acc']}</Hit_accession><Hit_len>{h.get('hlen', qlen)}</Hit_len>
<Hit_hsps><Hsp><Hsp_num>1</Hsp_num><Hsp_bit-score>{h.get('bits', 1000)}</Hsp_bit-score><Hsp_score>500</Hsp_score>
<Hsp_evalue>0</Hsp_evalue><Hsp_query-from>1</Hsp_query-from><Hsp_query-to>{alen}</Hsp_query-to>
<Hsp_hit-from>1</Hsp_hit-from><Hsp_hit-to>{alen}</Hsp_hit-to><Hsp_query-frame>1</Hsp_query-frame>
<Hsp_hit-frame>1</Hsp_hit-frame><Hsp_identity>{ident}</Hsp_identity><Hsp_positive>{ident}</Hsp_positive>
<Hsp_gaps>0</Hsp_gaps><Hsp_align-len>{alen}</Hsp_align-len><Hsp_qseq>AAAAA</Hsp_qseq><Hsp_hseq>AAAAA</Hsp_hseq>
<Hsp_midline>|||||</Hsp_midline></Hsp></Hit_hsps></Hit>"""
    return out


def make_xml(per_query: dict[str, list[dict]], qlen: int = 400) -> str:
    its = ""
    for i, (qid, hits) in enumerate(per_query.items(), start=1):
        its += f"""<Iteration><Iteration_iter-num>{i}</Iteration_iter-num><Iteration_query-ID>Query_{i}</Iteration_query-ID>
<Iteration_query-def>{qid}</Iteration_query-def><Iteration_query-len>{qlen}</Iteration_query-len>
<Iteration_hits>{hit_xml(hits, qlen)}</Iteration_hits>
<Iteration_stat><Statistics><Statistics_db-num>1</Statistics_db-num><Statistics_db-len>1</Statistics_db-len>
<Statistics_hsp-len>0</Statistics_hsp-len><Statistics_eff-space>0</Statistics_eff-space><Statistics_kappa>0.1</Statistics_kappa>
<Statistics_lambda>1.0</Statistics_lambda><Statistics_entropy>0.1</Statistics_entropy></Statistics></Iteration_stat></Iteration>"""
    return f"""<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput><BlastOutput_program>blastn</BlastOutput_program><BlastOutput_version>BLASTN 2.16.0+</BlastOutput_version>
<BlastOutput_reference>x</BlastOutput_reference><BlastOutput_db>nt</BlastOutput_db><BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
<BlastOutput_query-def>x</BlastOutput_query-def><BlastOutput_query-len>1</BlastOutput_query-len>
<BlastOutput_param><Parameters><Parameters_expect>10</Parameters_expect></Parameters></BlastOutput_param>
<BlastOutput_iterations>{its}</BlastOutput_iterations></BlastOutput>"""


def recs(spec: dict[str, str]):
    return [SeqRecord(Seq(seq), id=qid, description="") for qid, seq in spec.items()]


def write_fasta(path: Path, spec: dict[str, str]):
    path.write_text("".join(f">{k}\n{v}\n" for k, v in spec.items()))


def put_xml(out: Path, sub: str, marker: str, per_query: dict[str, list[dict]], spec: dict[str, str], hitlist=1000,
            entrez="(none)"):
    d = out / sub
    d.mkdir(parents=True, exist_ok=True)
    xml = d / f"{marker}.xml"
    xml.write_text(make_xml(per_query))
    nb.write_xml_meta(xml, recs(spec), {"database": "nt", "hitlist_size": hitlist, "entrez_query": entrez})


def cfg(**kw):
    c = dict(CFG["blast"])
    c.update(delay_between_markers_sec=0, delay_between_orders_sec=0, taxonomy_lookup=False)
    c.update(kw)
    return c


def read(path: Path):
    return list(csv.DictReader(path.open(encoding="utf-8-sig")))


QID = "s1|ITS|p"
TYPE_Q = "sequence from type[filter]"


def run_offline(tmp_path, full_hits, type_hits, c=None):
    out = tmp_path / "blast"
    spec = {QID: SEQ}
    put_xml(out, "raw_blast_xml", "ITS", {QID: full_hits}, spec)
    if type_hits is not None:
        put_xml(out, "raw_blast_xml_type", "ITS", {QID: type_hits}, spec, entrez=TYPE_Q)
    fasta = tmp_path / "c.fasta"
    write_fasta(fasta, spec)
    nb.run_blast(fasta, out, c or cfg(), offline=True)
    return out


# ------------------------------------------------------------ 🔴-5 type 第二輪
def test_type_ok_is_primary_and_full_db_is_kept_as_auxiliary(tmp_path):
    # 全庫：誤標的 Trichoderma lixii 與正確的 Penicillium chrysogenum 同分；type：只有 P. chrysogenum
    full = [{"acc": "A1", "defn": "Penicillium chrysogenum ITS"}, {"acc": "A2", "defn": "Trichoderma lixii ITS"}]
    typ = [{"acc": "T1", "defn": "Penicillium chrysogenum strain CBS 306.48 ITS"}]
    out = run_offline(tmp_path, full, typ)
    row = read(out / "taxonomy_decision_table.csv")[0]
    assert row["type_status"] == "TYPE_OK" and row["candidate_source"] == "type"
    assert row["high_quality_species_all"] == "Penicillium chrysogenum:1"          # 主判＝type
    assert "Trichoderma lixii" in row["full_species_best_pident_all"]              # 全庫為輔，另存
    assert [r["ref_source"] for r in read(out / "blast_type_hits.csv")] == ["type"]
    assert [r["ref_source"] for r in read(out / "blast_top100_hits.csv")] == ["all", "all"]


def test_type_second_round_no_hits_only_lowers_confidence(tmp_path):
    full = [{"acc": "A1", "defn": "Aspergillus niger ITS"}]
    typ = [{"acc": "T1", "defn": "Aspergillus tubingensis ITS", "ident": 380}]          # 95% -> 非高品質
    out = run_offline(tmp_path, full, typ)
    row = read(out / "taxonomy_decision_table.csv")[0]
    assert row["type_status"] == "TYPE_NO_HITS" and row["candidate_source"] == "全庫"
    assert row["high_quality_species_all"] == "Aspergillus niger:1"              # 沿用全庫候選，不判錯
    assert "NO_TYPE_SUPPORT" in row["review_flag"]
    assert row["suggested_action"] == "TYPE_SUPPORT_MISSING_REVIEW"


def test_type_not_run_when_no_type_xml_offline(tmp_path):
    out = run_offline(tmp_path, [{"acc": "A1", "defn": "Aspergillus niger ITS"}], None)
    row = read(out / "taxonomy_decision_table.csv")[0]
    assert row["type_status"] == "TYPE_NOT_RUN" and "TYPE_NOT_RUN" in row["review_flag"]


def test_type_round_can_be_disabled(tmp_path):
    out = run_offline(tmp_path, [{"acc": "A1", "defn": "Aspergillus niger ITS"}], None, cfg(type_round=False))
    row = read(out / "taxonomy_decision_table.csv")[0]
    assert row["type_status"] == "TYPE_DISABLED" and row["review_flag"] == "OK"
    assert not (out / "blast_type_hits.csv").exists()


class FakeHandle:
    def __init__(self, text):
        self.text = text

    def read(self):
        return self.text

    def close(self):
        pass


def test_online_sends_two_rounds_with_type_filter(tmp_path, monkeypatch):
    sent = []

    def fake_qblast(program, db, fasta, **kw):
        sent.append(kw["entrez_query"])
        defn = "Penicillium chrysogenum ITS"
        return FakeHandle(make_xml({QID: [{"acc": "Z1", "defn": defn}]}))

    monkeypatch.setattr(nb.NCBIWWW, "qblast", fake_qblast)
    monkeypatch.setattr(nb.time, "sleep", lambda s: None)
    out = tmp_path / "blast"
    fasta = tmp_path / "c.fasta"
    write_fasta(fasta, {QID: SEQ})
    nb.run_blast(fasta, out, cfg(email="lab@example.org"), offline=False)
    assert sent == ["(none)", TYPE_Q]
    assert read(out / "taxonomy_decision_table.csv")[0]["type_status"] == "TYPE_OK"
    meta = json.loads((out / "raw_blast_xml_type" / "ITS.xml.meta.json").read_text())
    assert meta["entrez_query"] == TYPE_Q and meta["hitlist_size"] == cfg()["hitlist_size"]


def test_type_round_failure_does_not_abort_whole_run(tmp_path, monkeypatch):
    calls = []

    def fake_qblast(program, db, fasta, **kw):
        calls.append(kw["entrez_query"])
        if kw["entrez_query"] == TYPE_Q:
            raise RuntimeError("NCBI busy")
        return FakeHandle(make_xml({QID: [{"acc": "Z1", "defn": "Aspergillus niger ITS"}]}))

    monkeypatch.setattr(nb.NCBIWWW, "qblast", fake_qblast)
    monkeypatch.setattr(nb.time, "sleep", lambda s: None)
    out = tmp_path / "blast"
    fasta = tmp_path / "c.fasta"
    write_fasta(fasta, {QID: SEQ})
    nb.run_blast(fasta, out, cfg(email="lab@example.org"), offline=False)
    row = read(out / "taxonomy_decision_table.csv")[0]
    assert row["type_status"] == "TYPE_NOT_RUN" and row["high_quality_species_all"] == "Aspergillus niger:1"
    assert "NCBI busy" in read(out / "blast_run_metadata.csv")[0]["type_round_errors"]


def test_legacy_xml_without_sidecar_is_resent_online_but_kept_offline(tmp_path, monkeypatch):
    out = tmp_path / "blast"
    d = out / "raw_blast_xml"
    d.mkdir(parents=True)
    (d / "ITS.xml").write_text(make_xml({QID: [{"acc": "A1", "defn": "Aspergillus niger ITS"}]}))
    r = recs({QID: SEQ})
    path, ids = nb.run_blast_batch("ITS", r, out, "nt", 1000, "(none)", None, "t", offline=True)
    assert ids == {QID}                                                    # 離線仍沿用
    sent = []
    monkeypatch.setattr(nb.NCBIWWW, "qblast",
                        lambda *a, **k: sent.append(1) or FakeHandle(make_xml({QID: [{"acc": "B", "defn": "X y"}]})))
    nb.run_blast_batch("ITS", r, out, "nt", 1000, "(none)", "lab@example.org", "t", offline=False)
    assert sent                                                            # 線上：無法確認當時的 hitlist／篩選 -> 重送
    assert list(d.glob("ITS.xml.stale-*"))


# ---- 交集：type 主判、全庫為輔
def crow(marker, hq, rank, type_status="TYPE_OK", src="type", full_rank="", flag="OK"):
    return {
        "sample": "X1", "marker": marker, "high_quality_species_all": hq, "species_best_pident_all": rank,
        "candidate_genus": "Penicillium", "review_flag": flag, "assembly_status": "PASS",
        "type_status": type_status, "candidate_source": src, "full_species_best_pident_all": full_rank,
    }


EXCL = CFG["consensus"]["exclude_keywords"]
LOOKUP = cc.load_lookup(gene_lookup_path(CFG))


def call(rows):
    return cc.call_sample("S1", "X1", rows, EXCL, LOOKUP, rank_cfg=CFG["consensus"])


def test_consensus_unique_with_type_support_on_all_markers_is_green():
    r = call([crow("ITS", "Penicillium chrysogenum:2", "Penicillium chrysogenum:100.000"),
              crow("CF", "Penicillium chrysogenum:1", "Penicillium chrysogenum:100.000; Penicillium rubens:98.500")])
    assert r["判讀狀態"] == cc.STATUS_UNIQUE
    assert r["判讀來源"] == "CF:type；ITS:type"


@pytest.mark.parametrize("state,why", [("TYPE_NO_HITS", "無 type 序列支持"), ("TYPE_NOT_RUN", "未做 type 第二輪")])
def test_consensus_missing_type_support_downgrades_to_pending_not_wrong(state, why):
    r = call([crow("ITS", "Penicillium chrysogenum:2", "Penicillium chrysogenum:100.000"),
              crow("CF", "Penicillium chrysogenum:1", "Penicillium chrysogenum:100.000; Penicillium rubens:98.500",
                   type_status=state, src="全庫")])
    assert r["判讀狀態"] == cc.STATUS_INTERSECT_PENDING            # 降信心，不是判錯（仍給出唯一物種名）
    assert r["唯一物種"] == "Penicillium chrysogenum" and why in r["待確認原因"]
    assert "CF:全庫" in r["判讀來源"]


def test_consensus_legacy_rows_without_type_columns_are_not_downgraded():
    legacy = [{"sample": "X1", "marker": m, "high_quality_species_all": "Penicillium chrysogenum:2",
               "species_best_pident_all": "Penicillium chrysogenum:100.000; Penicillium rubens:98.000",
               "candidate_genus": "Penicillium", "review_flag": "OK", "assembly_status": "PASS"}
              for m in ("ITS", "CF")]
    r = call(legacy)
    assert r["判讀狀態"] == cc.STATUS_UNIQUE and "無 type 欄位" in r["判讀來源"]


def test_consensus_notes_when_full_db_has_higher_non_type_species():
    r = call([crow("ITS", "Penicillium chrysogenum:2", "Penicillium chrysogenum:99.500",
                   full_rank="Trichoderma lixii:100.000; Penicillium chrysogenum:99.500"),
              crow("CF", "Penicillium chrysogenum:1", "Penicillium chrysogenum:100.000; Penicillium rubens:98.500")])
    assert r["判讀狀態"] == cc.STATUS_UNIQUE
    assert "全庫另有更高分 Trichoderma lixii" in r["備註"]


# ------------------------------------------------------------ 🔴-8 hitlist／同分合併
def test_default_hitlist_is_1000_and_sent_to_ncbi(tmp_path, monkeypatch):
    assert CFG["blast"]["hitlist_size"] == 1000
    seen = {}

    def fake_qblast(program, db, fasta, **kw):
        seen.update(kw)
        return FakeHandle(make_xml({QID: [{"acc": "Z1", "defn": "Aspergillus niger ITS"}]}))

    monkeypatch.setattr(nb.NCBIWWW, "qblast", fake_qblast)
    monkeypatch.setattr(nb.time, "sleep", lambda s: None)
    fasta = tmp_path / "c.fasta"
    write_fasta(fasta, {QID: SEQ})
    nb.run_blast(fasta, tmp_path / "o", cfg(email="lab@example.org", type_round=False), offline=False)
    assert seen["hitlist_size"] == seen["descriptions"] == seen["alignments"] == 1000


def test_candidate_set_is_not_cut_at_100_hits(tmp_path):
    hits = [{"acc": f"A{i}", "defn": f"Aspergillus species{chr(97 + i % 26)}{chr(97 + i // 26)} ITS"} for i in range(150)]
    out = run_offline(tmp_path, hits, None, cfg(type_round=False))
    row = read(out / "taxonomy_decision_table.csv")[0]
    assert len(cc.parse_counts(row["high_quality_species_all"])) == 150
    assert "HITLIST_SATURATED" not in row["review_flag"]              # 清單 1000 沒塞滿


def test_hit_records_parses_old_and_new_defline_formats():
    class Al:
        accession = "OW984504"
        hit_id = "gi|2259517814|emb|OW984504.1|"
        hit_def = ("Aspergillus udagawae genomic DNA >gi|2259531846|emb|OW988077.1| Aspergillus fumigatus ITS"
                   " >MK123456.1 Penicillium rubens ITS")

    recs_ = nb.hit_records(Al())
    assert recs_[0] == ("OW984504.1", "Aspergillus udagawae genomic DNA")
    assert recs_[1] == ("OW988077.1", "Aspergillus fumigatus ITS")
    assert recs_[2] == ("MK123456.1", "Penicillium rubens ITS")


def test_identical_sequence_records_merged_into_one_hit_all_count_as_candidates(tmp_path):
    # BLAST 把序列相同的 3 筆紀錄收成一個 hit：A. niger／A. tubingensis／A. welwitschiae 同分，都要進候選集合
    defn = "Aspergillus niger ITS &gt;gb|AB2.1| Aspergillus tubingensis ITS &gt;gb|AB3.1| Aspergillus welwitschiae ITS"
    out = run_offline(tmp_path, [{"acc": "AB1", "defn": defn}], None, cfg(type_round=False))
    row = read(out / "taxonomy_decision_table.csv")[0]
    assert set(cc.parse_counts(row["high_quality_species_all"])) == {
        "Aspergillus niger", "Aspergillus tubingensis", "Aspergillus welwitschiae"}
    best = cc.parse_scores(row["species_best_pident_all"])
    assert best["Aspergillus tubingensis"] == best["Aspergillus niger"] == 100.0
    hit = read(out / "blast_top100_hits.csv")[0]
    assert hit["merged_records"] == "2" and "_records" not in hit


def test_saturated_hitlist_flag_blocks_green():
    hits = [{"acc": f"A{i}", "defn": f"Aspergillus niger strain{i}"} for i in range(3)]
    sample = []
    for h in hits:
        sample.append({"query_id": QID, "sample": "s1", "marker": "ITS", "primers": "p", "query_length": "400",
                       "rank": str(len(sample) + 1), "species_guess": "Aspergillus niger", "genus_guess": "Aspergillus",
                       "organism_guess": "Aspergillus niger", "accession": h["acc"], "pident": "100.000",
                       "qcov_pct": "100.000", "align_len": "400", "evalue": "0", "merged_species": ""})
    full_sat = nb.summarize_query(QID, sample, 99.0, 95.0, 300, hitlist_size=3)
    assert "HITLIST_SATURATED" in full_sat["review_flag"]
    assert "HITLIST_SATURATED" not in nb.summarize_query(QID, sample, 99.0, 95.0, 300, hitlist_size=1000)["review_flag"]
    assert "HITLIST_SATURATED" in cc.DEFAULT_BLOCKING_FLAGS
    assert "HITLIST_SATURATED" in CFG["consensus"]["blocking_flags"]


# ------------------------------------------------------------ 🟡-1＋🔴-5＋🟢-1 taxid 正規化
def tnode(taxid, name, rank, lineage=(), aka=()):
    lin = "".join(f"<Taxon><TaxId>{i}</TaxId><ScientificName>{n}</ScientificName><Rank>{r}</Rank></Taxon>"
                  for i, n, r in lineage)
    akas = "<AkaTaxIds>" + "".join(f"<TaxId>{a}</TaxId>" for a in aka) + "</AkaTaxIds>" if aka else ""
    return (f"<Taxon><TaxId>{taxid}</TaxId><ScientificName>{name}</ScientificName><Rank>{rank}</Rank>"
            f"{akas}<LineageEx>{lin}</LineageEx></Taxon>")


FUNGI = [("4751", "Fungi", "kingdom")]
TRICH = ("5543", "Trichoderma", "genus")
ASP = ("5052", "Aspergillus", "genus")
TAXA = {
    # Hypocrea lixii 的舊 taxid 已併入 Trichoderma lixii（AkaTaxIds），兩個 accession 查到同一個現行名稱
    "111": tnode("222", "Trichoderma lixii", "species", FUNGI + [TRICH], aka=["111"]),
    "222": tnode("222", "Trichoderma lixii", "species", FUNGI + [TRICH], aka=["111"]),
    "300": tnode("300", "Aspergillus niger CBS 513.88", "strain",
                 FUNGI + [ASP, ("5061", "Aspergillus niger", "species")]),
    "400": tnode("400", "Aspergillus section Nigri", "section", FUNGI + [ASP]),
    "500": tnode("500", "Aspergillus niger species complex", "species group", FUNGI + [ASP]),
    "600": tnode("600", "Aspergillus flavus x Aspergillus oryzae", "species", FUNGI + [ASP]),
    "700": tnode("700", "synthetic construct", "species", [("28384", "other sequences", "no rank")]),
    "800": tnode("800", "Penicillium rubens", "species", FUNGI + [("5073", "Penicillium", "genus")]),
    "900": tnode("900", "Aspergillus sp.", "species", FUNGI + [ASP]),
}
ACC2TAX = {"OLD1.1": "111", "NEW1.1": "222", "STR1.1": "300", "SEC1.1": "400", "CPX1.1": "500", "HYB1.1": "600",
           "SYN1.1": "700", "RUB1.1": "800", "SPP1.1": "900"}


class FakeEutils:
    def __init__(self):
        self.calls = []

    def __call__(self, url, params):
        self.calls.append((url.rsplit("/", 1)[-1], dict(params)))
        if "esummary" in url:
            ids = params["id"].split(",")
            res = {"uids": [], }
            for i, acc in enumerate(ids):
                if acc in ACC2TAX:
                    res["uids"].append(str(i))
                    res[str(i)] = {"accessionversion": acc, "caption": acc.split(".")[0], "taxid": int(ACC2TAX[acc])}
            return json.dumps({"result": res}).encode()
        wanted = params["id"].split(",")
        seen, body = set(), ""
        for t in wanted:
            xml = TAXA[t]
            if xml not in seen:
                seen.add(xml)
                body += xml
        return f"<?xml version='1.0'?><TaxaSet>{body}</TaxaSet>".encode()


def mk_resolver(tmp_path, offline=False, http=None, **kw):
    sleeps = []
    r = tx.TaxonomyResolver(cache_path=tmp_path / "cache" / "tax.json", email="lab@example.org", offline=offline,
                            http=http or FakeEutils(), sleep=sleeps.append, **kw)
    r.sleeps = sleeps
    return r


def test_classify_taxon_species_strain_section_complex_hybrid_synthetic(tmp_path):
    r = mk_resolver(tmp_path)
    info = r.resolve(["OLD1.1", "NEW1.1", "STR1.1", "SEC1.1", "CPX1.1", "HYB1.1", "SYN1.1", "RUB1.1"])
    assert info["OLD1.1"]["species"] == info["NEW1.1"]["species"] == "Trichoderma lixii"        # 舊名併入現行名
    assert info["OLD1.1"]["taxid"] == "222" or info["OLD1.1"]["species"] == "Trichoderma lixii"
    assert info["STR1.1"]["species"] == "Aspergillus niger" and info["STR1.1"]["genus"] == "Aspergillus"   # 品系 -> 種
    assert info["SEC1.1"]["species"] == "" and info["SEC1.1"]["genus"] == "Aspergillus"        # section 不當種名（🟢-1）
    assert info["CPX1.1"]["species"] == ""                                                      # species group／complex
    assert info["HYB1.1"]["species"] == "Aspergillus sp."                                       # hybrid -> 不當二名法種
    assert info["SYN1.1"]["species"] == ""                                                      # Synthetic construct
    assert info["RUB1.1"]["species"] == "Penicillium rubens"


def test_resolver_caches_to_file_and_second_run_makes_no_requests(tmp_path):
    http = FakeEutils()
    r = mk_resolver(tmp_path, http=http)
    r.resolve(["RUB1.1", "STR1.1"])
    n = len(http.calls)
    assert n == 2                                                   # 1 次 esummary＋1 次 taxonomy efetch（批次合併）
    cache_text = (tmp_path / "cache" / "tax.json").read_text()
    assert "lab@example.org" not in cache_text                      # 信箱不落檔
    assert all(c[1].get("email") == "lab@example.org" for c in http.calls)
    http2 = FakeEutils()
    r2 = mk_resolver(tmp_path, http=http2)                          # 新 resolver 讀同一份快取檔
    assert r2.resolve(["RUB1.1", "STR1.1"])["RUB1.1"]["species"] == "Penicillium rubens"
    assert http2.calls == []


def test_resolver_offline_never_calls_network_and_returns_none_for_unknown(tmp_path):
    http = FakeEutils()
    r = mk_resolver(tmp_path, offline=True, http=http)
    assert r.resolve(["RUB1.1"]) == {"RUB1.1": None}
    assert http.calls == [] and r.stats["unresolved"] == 1


def test_resolver_refetches_names_older_than_max_age_but_keeps_accession_map(tmp_path):
    http = FakeEutils()
    r = mk_resolver(tmp_path, http=http, max_age_days=90, today=lambda: date(2026, 1, 1))
    r.resolve(["RUB1.1"])
    http2 = FakeEutils()
    r2 = mk_resolver(tmp_path, http=http2, max_age_days=90, today=lambda: date(2026, 1, 1) + timedelta(days=120))
    r2.resolve(["RUB1.1"])
    assert [c[0] for c in http2.calls] == ["efetch.fcgi"]           # accession->taxid 不變，只重查 taxid->名稱


def test_resolver_rate_limited_and_batched(tmp_path):
    http = FakeEutils()
    r = mk_resolver(tmp_path, http=http)
    r.api_key = None
    r.interval = 0.4
    r.resolve([f"X{i}.1" for i in range(450)])                      # 450 個查不到的 accession -> 3 批 esummary
    assert [c[0] for c in http.calls] == ["esummary.fcgi"] * 3
    assert len(r.sleeps) >= 2 and all(0 < s <= 0.4 for s in r.sleeps)


def test_resolver_survives_network_error_and_falls_back(tmp_path):
    def boom(url, params):
        raise RuntimeError("NCBI 503")

    r = mk_resolver(tmp_path, http=boom)
    assert r.resolve(["RUB1.1"]) == {"RUB1.1": None} and "503" in r.errors[0]


def two_hit_run(tmp_path, resolver, hits=None, c=None):
    hits = hits or [{"acc": "OLD1", "hit_id": "gb|OLD1.1|", "defn": "Hypocrea lixii ITS"},
                    {"acc": "NEW1", "hit_id": "gb|NEW1.1|", "defn": "Trichoderma lixii ITS"}]
    out = tmp_path / "blast"
    spec = {QID: SEQ}
    put_xml(out, "raw_blast_xml", "ITS", {QID: hits}, spec)
    fasta = tmp_path / "c.fasta"
    write_fasta(fasta, spec)
    nb.run_blast(fasta, out, c or cfg(type_round=False, taxonomy_lookup=True), offline=True, resolver=resolver)
    return out


def test_synonyms_are_normalized_to_current_name_before_intersection(tmp_path):
    out = two_hit_run(tmp_path, mk_resolver(tmp_path))
    hits = read(out / "blast_top100_hits.csv")
    assert [h["species_guess"] for h in hits] == ["Trichoderma lixii"] * 2
    assert hits[0]["title_species_guess"] == "Hypocrea lixii" and hits[0]["name_source"] == "taxonomy"
    assert hits[0]["taxid"] == "222"
    row = read(out / "taxonomy_decision_table.csv")[0]
    assert row["high_quality_species_all"] == "Trichoderma lixii:2"      # 沒正規化會是兩個種
    res = read(out / "taxonomy_resolution.csv")
    assert {r["accession"] for r in res} == {"OLD1.1", "NEW1.1"} and res[0]["resolved"] == "yes"


def test_unresolved_accession_falls_back_to_title_guess(tmp_path):
    out = two_hit_run(tmp_path, mk_resolver(tmp_path, offline=True),
                      hits=[{"acc": "ZZZ9", "hit_id": "gb|ZZZ9.1|", "defn": "Hypocrea lixii ITS"}])
    h = read(out / "blast_top100_hits.csv")[0]
    assert h["species_guess"] == "Hypocrea lixii" and h["name_source"] == "title" and h["taxid"] == ""
    assert read(out / "blast_run_metadata.csv")[0]["taxonomy_unresolved"] == "1"


def test_section_and_complex_hits_do_not_become_species_candidates(tmp_path):
    hits = [{"acc": "SEC1", "hit_id": "gb|SEC1.1|", "defn": "Aspergillus section Nigri ITS"},
            {"acc": "CPX1", "hit_id": "gb|CPX1.1|", "defn": "Aspergillus niger species complex"},
            {"acc": "SYN1", "hit_id": "gb|SYN1.1|", "defn": "Synthetic construct"},
            {"acc": "RUB1", "hit_id": "gb|RUB1.1|", "defn": "Penicillium rubens ITS"}]
    out = two_hit_run(tmp_path, mk_resolver(tmp_path), hits=hits)
    row = read(out / "taxonomy_decision_table.csv")[0]
    assert row["high_quality_species_all"] == "Penicillium rubens:1"


def test_taxonomy_lookup_disabled_gives_p5a_title_behavior(tmp_path):
    out = tmp_path / "blast"
    spec = {QID: SEQ}
    put_xml(out, "raw_blast_xml", "ITS", {QID: [{"acc": "OLD1", "defn": "Hypocrea lixii ITS"}]}, spec)
    fasta = tmp_path / "c.fasta"
    write_fasta(fasta, spec)
    nb.run_blast(fasta, out, cfg(type_round=False, taxonomy_lookup=False), offline=True)
    assert read(out / "blast_top100_hits.csv")[0]["name_source"] == "title"


def test_title_guess_edge_cases_section_complex_construct_family_cf():
    g = lambda t: nb.organism_guess(t)  # noqa: E731
    assert g("Aspergillus section Nigri strain X")[2] == ""                     # 只到屬
    assert g("Aspergillus fumigatus complex isolate")[2] == "Aspergillus fumigatus"
    assert g("Synthetic construct clone ABC")[2] == ""
    assert g("Hypocreales sp. isolate 12")[1] == ""                             # 目名不當屬名
    assert g("Trichocomaceae environmental sample")[1] == ""                    # 科名不當屬名
    assert g("Aspergillus cf. flavus strain A")[2] == "Aspergillus cf. flavus"  # cf. 規則不再被吃掉
    assert g("Aspergillus flavus x Aspergillus oryzae hybrid")[2] == "Aspergillus flavus"


def test_genus_alias_maps_old_genus_to_lookup_key():
    lookup = cc.load_lookup(gene_lookup_path(CFG))
    rows = [{"sample": "X1", "marker": "ITS", "high_quality_species_all": "Parengyodontium album:3",
             "species_best_pident_all": "Parengyodontium album:100", "candidate_genus": "Parengyodontium",
             "review_flag": "OK", "assembly_status": "PASS"}]
    r = cc.call_sample("S1", "X1", rows, EXCL, lookup, rank_cfg=CFG["consensus"])
    assert "需人工查證" not in r["下次建議 marker"] and r["下次建議 marker"]


def test_taxid_that_is_only_genus_sp_keeps_specific_title_name_as_conservative_candidate(tmp_path):
    hits = [{"acc": "SPP1", "hit_id": "gb|SPP1.1|", "defn": "Aspergillus pseudosclerotiorum strain X ITS"}]
    out = two_hit_run(tmp_path, mk_resolver(tmp_path), hits=hits)
    h = read(out / "blast_top100_hits.csv")[0]
    assert h["species_guess"] == "Aspergillus pseudosclerotiorum" and h["name_source"] == "title"
    assert h["taxid"] == "900"                                           # taxid 仍記在稽核欄
    assert read(out / "taxonomy_decision_table.csv")[0]["high_quality_species_all"] == "Aspergillus pseudosclerotiorum:1"


def test_merged_record_of_another_genus_is_flagged_not_hidden(tmp_path):
    # 序列相同的紀錄被併成一筆 hit，其中一筆是外屬誤標（例：Aspergillus Btub 序列被標成 Cercospora sojina）
    defn = "Aspergillus pseudoglaucus Btub &gt;gb|AB2.1| Cercospora sojina Btub"
    out = run_offline(tmp_path, [{"acc": "AB1", "defn": defn}], None, cfg(type_round=False))
    row = read(out / "taxonomy_decision_table.csv")[0]
    assert "MULTIPLE_GENERA_IN_TOP20" in row["review_flag"]
    assert set(cc.parse_counts(row["high_quality_species_all"])) == {"Aspergillus pseudoglaucus", "Cercospora sojina"}


def test_non_species_taxon_name_does_not_leak_into_organism_guess(tmp_path):
    out = two_hit_run(tmp_path, mk_resolver(tmp_path), hits=[
        {"acc": "SEC1", "hit_id": "gb|SEC1.1|", "defn": "Aspergillus section Nigri ITS"}])
    h = read(out / "blast_top100_hits.csv")[0]
    assert h["species_guess"] == "" and h["organism_guess"] == "Aspergillus"
