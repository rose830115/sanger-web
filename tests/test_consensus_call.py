"""Step 6 交集判讀規則：固定輸入 -> 預期輸出。"""

from pipeline.consensus_call import (
    NO_CANDIDATE,
    STATUS_GENUS_ONLY,
    STATUS_INTERSECT_PENDING,
    STATUS_UNIQUE_PENDING,
    STATUS_MULTI,
    STATUS_NONE,
    STATUS_SINGLE,
    STATUS_UNIQUE,
    call_sample,
    is_species_binomial,
    load_lookup,
    parse_counts,
)
from pipeline.config import gene_lookup_path, load_config

EXCL = ["clone", "gene", "genes", "genomic"]
CFG = load_config()
LOOKUP = load_lookup(gene_lookup_path(CFG))


def _default_rank(hq):
    # 預設排名：先列者 100%，其後每個遞減 1%（讓與排名無關的舊測試保持原意）
    names = parse_counts(hq)
    return "; ".join(f"{n}:{100 - i * 1.0:.3f}" for i, n in enumerate(names))


def row(marker, hq="", genus="Aspergillus", flag="OK", asm="PASS", rank=None):
    return {
        "species_best_pident_all": _default_rank(hq) if rank is None else rank,
        "assembly_status": asm,
        "sample": "X1",
        "marker": marker,
        "high_quality_species_counts": hq,
        "candidate_genus": genus,
        "review_flag": flag,
    }


def call(rows, lookup=None):
    return call_sample("S1", "X1", rows, EXCL, LOOKUP if lookup is None else lookup)


# ---------- 基本解析 ----------

def test_parse_counts():
    assert parse_counts("Aspergillus flavus:12; Aspergillus oryzae:3") == [
        "Aspergillus flavus",
        "Aspergillus oryzae",
    ]
    assert parse_counts("") == []


def test_binomial_rules():
    assert is_species_binomial("Aspergillus flavus", EXCL)
    assert is_species_binomial("Penicillium allii-sativi", EXCL)
    for bad in [
        "Aspergillus sp.", "Aspergillus sp", "Penicillium cf", "Fusarium aff.",
        "Aspergillus clone", "Penicillium genomic", "Aspergillus", "uncultured fungus",
    ]:
        assert not is_species_binomial(bad, EXCL), bad


# ---------- 交集 0 / 1 / >1 ----------

def test_intersection_one_is_unique():
    r = call([
        row("Btub", "Aspergillus sydowii:10; Aspergillus nidulans:2"),
        row("CF", "Aspergillus sydowii:9"),
    ])
    assert r["判讀狀態"] == STATUS_UNIQUE
    assert r["是否能比對到唯一物種"] == "是"
    assert r["唯一物種"] == "Aspergillus sydowii"
    assert r["Marker 數"] == "2"
    assert r["各 marker 高品質候選"] == (
        "Btub: Aspergillus sydowii | Aspergillus nidulans || CF: Aspergillus sydowii"
    )
    assert r["下次建議 marker"] == ""


def test_intersection_zero():
    r = call([
        row("ITS", "Aspergillus flavus:5"),
        row("Btub", "Aspergillus oryzae:5"),
    ])
    assert r["判讀狀態"] == STATUS_NONE
    assert r["是否能比對到唯一物種"] == "否"
    assert r["唯一物種"] == ""
    assert r["共同候選交集"] == ""


def test_intersection_multiple():
    r = call([
        row("Btub", "Aspergillus ochraceopetaliformis:4; Aspergillus flocculosus:3"),
        row("CF", "Aspergillus flocculosus:4; Aspergillus ochraceopetaliformis:2"),
    ])
    assert r["判讀狀態"] == STATUS_MULTI
    assert r["唯一物種"] == ""
    assert r["共同候選交集"] == "Aspergillus flocculosus | Aspergillus ochraceopetaliformis"


def test_three_markers_intersection():
    r = call([
        row("ITS", "Aspergillus flavus:5; Aspergillus oryzae:5"),
        row("Btub", "Aspergillus flavus:5; Aspergillus oryzae:5"),
        row("CF", "Aspergillus flavus:5"),
    ])
    assert r["判讀狀態"] == STATUS_UNIQUE
    assert r["唯一物種"] == "Aspergillus flavus"


def test_single_marker():
    r = call([row("Btub", "Aspergillus pseudoglaucus:7")])
    assert r["判讀狀態"] == STATUS_SINGLE
    assert r["是否能比對到唯一物種"] == "否"
    assert r["唯一物種"] == ""


# ---------- sp. / clone / gene 排除 ----------

def test_sp_and_clone_excluded_from_unique_call():
    # 兩 marker 都有 sp./clone，但只有 A. niger 是真的二名法；sp./clone 只進備註
    r = call([
        row("ITS", "Aspergillus niger:8; Aspergillus sp.:5; Aspergillus clone:2"),
        row("Btub", "Aspergillus niger:6; Aspergillus sp.:4"),
    ])
    assert r["判讀狀態"] == STATUS_UNIQUE
    assert r["唯一物種"] == "Aspergillus niger"
    assert "Aspergillus sp." in r["備註"]
    assert "Aspergillus clone" in r["備註"]


def test_sp_only_overlap_does_not_make_unique():
    # 兩 marker 唯一共同名稱是 sp. -> 交集為 0，不可當唯一物種
    r = call([
        row("ITS", "Penicillium sp.:9; Penicillium chrysogenum:3"),
        row("Btub", "Penicillium sp.:9; Penicillium rubens:3"),
    ])
    assert r["判讀狀態"] == STATUS_NONE
    assert r["唯一物種"] == ""


# ---------- 屬級建議 marker（gene_lookup）四類 ----------

def test_genus_only_suggests_top_priority_marker_with_primers():
    # (1) 屬在表內、只做過 ITS、BLAST 只到屬 -> 建議 priority 最高的 CF(CaM)，含 CMD5/CMD6 與備用 CF1/CF4
    r = call([row("ITS", "Aspergillus sp.:20", genus="Aspergillus", flag="NO_HIGH_QUALITY_HITS")])
    assert r["判讀狀態"] == STATUS_GENUS_ONLY      # 🔴-3：Aspergillus 的 ITS 無種級解析力
    assert "CF" in r["下次建議 marker"]
    assert "CMD5" in r["建議引子"] and "CMD6" in r["建議引子"]
    assert "CF1" in r["建議引子"] and "CF4" in r["建議引子"]
    assert "BLAST 只到屬：Aspergillus" in r["備註"]


def test_genus_only_skips_markers_already_done():
    # (2) Penicillium 第一建議 Btub 已做過 -> 改建議 CF
    r = call([
        row("ITS", "Penicillium sp.:20", genus="Penicillium"),
        row("Btub", "", genus="Penicillium", flag="NO_HIGH_QUALITY_HITS"),
    ])
    assert r["下次建議 marker"].startswith("CF")
    assert "CMD5" in r["建議引子"]
    assert "Btub" not in r["下次建議 marker"]


def test_genus_only_unknown_genus_asks_for_manual_lookup():
    # (3) 屬不在 gene_lookup -> 不猜，標人工查證
    r = call([row("ITS", "Zzzgenus sp.:20", genus="Zzzgenus")])
    assert "人工查證" in r["下次建議 marker"]
    assert r["建議引子"] == ""


def test_genus_only_all_markers_done():
    # (4) 表內所有 marker 都做過 -> 標人工判斷
    lookup = {
        "Testgenus": {
            "its_sufficient": False,
            "markers": [
                {"priority": 1, "marker": "Btub", "gene": "BenA",
                 "primers": {"fwd": {"name": "F", "seq": "AAA"}, "rev": {"name": "R", "seq": "TTT"}}}
            ],
        }
    }
    r = call([row("Btub", "Testgenus sp.:9", genus="Testgenus")], lookup=lookup)
    assert "人工判斷" in r["下次建議 marker"]


def test_species_level_single_marker_now_gets_suggestion():
    # 🟡-3：以前「ITS 有種名候選」會讓黃框空白；現在只要沒有唯一物種就給建議
    r = call([row("ITS", "Aspergillus flavus:9", genus="Aspergillus")])
    assert r["下次建議 marker"].startswith("CF")
    assert "CMD5" in r["建議引子"]


def test_empty_marker_candidates_counted_and_labelled():
    r = call([
        row("ITS", "Aspergillus flavus:9"),
        row("Btub", "", genus="Aspergillus", flag="NO_HIGH_QUALITY_HITS"),
    ])
    assert r["Marker 數"] == "2"
    assert r["判讀狀態"] == STATUS_NONE
    assert NO_CANDIDATE in r["各 marker 高品質候選"]


# ---------- 🔴-1：交集用完整候選欄，不被「前 8 名顯示欄」截斷 ----------

def test_intersection_uses_full_candidate_column_beyond_8():
    # ITS 有 9 個候選，第 9 名 (luchuensis) 在顯示欄被截掉；Btub 只有 luchuensis 與 costaricensis
    its_all = ["Aspergillus niger", "Aspergillus tubingensis", "Aspergillus welwitschiae",
               "Aspergillus awamori", "Aspergillus brasiliensis", "Aspergillus piperis",
               "Aspergillus vadensis", "Aspergillus costaricensis", "Aspergillus luchuensis"]
    its_display = "; ".join(f"{n}:3" for n in its_all[:8])      # 顯示欄只有前 8 種
    its_full = "; ".join(f"{n}:3" for n in its_all)
    btub = "Aspergillus luchuensis:4; Aspergillus costaricensis:4"
    rows_ = [
        {**row("ITS", its_display), "high_quality_species_all": its_full},
        {**row("Btub", btub), "high_quality_species_all": btub},
    ]
    r = call(rows_)
    assert r["判讀狀態"] == STATUS_MULTI
    assert r["唯一物種"] == ""
    assert r["共同候選交集"] == "Aspergillus costaricensis | Aspergillus luchuensis"


def test_old_format_without_full_column_warns_when_8_or_more():
    display = "; ".join(f"{n}:1" for n in [f"Aspergillus s{chr(97 + i)}" for i in range(8)])
    r = call([row("ITS", display), row("Btub", "Aspergillus sa:2")])
    assert "舊格式" in r["備註"]


def test_summarize_query_writes_full_and_display_columns():
    from pipeline.sanger_ncbi_blast import summarize_query

    hits = []
    for i in range(10):
        hits.append({
            "rank": str(i + 1), "pident": "100.0", "qcov_pct": "100.0",
            "species_guess": f"Aspergillus s{chr(97 + i)}", "genus_guess": "Aspergillus",
            "organism_guess": f"Aspergillus s{chr(97 + i)}", "accession": f"A{i}", "evalue": "0",
        })
    out = summarize_query("S|ITS|p", hits, 99.0, 95.0)
    assert len(parse_counts(out["high_quality_species_counts"])) == 8
    assert len(parse_counts(out["high_quality_species_all"])) == 10


# ---------- 🔴-2：綠燈閘門 ----------

def test_unique_but_assembly_warn_is_pending_not_green():
    r = call([
        row("Btub", "Aspergillus clavatus:5"),
        row("CF", "Aspergillus clavatus:5", asm="WARN_CONSENSUS_CHECK"),
    ])
    assert r["判讀狀態"] == STATUS_UNIQUE_PENDING
    assert r["唯一物種"] == "Aspergillus clavatus"        # 名字仍顯示，但不是綠燈
    assert r["是否能比對到唯一物種"] == "待確認"
    assert "CF 組裝 WARN_CONSENSUS_CHECK" in r["待確認原因"]
    assert "色譜圖" in r["備註"]


def test_unique_with_blocking_flag_is_pending():
    r = call([
        row("Btub", "Aspergillus clavatus:5", flag="LOW_TOP_QCOV"),
        row("CF", "Aspergillus clavatus:5"),
    ])
    assert r["判讀狀態"] == STATUS_UNIQUE_PENDING
    assert "Btub 旗標 LOW_TOP_QCOV" in r["待確認原因"]


def test_non_blocking_flag_stays_green():
    r = call([
        row("Btub", "Aspergillus clavatus:5", flag="MULTIPLE_SPECIES_IN_TOP20"),
        row("CF", "Aspergillus clavatus:5"),
    ])
    assert r["判讀狀態"] == STATUS_UNIQUE


def test_missing_assembly_status_is_not_green():
    rows_ = [row("Btub", "Aspergillus clavatus:5"), row("CF", "Aspergillus clavatus:5")]
    rows_[0]["assembly_status"] = ""
    assert call(rows_)["判讀狀態"] == STATUS_UNIQUE_PENDING


def test_unrun_marker_noted_not_hidden():
    from pipeline.consensus_call import call_sample

    r = call_sample("S1", "X1", [row("ITS", "Aspergillus flavus:5")], EXCL, LOOKUP,
                    unrun={"Btub": "組裝 FAIL_PAIRING"})
    assert "Btub 未進判讀：組裝 FAIL_PAIRING" in r["備註"]


# ---------- 🔴-3：ITS 分不開的屬，只有 ITS 不給種名 ----------

def test_its_only_in_non_sufficient_genus_gives_genus_only():
    r = call([row("ITS", "Penicillium chrysogenum:20", genus="Penicillium")])
    assert r["判讀狀態"] == STATUS_GENUS_ONLY
    assert r["唯一物種"] == "" and r["共同候選交集"] == ""
    assert r["判讀屬"] == "Penicillium"
    assert "its_sufficient=false" in r["待確認原因"]


def test_its_only_in_sufficient_genus_keeps_single_status():
    r = call([row("ITS", "Lichtheimia corymbifera:9", genus="Lichtheimia")])
    assert r["判讀狀態"] == STATUS_SINGLE
    assert r["共同候選交集"] == "Lichtheimia corymbifera"


def test_non_its_single_marker_keeps_single_status():
    r = call([row("Btub", "Penicillium rubens:7", genus="Penicillium")])
    assert r["判讀狀態"] == STATUS_SINGLE


# ---------- 🔴-4：唯一物種三條件（最高分＋與次佳種差距） ----------

def test_unique_requires_top_score_in_every_marker():
    # 交集只有 niger，但 Btub 最高分其實是 tubingensis（99.9 vs 99.1）
    r = call([
        row("ITS", "Aspergillus niger:5", rank="Aspergillus niger:100; Aspergillus tubingensis:98"),
        row("Btub", "Aspergillus niger:3; Aspergillus tubingensis:3",
            rank="Aspergillus tubingensis:99.9; Aspergillus niger:99.1"),
    ])
    assert r["判讀狀態"] == STATUS_INTERSECT_PENDING
    assert r["唯一物種"] == "Aspergillus niger"
    assert r["是否能比對到唯一物種"] == "待覆核"
    assert "最高分是 Aspergillus tubingensis" in r["待確認原因"]


def test_unique_within_tolerance_of_top_is_ok():
    r = call([
        row("ITS", "Aspergillus niger:5", rank="Aspergillus niger:99.9; Aspergillus awamori:97"),
        row("Btub", "Aspergillus niger:3", rank="Aspergillus niger:99.7; Aspergillus tubingensis:99.8; Aspergillus x:97"),
    ])
    assert r["判讀狀態"] == STATUS_UNIQUE      # 差 0.1 <= 0.2；ITS 差距 2.9 >= 0.8


def test_unique_needs_gap_to_second_species():
    r = call([
        row("ITS", "Aspergillus niger:5", rank="Aspergillus niger:99.8; Aspergillus awamori:99.5"),
        row("Btub", "Aspergillus niger:3", rank="Aspergillus niger:99.6; Aspergillus tubingensis:99.2"),
    ])
    assert r["判讀狀態"] == STATUS_INTERSECT_PENDING
    assert "差距均 < 0.8%" in r["待確認原因"]


def test_gap_from_one_discriminating_marker_is_enough():
    # ITS 分不開 flavus/oryzae（差 0），CaM 差 3% -> 整體有鑑別力
    r = call([
        row("ITS", "Aspergillus flavus:5; Aspergillus oryzae:5",
            rank="Aspergillus flavus:100; Aspergillus oryzae:100"),
        row("CF", "Aspergillus flavus:4", rank="Aspergillus flavus:100; Aspergillus oryzae:97"),
    ])
    assert r["判讀狀態"] == STATUS_UNIQUE


def test_synonym_is_not_counted_as_runner_up():
    # Hypocrea lixii 與 Trichoderma lixii 同種小名不同屬 -> 視為同義，不算次佳競爭種
    r = call([
        row("TEF1", "Trichoderma lixii:5", genus="Trichoderma",
            rank="Trichoderma lixii:100; Hypocrea lixii:100; Trichoderma atroviride:97"),
        row("RPB2", "Trichoderma lixii:5", genus="Trichoderma",
            rank="Trichoderma lixii:100; Trichoderma atroviride:98"),
    ])
    assert r["判讀狀態"] == STATUS_UNIQUE


def test_missing_rank_data_cannot_be_green():
    rows_ = [row("Btub", "Aspergillus clavatus:5"), row("CF", "Aspergillus clavatus:5")]
    rows_[0]["species_best_pident_all"] = ""
    r = call(rows_)
    assert r["判讀狀態"] == STATUS_INTERSECT_PENDING
    assert "無法驗證" in r["待確認原因"]


def test_gate_and_rank_reasons_combine_into_pending_confirm():
    r = call([
        row("Btub", "Aspergillus niger:3", asm="WARN_CONSENSUS_CHECK",
            rank="Aspergillus niger:99.0; Aspergillus awamori:99.0"),
        row("CF", "Aspergillus niger:3", rank="Aspergillus niger:99.0; Aspergillus awamori:99.0"),
    ])
    assert r["判讀狀態"] == STATUS_UNIQUE_PENDING
    assert "組裝" in r["待確認原因"] and "差距" in r["待確認原因"]


def test_thresholds_come_from_config():
    c = load_config()["consensus"]
    assert c["top_tolerance_pct"] == 0.2 and c["min_gap_pct"] == 0.8
    cfg_text = (open(__import__("pipeline.config", fromlist=["x"]).DEFAULT_CONFIG_PATH, encoding="utf-8").read())
    assert "內部暫定" in cfg_text


# ---------- 🟡-14：同一株跨訂單合併判讀 ----------

def _tax_csv(path, rows_):
    import csv as _csv
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = ["sample", "marker", "assembly_status", "candidate_genus", "high_quality_species_all",
            "species_best_pident_all", "review_flag"]
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        w = _csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows_)


def _tax_row(sample, marker, hq, rank, asm="PASS"):
    return {"sample": sample, "marker": marker, "assembly_status": asm, "candidate_genus": "Aspergillus",
            "high_quality_species_all": hq, "species_best_pident_all": rank, "review_flag": "OK"}


def test_same_sample_across_orders_is_merged(tmp_path):
    from pipeline.consensus_call import build_rows

    _tax_csv(tmp_path / "S1_ncbi_blast" / "taxonomy_decision_table.csv", [
        _tax_row("63aSID322a", "ITS", "Aspergillus welwitschiae:5; Aspergillus niger:5",
                 "Aspergillus welwitschiae:100; Aspergillus niger:100"),
        _tax_row("other", "ITS", "Aspergillus flavus:3", "Aspergillus flavus:100"),
    ])
    _tax_csv(tmp_path / "S2_ncbi_blast" / "taxonomy_decision_table.csv", [
        _tax_row("63aSID322a", "Btub", "Aspergillus welwitschiae:4", "Aspergillus welwitschiae:100; Aspergillus niger:97"),
    ])
    rows_ = build_rows(tmp_path, ["S2"], CFG)           # 只處理 S2，S1 是歷史單
    assert len(rows_) == 1
    r = rows_[0]
    assert r["樣本"] == "63aSID322a" and r["訂單"] == "S2"
    assert r["判讀狀態"] == STATUS_UNIQUE and r["唯一物種"] == "Aspergillus welwitschiae"
    assert r["合併訂單"] == "S1、S2" and "跨訂單合併判讀" in r["備註"]
    assert r["Marker 數"] == "2"
    # 沒有補做新 marker 的歷史樣本 "other" 不會被拉進來
    both = build_rows(tmp_path, ["S1", "S2"], CFG)
    assert [x["樣本"] for x in both] == ["63aSID322a", "other"]
    assert [x["訂單"] for x in both] == ["S1", "S1"]


def test_merge_can_be_switched_off(tmp_path):
    from pipeline.consensus_call import build_rows

    _tax_csv(tmp_path / "S1_ncbi_blast" / "taxonomy_decision_table.csv",
             [_tax_row("A1", "ITS", "Aspergillus niger:5", "Aspergillus niger:100")])
    _tax_csv(tmp_path / "S2_ncbi_blast" / "taxonomy_decision_table.csv",
             [_tax_row("A1", "Btub", "Aspergillus niger:5", "Aspergillus niger:100")])
    cfg = {**CFG, "consensus": {**CFG["consensus"], "merge_across_orders": False}}
    assert build_rows(tmp_path, ["S2"], cfg)[0]["判讀狀態"] == STATUS_SINGLE
    assert build_rows(tmp_path, ["S2"], CFG)[0]["Marker 數"] == "2"


# ---------- 🟡-3：沒有唯一物種就給建議，排除已做過的 marker ----------

def test_suggestion_given_when_its_has_many_species():
    r = call([row("ITS", "Penicillium rubens:5; Penicillium chrysogenum:5; Penicillium allii:3", genus="Penicillium")])
    assert r["判讀狀態"] == STATUS_GENUS_ONLY
    assert r["下次建議 marker"].startswith("Btub")


def test_suggestion_when_intersection_multiple_excludes_done_markers():
    r = call([
        row("Btub", "Penicillium rubens:5; Penicillium chrysogenum:5", genus="Penicillium"),
        row("CF", "Penicillium rubens:5; Penicillium chrysogenum:5", genus="Penicillium"),
    ])
    assert r["判讀狀態"] == STATUS_MULTI
    assert "Btub" not in r["下次建議 marker"] and "CF" not in r["下次建議 marker"]
    assert r["下次建議 marker"].startswith("RPB2")


def test_suggestion_when_intersection_empty():
    r = call([
        row("ITS", "Aspergillus flavus:5", genus="Aspergillus"),
        row("Btub", "Aspergillus oryzae:5", genus="Aspergillus"),
    ])
    assert r["判讀狀態"] == STATUS_NONE
    assert r["下次建議 marker"].startswith("CF")          # Aspergillus 第一優先 CF，ITS／Btub 已做過


def test_no_suggestion_when_green_or_pending_assembly():
    green = call([row("Btub", "Aspergillus clavatus:5"), row("CF", "Aspergillus clavatus:5")])
    assert green["判讀狀態"] == STATUS_UNIQUE and green["下次建議 marker"] == ""
    pend = call([row("Btub", "Aspergillus clavatus:5", asm="WARN_CONSENSUS_CHECK"),
                 row("CF", "Aspergillus clavatus:5")])
    assert pend["判讀狀態"] == STATUS_UNIQUE_PENDING and pend["下次建議 marker"] == ""


def test_suggestion_for_intersect_pending_gap_failure():
    r = call([
        row("ITS", "Aspergillus niger:5", rank="Aspergillus niger:99.8; Aspergillus awamori:99.5"),
        row("Btub", "Aspergillus niger:3", rank="Aspergillus niger:99.6; Aspergillus tubingensis:99.2"),
    ])
    assert r["判讀狀態"] == STATUS_INTERSECT_PENDING
    assert r["下次建議 marker"].startswith("CF")


def test_not_blasted_marker_counts_as_done_for_suggestion():
    from pipeline.consensus_call import call_sample

    r = call_sample("S1", "X1", [row("ITS", "Penicillium rubens:5; Penicillium allii:3", genus="Penicillium")],
                    EXCL, LOOKUP, unrun={"Btub": "尚未送 BLAST（NOT_BLASTED）"})
    assert not r["下次建議 marker"].startswith("Btub")


def test_unverified_primers_never_print_none():
    # 知識庫引子為 null（Neurospora 的多基因座）-> 不可輸出 "None None"，也不能被當成可下單
    r = call([row("ITS", "Neurospora crassa:5", genus="Neurospora")])
    assert "None" not in r["建議引子"]
    assert "引子待查，不可下單" in r["建議引子"]
