"""引子名 -> marker 映射測試：涵蓋 gene_lookup 全部引子名。"""

import json

import pytest

from pipeline.markers import (
    LOOKUP_PATH,
    UNKNOWN_MARKER,
    is_forward,
    load_primer_table,
    lookup_primer_names,
    marker_for_primer,
)
from pipeline.sanger_pair_assemble import primer_order
from pipeline.sanger_qc import parse_name


@pytest.mark.parametrize(
    "primer,marker",
    [
        ("CMD5", "CF"), ("cmd6", "CF"), ("CF1", "CF"), ("CF4", "CF"), ("CL1", "CF"), ("CL2A", "CF"),
        ("5F", "RPB2"), ("5Feur", "RPB2"), ("5f2", "RPB2"), ("7cf", "RPB2"), ("7CR", "RPB2"),
        ("7CReur", "RPB2"), ("11ar", "RPB2"), ("rpb2AM-7R", "RPB2"),
        ("gpd1", "GAPDH"), ("GPD2", "GAPDH"),
        ("ACT-512F", "ACT"), ("act-783r", "ACT"),
        ("Bt2a", "Btub"), ("Bt2b", "Btub"), ("T1", "Btub"), ("T10", "Btub"), ("T22", "Btub"),
        ("Ben2f", "Btub"), ("TUB4Rd", "Btub"),
        ("EF1-728F", "TEF1"), ("EF-1", "TEF1"), ("TEF1LLErev", "TEF1"),
        ("NL-1", "LSU"), ("LR0R", "LSU"), ("LR5", "LSU"),
        ("ITS1", "ITS"), ("ITS4", "ITS"), ("ITS5", "ITS"), ("V9G", "ITS"), ("LS266", "ITS"),
    ],
)
def test_known_primers(primer, marker):
    assert marker_for_primer(primer) == marker


def test_every_lookup_primer_is_mapped():
    """gene_lookup 內每個引子名都要有對應 marker（不得掉進「其他」）。"""
    data = json.loads(LOOKUP_PATH.read_text(encoding="utf-8"))
    assert data
    names = lookup_primer_names()
    assert len(names) >= 40
    for name in names:
        assert marker_for_primer(name) != UNKNOWN_MARKER, name
        # 大小寫不敏感
        assert marker_for_primer(name.upper()) == marker_for_primer(name.lower())


def test_unknown_primer_goes_to_other_not_dropped():
    assert marker_for_primer("M13F") == UNKNOWN_MARKER
    assert marker_for_primer("") == UNKNOWN_MARKER
    meta = parse_name("S1_sampleA_M13F_01")
    assert meta["marker"] == UNKNOWN_MARKER
    assert meta["primer"] == "M13F"


def test_parse_name_new_primers():
    assert parse_name("S1_x_CMD5_01")["marker"] == "CF"
    assert parse_name("S1_x_7CR_01")["marker"] == "RPB2"
    assert parse_name("S1_x_gpd1_01")["marker"] == "GAPDH"


def test_legacy_markers_unchanged():
    assert parse_name("S1_a_ITS5_01")["marker"] == "ITS"
    assert parse_name("S1_a_Bt2b_01")["marker"] == "Btub"
    assert parse_name("S1_a_EF-1_01")["marker"] == "TEF1"
    assert parse_name("S1_a_CF4_01")["marker"] == "CF"


def test_primer_order_forward_first():
    for fwd in ["ITS5", "Bt2a", "EF-1", "CF1", "CMD5", "5F", "Ben2f"]:
        assert is_forward(fwd), fwd
    for rev in ["ITS4", "Bt2b", "CF4", "CMD6", "7CR", "EF2"]:
        assert not is_forward(rev), rev
    assert primer_order("CMD5") < primer_order("CMD6")
    assert primer_order("5F") < primer_order("7CR")
    assert primer_order("ITS5") < primer_order("ITS4")


def test_table_has_roles():
    table = load_primer_table()
    assert table["cmd5"]["role"] == "fwd" and table["cmd6"]["role"] == "rev"
