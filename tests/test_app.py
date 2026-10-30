"""Streamlit AppTest 煙霧測試 + 上傳驗證單元測試。

需要 NAS 的測試（回歸批次 S2604280595 離線整流程）在 NAS 沒掛載時自動跳過。
NAS 只讀；所有輸出都進專案 test_output/。
"""

import io
import time
import zipfile
from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

from pipeline.config import PROJECT_ROOT, load_config, rawdata_dir
from webapp import runner

APP = str(PROJECT_ROOT / "app.py")
PW = "pytest-pw-123"
CFG = load_config()
NAS_RAW = rawdata_dir(CFG)
NAS_OUT = NAS_RAW.parent / CFG["paths"]["output_subdir"]
ORDER = "S2604280595"
HAS_NAS = (NAS_RAW / f"{ORDER}_S16133").is_dir() and (NAS_OUT / f"{ORDER}_ncbi_blast" / "raw_blast_xml").is_dir()


def make_app(timeout=30):
    at = AppTest.from_file(APP, default_timeout=timeout)
    at.secrets["APP_PASSWORD"] = PW
    return at


def fake_zip(names, folder=""):
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        for n in names:
            zf.writestr(f"{folder}{n}", b"x")
    return buf.getvalue()


# ---------------------------------------------------------------- 密碼閘
def test_gate_blocks_without_password():
    at = make_app().run()
    assert not at.exception
    assert len(at.text_input) == 1 and at.text_input[0].label == "請輸入密碼"
    assert len(at.file_uploader) == 0          # 沒登入看不到上傳區


def test_gate_wrong_password():
    at = make_app().run()
    at.text_input[0].set_value("wrong").run()
    assert len(at.error) == 1 and "密碼不對" in at.error[0].value
    assert len(at.file_uploader) == 0


def test_gate_correct_password_shows_uploader():
    at = make_app().run()
    at.text_input[0].set_value(PW).run()
    assert not at.exception
    assert len(at.file_uploader) == 1


def test_gate_locked_when_secret_missing():
    at = AppTest.from_file(APP, default_timeout=30)
    at.secrets["APP_PASSWORD"] = ""
    at.run()
    assert "已鎖定" in at.error[0].value
    assert len(at.text_input) == 0 and len(at.file_uploader) == 0


# ---------------------------------------------------------------- 上傳驗證
def test_bad_filename_message_lists_file(tmp_path):
    data = fake_zip(["S1_sampleA_ITS4_01.ab1", "badname.ab1", "S1_x__01.ab1"], folder="S1_B1/")
    with pytest.raises(runner.UploadError) as exc:
        runner.extract_uploads([("res.zip", data)], tmp_path)
    msg = str(exc.value)
    assert "badname.ab1" in msg and "S1_x__01.ab1" in msg
    assert "S1_sampleA_ITS4_01.ab1" not in msg.split("正確格式")[1].split("\n", 1)[1]


def test_zip_without_ab1(tmp_path):
    with pytest.raises(runner.UploadError, match="找不到任何 .ab1"):
        runner.extract_uploads([("a.zip", fake_zip(["readme.txt"]))], tmp_path)


def test_not_a_zip(tmp_path):
    with pytest.raises(runner.UploadError, match="不是有效的 zip"):
        runner.extract_uploads([("a.zip", b"not zip")], tmp_path)


def test_zip_slip_ignored(tmp_path):
    data = fake_zip(["../../evil_x_ITS4_01.ab1", "S1_A_ITS4_01.ab1"], folder="S1_B/")
    rep = runner.extract_uploads([("a.zip", data)], tmp_path / "raw")
    assert rep.n_files == 1
    assert not (tmp_path / "evil_x_ITS4_01.ab1").exists()


def test_report_groups_and_unknown_primer(tmp_path):
    data = fake_zip(["S1A_x_ITS4_01.ab1", "S1B_x_M13F_02.ab1", "._junk.ab1"], folder="S1_B1/")
    rep = runner.extract_uploads([("a.zip", data)], tmp_path)
    assert rep.orders == ["S1"] and rep.n_files == 2
    assert rep.unknown_primers == {"M13F": 1}
    assert (tmp_path / "S1_B1" / "S1A_x_ITS4_01.ab1").is_file()


def test_output_under_nas_is_refused(tmp_path):
    cfg = {**CFG, "paths": {**CFG["paths"], "nas_root": str(tmp_path / "nas")}}
    with pytest.raises(runner.UploadError):
        runner.assert_not_under_nas(tmp_path / "nas" / "rawdata" / "x", cfg)
    with pytest.raises(runner.UploadError):
        runner.assert_not_under_nas(Path("/Volumes/whatever/out"), cfg)
    runner.assert_not_under_nas(PROJECT_ROOT / "test_output" / "x", cfg)


def test_output_guard_with_empty_nas_root():
    cfg = {**CFG, "paths": {**CFG["paths"], "nas_root": ""}}
    with pytest.raises(runner.UploadError):
        runner.assert_not_under_nas(Path("/Volumes/whatever/out"), cfg)
    runner.assert_not_under_nas(PROJECT_ROOT / "test_output" / "x", cfg)


def test_format_primer_pair():
    assert runner.format_primer_pair("rpb2: 5F GAYGAY / 7CR CCCATR") == "5F / 7CR"
    assert runner.format_primer_pair("CaM: CMD5 AAA / CMD6 CCC；CaM(備用): CF1 GGG / CF4 TTT") == "CMD5 / CMD6；備用 CF1 / CF4"


# ---------------------------------------------------------------- 整流程（離線，回歸批次）
@pytest.mark.skipif(not HAS_NAS, reason="NAS 未掛載")
def test_full_upload_flow_offline(monkeypatch, tmp_path):
    # 修正 🟡-18：每次用全新臨時目錄，確保 QC／Assembly／BLAST 解析真的重算，不吃上次留下的輸出
    work = tmp_path / "jobs"
    work.mkdir()
    started = time.time() - 1
    monkeypatch.setattr(runner, "WORK_BASE", work)

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        for sub in sorted(NAS_RAW.glob(f"{ORDER}_*")):       # NAS 只讀
            for f in sorted(sub.glob("*.ab1")):
                zf.write(f, f"{sub.name}/{f.name}")
    zip_bytes = buf.getvalue()
    xmls = [(p.name, p.read_bytes()) for p in sorted((NAS_OUT / f"{ORDER}_ncbi_blast" / "raw_blast_xml").glob("*.xml"))]
    assert xmls

    at = make_app(timeout=600).run()
    at.text_input[0].set_value(PW).run(timeout=600)
    at.file_uploader(key="zips").upload("s2604280595.zip", zip_bytes).run(timeout=600)
    assert not at.exception, at.exception
    assert any("收到" in s.value for s in at.success)
    assert at.button(key="go_upload")

    at.checkbox(key="offline_up").check().run(timeout=600)
    at.file_uploader(key="xmls").set_value([(n, d, "text/xml") for n, d in xmls]).run(timeout=600)
    at.button(key="go_upload").click().run(timeout=600)
    assert not at.exception, at.exception
    last = at.session_state["last_run"]
    assert last["failures"] == [], last["failures"]
    assert last["orders"] == [ORDER]

    html_all = " ".join(m.value for m in at.markdown)
    assert "下次建議 Marker＋引子對" in html_all
    assert "51aIMCN3c" in html_all
    assert "RPB2" in html_all and "5F / 7CR" in html_all      # 屬 Trichoderma -> 下次建議 RPB2
    assert len(at.download_button) == 1

    out = Path(last["out"])
    assert str(out).startswith(str(work))                 # 輸出在本次臨時夾內
    assert not str(out).startswith(str(NAS_RAW.parent))   # 絕不在 NAS
    # 證明是真的重算：QC／Assembly 輸出的修改時間都晚於本測試開始時間
    asm = out / f"{ORDER}_combined_assembly" / "assembly_summary.csv"
    assert asm.is_file() and asm.stat().st_mtime >= started
    res = runner.build_result_rows(out, [ORDER])
    assert res["counts"]["samples"] >= 14
    z = zipfile.ZipFile(io.BytesIO(runner.build_result_zip(out, res["rows"])))
    names = z.namelist()
    assert "結果總表.csv" in names
    assert any(n.endswith("paired_consensus.fasta") for n in names)
    assert any(n.endswith("taxonomy_decision_table.csv") for n in names)


def test_local_mode_hidden_by_default_and_shown_by_switch(monkeypatch):
    at = make_app().run()
    at.text_input[0].set_value(PW).run()
    assert len(at.tabs) == 0                       # 預設（雲端版）沒有本機模式分頁
    monkeypatch.setenv("SANGER_ENABLE_LOCAL_MODE", "1")
    at2 = make_app().run()
    at2.text_input[0].set_value(PW).run()
    assert not at2.exception
    assert len(at2.tabs) == 2


# ---------------------------------------------------------------- 結果頁燈號（🔴-2／🔴-3）
def _write(path, header, rows):
    import csv
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()
        w.writerows(rows)


def test_result_levels_green_only_for_clean_unique(tmp_path):
    from pipeline import consensus_call as cc

    def crow(sample, status, unique="", genus="", cand=""):
        r = {h: "" for h in cc.HEADER}
        r.update({"訂單": "S1", "樣本": sample, "判讀狀態": status, "唯一物種": unique,
                  "判讀屬": genus, "共同候選交集": cand, "Marker 數": "2"})
        return r

    _write(tmp_path / "2026-01-01_new_orders_marker_intersection_summary.csv", cc.HEADER, [
        crow("a", cc.STATUS_UNIQUE, unique="Aspergillus niger"),
        crow("b", cc.STATUS_UNIQUE_PENDING, unique="Aspergillus niger"),
        crow("c", cc.STATUS_GENUS_ONLY, genus="Penicillium"),
        crow("d", cc.STATUS_NOT_BLASTED),
    ])
    rows = {r["樣本"]: r for r in runner.build_result_rows(tmp_path, ["S1"])["rows"]}
    assert rows["a"]["等級"] == "ok" and rows["a"]["物種／候選"] == "Aspergillus niger"
    assert rows["b"]["等級"] == "warn" and "待確認" in rows["b"]["物種／候選"]
    assert rows["c"]["等級"] == "warn" and "僅到屬" in rows["c"]["物種／候選"]
    assert "Penicillium chrysogenum" not in rows["c"]["物種／候選"]
    assert rows["d"]["等級"] == "warn" and "尚未比對" in rows["d"]["物種／候選"]
