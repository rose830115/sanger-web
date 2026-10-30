"""Sanger 定序比對工具：網頁介面（Streamlit）。

本機啟動：  .venv/bin/streamlit run app.py
密碼：      .streamlit/secrets.toml 的 APP_PASSWORD（範本見 secrets.toml.example；雲端在 App settings -> Secrets 設定）
"""

from __future__ import annotations

import hmac
import html
from pathlib import Path

import streamlit as st

from pipeline.config import PROJECT_ROOT, load_config, rawdata_dir
from pipeline.scan_batches import scan
from webapp import runner

st.set_page_config(page_title="Sanger 定序比對", page_icon="🧬", layout="wide")

CSS = """
<style>
html, body, [class*="css"] { font-size: 19px; }
h1 { font-size: 2.2rem !important; }
.stButton button, .stDownloadButton button { font-size: 1.25rem; padding: 0.6rem 1.6rem; }
table.res { border-collapse: collapse; width: 100%; font-size: 1.05rem; }
table.res th { background: #2f4f4f; color: #fff; padding: 10px; text-align: left; }
table.res td { padding: 10px; border-bottom: 1px solid #ddd; vertical-align: top; }
table.res tr.ok td.stat { background: #e3f4e3; }
table.res tr.warn td.stat { background: #fff3cd; }
table.res tr.fail td.stat { background: #f8d7da; }
table.res td.next { background: #fff7d6; border: 3px solid #e0a800; font-size: 1.25rem; font-weight: 700; min-width: 240px; }
table.res td.next .pair { display: block; font-size: 1.15rem; color: #7a4b00; margin-top: 4px; }
table.res td.next.none { background: #f3f3f3; border: 1px solid #ddd; color: #888; font-weight: 400; }
.stepbox { background: #f4f8f8; border-left: 6px solid #2f4f4f; padding: 12px 18px; margin: 10px 0; border-radius: 6px; }
</style>
"""
st.markdown(CSS, unsafe_allow_html=True)


# ------------------------------------------------------------------ 密碼閘
def password_ok() -> bool:
    """密碼存在 st.secrets["APP_PASSWORD"]；沒設定就鎖死（寧可打不開，不可裸奔）。"""
    if st.session_state.get("authed"):
        return True
    try:
        secret = st.secrets["APP_PASSWORD"]
    except Exception:  # noqa: BLE001  (沒有 secrets 檔或沒這個 key)
        secret = None
    st.title("🧬 Sanger 定序比對")
    if not secret:
        st.error("網站尚未設定密碼（APP_PASSWORD），為了安全已鎖定。請管理員先完成設定。")
        return False
    entered = st.text_input("請輸入密碼", type="password", key="pw_input")
    if entered:
        if hmac.compare_digest(str(entered).encode(), str(secret).encode()):
            st.session_state["authed"] = True
            st.rerun()
        st.error("密碼不對，請再試一次。")
    return False


# ------------------------------------------------------------------ 結果頁
def render_table(rows: list[dict[str, str]]) -> None:
    e = html.escape
    head = "".join(
        f"<th>{h}</th>" for h in ["樣本", "Marker", "判讀狀態", "物種／候選", "下次建議 Marker＋引子對"]
    )
    body = []
    for r in rows:
        has_next = r["下次建議 Marker"] and r["下次建議 Marker"] != runner.NO_SUGGESTION
        nxt = e(r["下次建議 Marker"])
        if has_next and r["建議引子對"]:
            nxt += f'<span class="pair">引子：{e(r["建議引子對"])}</span>'
        cls = "next" if has_next else "next none"
        note = f'<br><small>{e(r["備註"])}</small>' if r["備註"] else ""
        body.append(
            f'<tr class="{r["等級"]}"><td><b>{e(r["樣本"])}</b><br><small>{e(r["訂單"])}</small></td>'
            f'<td>{e(r["Marker"])}</td><td class="stat">{e(r["判讀狀態"])}{note}</td>'
            f'<td>{e(r["物種／候選"])}</td><td class="{cls}">{nxt}</td></tr>'
        )
    st.markdown(
        f'<table class="res"><thead><tr>{head}</tr></thead><tbody>{"".join(body)}</tbody></table>',
        unsafe_allow_html=True,
    )


def show_results(out: Path, orders: list[str], failures: list[str], stamp: str) -> None:
    res = runner.build_result_rows(out, orders)
    c = res["counts"]
    st.header("② 結果")
    a, b, d, f = st.columns(4)
    a.metric("樣本數", c["samples"])
    b.metric("已定到種", c["ok"])
    d.metric("需注意", c["warn"])
    f.metric("未能比對", c["fail"])
    st.caption(
        "「已定到種」＝兩個以上 marker 交集唯一、各 marker 組裝都 PASS、沒有阻斷旗標。"
        "組裝有疑慮或旗標未過的唯一物種會標「待確認」、最高分／差距條件未過的標「待覆核」（黃），請先看色譜圖與 BLAST 排名；"
        "ITS 分不開種的屬只報到屬／section。"
    )
    for w in res["warnings"]:
        st.warning(w)
    if failures:
        st.error("以下步驟失敗（已完成的部分仍可看）：\n\n" + "\n".join(f"- {x}" for x in failures))
    st.subheader("黃框欄位＝下次要送定序的 Marker 與引子，回填 Notion「下次送定序 Marker」用")
    if res["rows"]:
        render_table(res["rows"])
        with st.expander("引子完整序列"):
            for r in res["rows"]:
                if r["引子序列"]:
                    st.text(f"{r['樣本']}：{r['引子序列']}")
    else:
        st.info("沒有可顯示的結果。")
    st.download_button(
        "⬇️ 下載全部結果（zip：結果總表＋所有 CSV／FASTA）",
        data=runner.build_result_zip(out, res["rows"]),
        file_name=f"sanger_results_{stamp}.zip",
        mime="application/zip",
        type="primary",
    )


def store_run(out: Path, orders: list[str], result: dict) -> None:
    st.session_state["last_run"] = {
        "out": str(out),
        "orders": orders,
        "failures": result.get("failures", []),
    }


def run_with_progress(cfg, raw, out, orders, offline) -> None:
    bar = st.progress(0.0, text="準備中…")

    def on_progress(frac: float, msg: str) -> None:
        bar.progress(frac, text=msg)

    result = runner.run_job(cfg, raw, out, orders, offline, on_progress=on_progress)
    store_run(out, orders, result)


# ------------------------------------------------------------------ 模式：上傳
def secret_email(cfg: dict) -> str:
    try:
        return str(st.secrets.get("NCBI_EMAIL", "")) or cfg["blast"].get("email", "")
    except Exception:  # noqa: BLE001
        return cfg["blast"].get("email", "")


def upload_mode(cfg: dict) -> None:
    st.markdown(
        '<div class="stepbox"><b>步驟 1</b>　把定序公司寄來的結果 zip 拖進下面（可一次放好幾個）<br>'
        "<b>步驟 2</b>　按「開始比對」，等它跑完<br><b>步驟 3</b>　看結果、下載</div>",
        unsafe_allow_html=True,
    )
    files = st.file_uploader("定序結果 zip", type=["zip"], accept_multiple_files=True, key="zips")
    if not files:
        return
    payload = [(f.name, f.getvalue()) for f in files]
    job_id = runner.job_id_for(payload)
    job_dir = runner.new_job_dir(job_id)
    prep = st.session_state.setdefault("prep", {})
    if job_id not in prep or not (job_dir / "raw").is_dir():
        try:
            prep[job_id] = runner.extract_uploads(payload, job_dir / "raw")
        except runner.UploadError as exc:
            st.error(str(exc))
            return
    report = prep[job_id]
    st.success(f"收到 {report.n_files} 個定序檔，訂單：{'、'.join(report.orders)}")
    if report.unknown_primers:
        names = "、".join(f"{k}（{v} 檔）" for k, v in report.unknown_primers.items())
        st.warning(f"有程式不認得的引子名：{names}。這些資料不會被丟掉，marker 會標成「其他」，結果請人工確認。")

    with st.expander("進階選項"):
        if st.button("清除重跑（刪掉這批的暫存結果，從頭算）", key="clear_job"):
            runner.clear_job(job_id)
            st.session_state.get("prep", {}).pop(job_id, None)
            st.session_state.pop("last_run", None)
            st.success("已清除這批的暫存結果，請重新按「開始比對」。")
            st.stop()
        offline = st.checkbox(
            "離線重解析（不連 NCBI，只用已有的 BLAST 原始檔 .xml 重新判讀）", key="offline_up"
        )
        xmls = (
            st.file_uploader("BLAST 原始檔（.xml，可多選；僅單一訂單時可用）", type=["xml"],
                             accept_multiple_files=True, key="xmls")
            if offline
            else []
        )
    minutes = runner.estimate_minutes(report, offline)
    if offline:
        st.info("離線模式：不會連線 NCBI，約 1 分鐘內完成。")
    else:
        st.info(
            f"NCBI 有限速，這批大約需要 **{minutes} 分鐘**（每個 marker 約 {runner.MINUTES_PER_MARKER_ONLINE} 分鐘）。"
            "**請不要關閉或重新整理網頁。** 萬一中斷，重新上傳同樣的 zip 會從斷掉的地方接著跑。"
        )
    if st.button("開始比對", type="primary", key="go_upload"):
        run_cfg = {**cfg, "blast": {**cfg["blast"], "email": secret_email(cfg)}}
        if not offline and not run_cfg["blast"]["email"]:
            st.error("管理員尚未設定 NCBI 聯絡信箱（secrets 的 NCBI_EMAIL），目前只能用離線重解析。")
            return
        if offline and xmls:
            if len(report.orders) != 1:
                st.error("上傳 BLAST 原始檔時，一次只能處理單一訂單。")
                return
            runner.seed_blast_xml(job_dir / "out", report.orders, [(x.name, x.getvalue()) for x in xmls])
        run_with_progress(run_cfg, job_dir / "raw", job_dir / "out", report.orders, offline)


# ------------------------------------------------------------------ 模式：本機資料夾
def local_mode(cfg: dict) -> None:
    st.markdown(
        '<div class="stepbox">直接讀實驗室 NAS 的 rawdata（<b>只讀、不會寫入 NAS</b>），結果存到本機資料夾。</div>',
        unsafe_allow_html=True,
    )
    raw_text = st.text_input("rawdata 資料夾", value=str(rawdata_dir(cfg)), key="local_raw")
    out_text = st.text_input(
        "結果輸出夾（不可在 NAS 底下）",
        value=str(PROJECT_ROOT / runner.app_settings(cfg)["local_output_dir"]),
        key="local_out",
    )
    seed_text = st.text_input("既有結果夾（唯讀，離線重解析用，可留空）", value="", key="local_seed")
    raw = Path(raw_text).expanduser()
    if not raw.is_dir():
        st.warning("找不到 rawdata 資料夾（NAS 是否已掛載？）")
        return
    info = scan(raw, Path(out_text).expanduser())
    names = sorted(info["orders"])
    picked = st.multiselect("要處理的訂單", names, key="local_orders")
    offline = st.checkbox("離線重解析（不連 NCBI）", value=True, key="offline_local")
    if st.button("開始比對", type="primary", key="go_local") and picked:
        out = Path(out_text).expanduser()
        try:
            runner.assert_not_under_nas(out, cfg)
        except runner.UploadError as exc:
            st.error(str(exc))
            return
        if seed_text.strip():
            runner.copy_blast_xml_from(Path(seed_text).expanduser(), out, picked)
        run_cfg = {**cfg, "blast": {**cfg["blast"], "email": secret_email(cfg)}}
        run_with_progress(run_cfg, raw, out, picked, offline)


# ------------------------------------------------------------------ 主畫面
def main() -> None:
    if not password_ok():
        return
    cfg = load_config()
    runner.cleanup_old_jobs()        # 自動清掉 24 小時前的工作夾（README 的承諾）
    st.title("🧬 Sanger 定序比對")
    local_on = bool(runner.app_settings(cfg)["enable_local_mode"])
    st.header("① 上傳與比對")
    if local_on:
        tab_up, tab_local = st.tabs(["上傳 zip（一般用這個）", "本機資料夾（實驗室電腦）"])
        with tab_up:
            upload_mode(cfg)
        with tab_local:
            local_mode(cfg)
    else:
        upload_mode(cfg)
    last = st.session_state.get("last_run")
    if last:
        show_results(Path(last["out"]), last["orders"], last["failures"], last["orders"][0] if last["orders"] else "run")


main()
