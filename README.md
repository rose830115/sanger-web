# sanger-web：Sanger 定序比對流程（Phase 2 腳本化）

把 `/sanger-pipeline` 七步流程寫成不需要 AI 的腳本（規則 SSOT：`000_Agent/skills/sanger-pipeline/SKILL.md`）。
Phase 3 在這之上加了 Streamlit 網頁介面（`app.py`），所以每一步都是「可 import 的函式 + 命令列入口」。

## 安裝

```bash
cd 600_Dev/sanger-web
python3 -m venv .venv                 # 本機 Python 3.14.7
.venv/bin/pip install -r requirements.txt
```

## 啟用前必做：填 NCBI 信箱

`config.yaml` 的 `blast.email` **預設留空**（專案內不放任何個人信箱）。
要跑線上 BLAST 前，填實驗室共用信箱（NCBI 要求留聯絡信箱）；也可用環境變數 `SANGER_NCBI_EMAIL`。
沒填就跑線上 BLAST 會直接停下來並提示，不會偷偷送件。

## 設定

所有路徑與判定門檻都在 `config.yaml`（QC、Assembly、BLAST、Step 6/7 參數）。
repo 是公開的，所以 `config.yaml` 裡 `paths.nas_root` 留空、不放實驗室內部路徑。
實驗室本機要用「本機資料夾」模式：把 `config.yaml` 複製成 `config.local.yaml`（已被
`.gitignore` 擋住、不會上傳），填入實際 NAS 掛載路徑，啟動時指定：
`SANGER_CONFIG=config.local.yaml .venv/bin/streamlit run app.py`；
或改用環境變數 `SANGER_NAS_ROOT` 覆蓋。雲端上傳 zip 模式用不到 NAS 路徑。
想把輸出寫到別處（測試）：`paths.output_dir`，或各腳本的 `--output-dir`。
`SANGER_CONFIG` 環境變數或 `--config` 可指定另一份設定檔。

## 七步對照

| SKILL.md 步驟 | 腳本 | 指令（在專案根目錄） |
| :-- | :-- | :-- |
| 1 掃描新批次 | `pipeline/scan_batches.py` | `python3 -m pipeline.scan_batches` |
| 2 QC | `pipeline/sanger_qc.py` | `python3 -m pipeline.sanger_qc RAW_SUBFOLDER -o OUT` |
| 3 Assembly | `pipeline/sanger_pair_assemble.py` | `python3 -m pipeline.sanger_pair_assemble DIR1 DIR2 -o OUT` |
| 4 BLAST | `pipeline/sanger_ncbi_blast.py` | `python3 -m pipeline.sanger_ncbi_blast FASTA -o OUT --assembly-summary CSV` |
| 5 結果摘要 | `pipeline/summarize.py` | `python3 -m pipeline.summarize --orders S...` |
| 6 多 marker 交集判讀 | `pipeline/consensus_call.py` | `python3 -m pipeline.consensus_call --orders S...` |
| 7 跨批次總表 | `pipeline/cross_batch_summary.py` | `python3 -m pipeline.cross_batch_summary` |
| 全部 | `run_all.py` | `python3 run_all.py [--from-step blast] [--orders S...] [--offline] [--dry-run]` |

`run_all.py` 每步「缺什麼補什麼」：既有輸出不覆蓋（`--force` 才重做 QC／Assembly），BLAST 有 XML 斷點續跑，
`--offline` 只重解析既有 XML、絕不連線 NCBI。從 summarize 之後續跑時要用 `--orders` 指定訂單。

## 交集判讀補充（consensus_call）

規則照 SKILL.md Step 6。實作時 SKILL 沒寫到、由腳本補的細節（請 Rose 確認）：

- 「二名法學名」判定：`Genus epithet`，epithet 為小寫；`sp.`／`cf.`／`aff.`、含 clone/gene/genes/genomic 的名稱都不算，只進備註。
- 某個 marker 的高品質候選為空 -> 該 marker 仍計入 Marker 數，交集因此為 0（輸出「無共同高品質物種」）。
- 只有一組 marker 時，「共同候選交集」欄填該 marker 自己的候選（與 2026-06-15 既有成品一致）。
- 「只到屬」定義：所有 marker 都沒有任何二名法候選、但知道屬（taxonomy 表 candidate_genus）。此時查 `data/gene_lookup.json`，
  輸出「下次建議 marker」（優先序最高且該樣本尚未做過的 marker，備用引子一併列出）與「建議引子」。
  屬不在 gene_lookup 時輸出「需人工查證」，不猜。

## 網頁介面（Phase 3）

```bash
cd 600_Dev/sanger-web
.venv/bin/pip install -r requirements.txt
cp .streamlit/secrets.toml.example .streamlit/secrets.toml   # 然後把 APP_PASSWORD 改成實驗室密碼
.venv/bin/streamlit run app.py                                # 瀏覽器開 http://localhost:8501
```

- 一般用法：登入 -> 把定序公司的結果 zip 拖進去 -> 按「開始比對」-> 看結果 -> 下載 zip。
- 結果頁黃框那欄「下次建議 Marker＋引子對」就是回填 Notion「下次送定序 Marker」用的。
- 「進階選項 -> 離線重解析」：不連 NCBI，只用已有的 BLAST 原始檔（.xml）重新判讀。
- 「本機資料夾」分頁（給實驗室電腦用，雲端版看不到）：`config.yaml` 的 `app.enable_local_mode: true`，
  或啟動前 `SANGER_ENABLE_LOCAL_MODE=1`。直接讀 NAS rawdata（只讀），結果寫本機 `local_output/`；輸出夾若在 NAS 底下程式會直接擋掉。
- 檔名格式：定序公司原始檔名 `流水號_樣本名_引子_孔位.ab1`（例：`S260428059517_63aSID322a_ITS4_01.ab1`），不合格會列出是哪個檔。
- 程式碼結構：`app.py`（畫面）、`webapp/runner.py`（驗證／執行／整理結果，不依賴 streamlit）、`pipeline/markers.py`（引子→marker）。

### 部署到 Streamlit Community Cloud（逐步教學）

> 這份教學由小玫瑰整理，還沒有實際部署過；第一次做請照順序，卡住就截圖給小玫瑰。

**前置：** 一個 GitHub 帳號、一個 Streamlit Community Cloud 帳號（用 GitHub 登入 https://share.streamlit.io 即可，免費）。

1. **把程式放上 GitHub。** ✅ 已由小玫瑰完成（2026-10-30）：私有 repo `rose830115/sanger-web`。
   - 上雲內容＝**乾淨快照分支**（本機 `cloud-main`，推到 GitHub 當 `main`）：只含目前程式狀態，**不帶本機開發歷史**（開發歷史中曾短暫出現個人信箱，已清掉但 git 歷史會留底，所以雲端一律用快照）。
   - 上傳前確認過：GitHub 上看不到 `.streamlit/secrets.toml`（被 `.gitignore` 擋住，密碼不會上傳）。
2. **建立 App。** 登入 https://share.streamlit.io -> 右上角 **Create app** -> **Deploy a public app from GitHub**（Private repo 第一次要按授權讓 Streamlit 讀取）。
   - Repository：選 `你的帳號/sanger-web`
   - Branch：`main`
   - Main file path：`app.py`
   - 展開 **Advanced settings**：Python version 選 **3.13**（或 3.12；本機開發用 3.14，雲端若沒有 3.14 就選最新的）。
3. **設定密碼（最重要）。** 同一個 Advanced settings 的 **Secrets** 欄貼上（密碼自己決定、告訴實驗室同事）：
   ```toml
   APP_PASSWORD = "你決定的實驗室密碼"
   NCBI_EMAIL = "實驗室共用信箱"
   ```
   `NCBI_EMAIL` 是 NCBI 要求留的聯絡信箱，請填實驗室共用信箱，不要填個人信箱。之後要改：App 右下角 Manage app -> Settings -> Secrets。
4. 按 **Deploy**，等 2–5 分鐘安裝套件。成功後會得到一個網址（`https://xxx.streamlit.app`），把網址和密碼給同事。
5. **驗收：** 開網址 -> 應該先看到「請輸入密碼」；輸入錯的密碼要被擋；輸入對的才看到上傳區。

**Community Cloud 的限制與對策**

| 限制 | 會發生什麼 | 對策 |
| :-- | :-- | :-- |
| 記憶體約 1GB | 一次上傳太多 zip 可能被強制重啟 | 一次放 1–2 個訂單；超過解壓後 400MB 程式會直接拒收 |
| 閒置約 12 小時休眠 | 隔天第一次開要按「Wake up」等 30 秒 | 正常現象，按鈕就好 |
| 臨時檔不保證保留 | 重開機後上次的結果會消失 | **跑完立刻按「下載全部結果」**；程式會自動清掉 24 小時前的舊工作夾 |
| NCBI 限速 | 線上 BLAST 每個 marker 要送兩次（全庫＋type material 第二輪，P5b），約 8 分鐘；另加學名查詢（E-utilities，每秒 ≤3 次） | 畫面會預估時間；請別關網頁。中斷後重傳同樣 zip 會從 BLAST 已存的進度續跑 |
| 雲端碰不到 NAS | 無法直接讀實驗室 NAS | 一律用「上傳 zip」；「本機資料夾」分頁只在實驗室電腦開 |

**更新程式：** 在本機 `main` 改好並驗收後，由小玫瑰做快照同步（不要直接 push 本機 main，開發歷史不上雲）：

```bash
tree=$(git rev-parse main^{tree})
prev=$(git rev-parse cloud-main)
new=$(git commit-tree "$tree" -p "$prev" -m "release: sync main@$(git rev-parse --short main)")
git branch -f cloud-main "$new"
git push origin cloud-main:main
```

Streamlit 偵測到 GitHub `main` 變動會自動重新部署（沒吃到新版就 Manage app -> Reboot）。

## 測試

```bash
.venv/bin/python -m pytest tests/ -v
.venv/bin/python -m pyflakes pipeline/*.py webapp/*.py *.py tests/*.py
```

## BLAST 層（P5b，2026-10-30）

對應盲審統整報告 🔴-5、🔴-8、🟡-1、🟢-1。所有新數值都在 `config.yaml`，並註明「內部暫定」。

| 項目 | 做法 | 設定 |
| :-- | :-- | :-- |
| type material 第二輪（🔴-5） | 每個 marker 除了全庫 nt，再送一次 `entrez_query="sequence from type[filter]"`。type 有高品質候選 -> type 主判（`candidate_source=type`），全庫結果另存 `full_*` 欄當輔助。第二輪沒有高品質 hit（`NO_TYPE_SUPPORT`）或沒跑成（`TYPE_NOT_RUN`）-> 沿用全庫候選，**只降信心**：唯一物種降為「交集唯一（待覆核）」，不判錯。第二輪網路失敗不會中斷整單。每個 BLAST request 有逾時保護：超過 `blast.request_timeout_sec`（預設 900 秒）自動放棄並同參數重發一次（log 記 `BLAST_RESENT`，網頁顯示「NCBI 佇列壅塞，已自動重送」）；重發後再逾時 `blast.resend_timeout_sec`（預設 1800 秒）：全庫輪該 marker 批次失敗（可斷點續跑），type 輪只降信心並在 `type_round_errors` 標 `TIMEOUT_ABANDONED`（與網路失敗區分）。 | `blast.type_round`、`blast.type_entrez_query` |
| 清單加大（🔴-8） | `hitlist_size` 100 -> 1000，解析時依「物種」去重（同一物種多筆只算一個候選）。NCBI 的 URL API 沒有依 taxid 去重的選項、XML 也不含 taxid，所以不能在伺服器端去重，採「取大清單＋解析時去重」。BLAST 會把序列完全相同的紀錄併成一筆 hit，這些併入紀錄各自可能是不同物種，現在都計入候選（`merged_species` 欄）。清單被高品質 hit 塞滿且最後一筆仍達標 -> `HITLIST_SATURATED`（列入阻斷旗標，不給綠燈）。 | `blast.hitlist_size` |
| 學名正規化（🟡-1、🟢-1） | 每筆 hit 的 accession -> taxid（E-utilities esummary）-> NCBI Taxonomy 現行名稱（efetch）。舊名／同義名（Hypocrea lixii、Emericella、Eurotium…）在交集前就換成現行名；品系／亞種歸回種；section、species complex、hybrid、Synthetic construct 不當種名。查不到 taxid 的退回標題猜測（`name_source=title`）。逐筆對照存 `taxonomy_resolution.csv`。 | `blast.taxonomy_lookup`、`taxonomy_cache`、`taxonomy_cache_days`、`taxonomy_min_pident`、`consensus.genus_aliases` |

**taxid 查詢的限速與快取：** 沒填 API key 每秒最多 3 次請求，程式固定間隔 0.4 秒（環境變數 `NCBI_API_KEY` 有值才放寬到 0.12 秒）；
esummary 每批 200 個 accession、taxonomy 每批 100 個 taxid。accession->taxid 與 taxid->名稱落檔在 `cache/taxonomy_cache.json`（不進 git、不含信箱），
名稱超過 90 天才重查。沒有信箱或 `taxonomy_lookup: false` 時只吃快取。只替 identity >= 95%（或排名前 20）的 hit 查。

**輸出變動：** `blast_top100_hits.csv`（名稱沿用，實際最多 1000 筆）多了 `ref_source`、`accession_version`、`taxid`、`taxon_name`、`taxon_rank`、`name_source`、`title_species_guess`、`merged_records`、`merged_species`；
新增 `blast_type_hits.csv`（type 第二輪的 hit）、`taxonomy_resolution.csv`、`raw_blast_xml_type/`；
`taxonomy_decision_table.csv` 多了 `candidate_source`、`type_status`、`full_*`；交集表多了「判讀來源」欄（例：`CF:type；ITS:全庫(無 type 支持)`）。

**已知限制：**
- NCBI Taxonomy 的「現行名稱」不等於 MycoBank／Index Fungorum 的正式名稱；入庫前仍建議對 MycoBank／Index Fungorum（P5 工單原則）。
- type 第二輪只能確認「這個種在 type 序列裡有高分 hit」，type 庫沒收的種（沒有 type 序列的種）會落到「無 type 支持，降信心」。
- 沒有側車檔的舊 XML：離線沿用；線上重跑會改名 `.stale` 並重送（無從得知當時的 hitlist 與篩選條件）。
- 補查學名（舊結果已經用標題猜測解析過）：刪掉該訂單 `_ncbi_blast/taxonomy_decision_table.csv` 再跑 `--from-step blast` 即可重新解析（XML 會逐條驗證，不會重送 NCBI）。

## SKILL.md 與實作的差異清單（🟢-7，供之後人工回寫 SKILL）

SKILL.md 在 `000_Agent/skills/sanger-pipeline/`，本專案不改它。以實作為準的差異如下：

| SKILL.md 寫法 | 實作 |
| :-- | :-- |
| 腳本與輸出放 NAS 同一層，`python3 sanger_xxx.py` 直接跑 | 腳本在 `600_Dev/sanger-web/pipeline/`，用 `python3 -m pipeline.xxx` 或 `run_all.py`；路徑與門檻全在 `config.yaml`，NAS 唯讀 |
| Step 1 以「輸出夾存在」判定做完 | 輸入雜湊清單 `{訂單}_input_manifest.json`：ab1／設定／知識庫一變，下游自動失效，舊輸出搬進 `_stale/` |
| Step 3 狀態只有 PASS／WARN_MANUAL_REVIEW／FAIL_ASSEMBLY_CHECK／FAIL_PAIRING_CHECK | 另有 `WARN_CONSENSUS_CHECK`（identity < 98% 或 N > 10）、`FAIL_READ_UNUSABLE`（壞 read 只擋該對，不炸整單）；QC 另有 `FAIL_READ_UNUSABLE` |
| Step 4 範例指令 `--hitlist-size 100 --email <個人信箱>`、`run_blast.py` | 預設 hitlist 1000；信箱只能填實驗室共用信箱（config／環境變數／網頁 secrets），不寫進檔案；`run_blast.py` 不存在，用 `run_all.py --from-step blast` |
| Step 4「XML 已存在（> 1000 bytes）自動跳過」 | XML 旁有側車檔，逐條比序列 SHA1，對不上的重送（舊 XML 沒側車檔：離線沿用、線上重送）；沒送的樣本標 `NOT_BLASTED`（與 `NO_HITS` 不同） |
| Step 4 只送全庫 nt | 每個 marker 另送 type material 第二輪（P5b） |
| Step 4 覆蓋率未定義 | qcov 分母取 query 與參考序列較短者；最短對齊長度 300 bp（`min_align_len`） |
| Step 5 的 suggested_action 五種 | 另有 `NOT_BLASTED_RERUN_BLAST`、`ITS_GENUS_LEVEL_ONLY`、`TYPE_SUPPORT_MISSING_REVIEW`、`HITLIST_SATURATED_CHECK_TYPE_SEQUENCES`；review_flag 另有 `NOT_BLASTED`、`NO_TYPE_SUPPORT`、`TYPE_NOT_RUN`、`HITLIST_SATURATED` |
| Step 6 取 `high_quality_species_counts`（顯示用，只留前 8 種） | 取 `high_quality_species_all`（完整不截斷）；有 type 結果時取 type 版本 |
| Step 6 判讀狀態四種 | 另有「唯一物種（待確認）」（組裝非 PASS／阻斷旗標）、「交集唯一（待覆核）」（最高分／差距／type 支持未過）、「僅到屬／section，種名待第二 marker」（`its_sufficient=false` 的屬只有 ITS）、「尚未比對（NOT_BLASTED）」；輸出多「待確認原因」「合併訂單」「判讀屬」「判讀來源」欄 |
| Step 6 以 `OrderID + sample` 分組 | 以樣本編號為鍵跨訂單合併（同一株永遠同編號，Rose 2026-10-30 確認），歷史單 marker 在新 marker 到貨時自動併入 |
| Step 6 唯一物種＝交集只有 1 個 | 另需：每個 marker 最高分（差 <= 0.2%）、與次佳非同義種差 >= 0.8%（兩值內部暫定，非 CLSI 查證值） |
| Step 6 下次建議 marker 只在「沒有任何候選」時給 | 只要沒有確定的唯一物種就給，排除該株已做過（含 NOT_BLASTED）的 marker；引子為 null 時標「引子待查，不可下單」 |
| Step 6 狀態字串寫「兩組」 | 規則是「2 組以上」，三個 marker 時字面仍沿用「兩組」（README 技術債 #7） |
| Step 7 內嵌 Python 腳本，以 `organism_guess` 分組 | 腳本在 `pipeline/cross_batch_summary.py`；`organism_guess` 現已是 taxid 正規化後的現行名稱，type 第二輪的 hit 另存 `blast_type_hits.csv`，不進這張備查表 |
| 注意事項：同一 sample+marker 恰好 2 條 reads | 同（盲審 🟡-9 尚未改：3 條 reads 或單條 read 仍整組 FAIL） |

### P5b 線上實測（2026-10-30，S2604280595 的 2 個樣本、3 個 marker）

樣本 61aISTW5c（Btub、CF）與 66aSCN231a（TEF1）；每個 marker 送全庫＋type 兩輪，hitlist 1000，信箱以臨時環境變數傳入（未寫入任何檔案）。

- 全庫輪：NCBI 回傳 1000 筆 alignment（Btub 其中 147 筆為合併了相同序列紀錄的 hit，最多併 21 筆），hitlist 1000 有效。
- type 輪（`sequence from type[filter]`）：可用，回傳 1000 筆；*61aISTW5c* Btub 的 type 序列 *P. copticola* 99.11% 達標（`TYPE_OK`，候選只剩 1 個種）；CF（98.75%）與 TEF1（*T. orientale* 95.56%）的 type 序列未達 99%，標 `NO_TYPE_SUPPORT`，沿用全庫候選並降信心。最終 61aISTW5c = 「交集唯一（待覆核）」，判讀來源 `Btub:type；CF:全庫(無 type 支持)`。
- 學名查詢：三個 marker 的 accession 全部對得到 taxid（0 個查不到）。
- **營運注意（NCBI 佇列）：** type 輪第一次送出的 4 個 request（Btub、CF、TEF1）等了 40–70 分鐘仍在 WAITING（全庫輪同一時段約 1.5–10 分鐘就回來）；放棄後重新送出同樣查詢，CF／TEF1／Btub 約 1 分鐘內就回來。結論：type 輪偶爾會卡在 NCBI 佇列，卡超過約 15 分鐘就中止、重送一次比空等有效。程式目前沒有自動逾時重送（Biopython 的 qblast 只會警告「超過 10 分鐘，建議重送」）。
- 另一次實測中途本機網路瞬斷（Connection refused）：type 輪失敗只降信心、整單不中斷（符合設計）；全庫輪失敗則該單中止，已存的 marker XML 在重跑時沿用。

## 技術債清單

1. **`Bio.pairwise2` 已被 Biopython 標示淘汰**（1.87 載入時即發 DeprecationWarning，官方說明未來版本會移除）。
   `sanger_pair_assemble.py` 的 `align_pair` 依賴它。`requirements.txt` 因此釘死 `biopython==1.87`；
   升級前必須先改用 `Bio.Align.PairwiseAligner`，並用既有舊批次做 consensus 一致性回歸（兩者在同分比對時的取捨可能不同）。
2. ~~Marker／引子辨識只涵蓋 ITS／Bt／EF／CF~~ **（Phase 3 已解）** 改由 `pipeline/markers.py` 從 `data/gene_lookup.json` 自動建映射；
   認不得的引子名歸「其他」並標 WARN、不丟資料。注意：同一樣本若有 2 條以上「其他」reads 會被當同一組 marker 互相配對，結果請人工確認。
3. **BLAST 走 NCBI qblast 線上 API**：慢、有限速、結果會隨 nt 資料庫更新而變；大量時考慮本機 BLAST+ 與固定版本資料庫。
   盲審延伸建議（尚未做）：存檔時記 `BlastOutput_db`／版本與判讀依據的 accession 清單；長期改本機 BLAST+ 配固定版本策展庫（P5b 的 type 第二輪與 taxonomy 快取是過渡方案）。
4. 腳本以 `python3 -m pipeline.xxx` 從專案根執行，尚未打包成可安裝套件（Phase 3 網頁介面前視需要處理）。
5. `run_all.py` 不記錄上次跑到哪（無 state 檔），從 summarize 起續跑需手動給 `--orders`。
6. 回歸驗證目前只比對輸出 CSV；Assembly 的 `alignment_checkpoints/*.txt` 等輔助檔未逐檔比對。
7. SKILL.md Step 4 提到 `run_blast.py`、Step 6 判讀狀態字串寫「兩組」但規則寫「2 組以上」（3 個 marker 時字面仍沿用「兩組」），文字待 SSOT 統一。

## 回歸驗證（2026-10-30，NAS 唯讀、輸出寫 `test_output/`）

批次 `S2604280595`（2 個子批次、15 個樣本），用新腳本從 rawdata 重跑，BLAST 離線重解析 NAS 既有 XML：

```bash
mkdir -p test_output/reg_run/S2604280595_ncbi_blast/raw_blast_xml
cp "$NAS/定序整理結果/S2604280595_ncbi_blast/raw_blast_xml/"*.xml test_output/reg_run/S2604280595_ncbi_blast/raw_blast_xml/
python3 run_all.py --rawdata-dir "$NAS/rawdata" --output-dir test_output/reg_run --orders S2604280595 --offline --date 2026-06-15
diff -rq "$NAS/定序整理結果/S2604280595_combined_assembly" test_output/reg_run/S2604280595_combined_assembly
```

結果：QC（兩個 `_processed`）與 Assembly（含 FASTA、alignment_checkpoints）與舊輸出 **逐檔位元組相同**；
BLAST 的 `blast_top100_hits.csv` 1500 列僅 `blast_date` 欄不同；`blast_species_summary` 13 個樣本與舊總表完全一致；
Step 5 摘要表對舊成品 `2026-06-15_new_orders_taxonomy_summary.csv`（4 訂單 37 列）逐位元組相同；
Step 6 對舊 `2026-06-15_new_orders_marker_intersection_summary.csv`（28 樣本）判讀狀態／唯一物種／交集全部一致，
差異只在備註欄新增「BLAST 只到屬」與下次建議 marker 欄。
舊 `all_orders` 版有「沒有可用 BLAST 結果」等 SKILL.md 沒定義的狀態（Assembly 失敗、未送 BLAST 的樣本），新腳本只處理 taxonomy 表內的樣本，這類樣本不會出現在 Step 6 表中。

### P5b 回歸（2026-10-30，三單離線：S2604280595／S2605070135／S2601100231）

以 P5a 驗收後的輸出為基準，從 BLAST 步驟重新解析既有 XML（`run_all.py --from-step blast --offline`，NAS 唯讀、輸出寫 `test_output/p5b_regression/`）。
學名查詢走真實 E-utilities（共約 3,200 個 accession、0 個查不到；BLAST 本身沒有重送）。

- **QC／Assembly**：逐檔位元組與 P5a 一致（5 個 `_processed`、3 個 `_combined_assembly` 全部 `diff -rq` 無差異）。
- **只關 type 第二輪（`type_round: false`，單看 taxid／hitlist／合併紀錄的影響）**：4,200 筆 hit 的 `pident`／`qcov`／排名完全不變；34 個樣本中 25 列判讀有文字級變動，**唯一物種的 7 個樣本結論不變**。變動歸類：
  - 431 筆 hit 的物種名改變：367 筆只是 `sp` -> `sp.` 寫法統一；42 筆舊名／同義名轉現行名（Hypocrea lixii -> Trichoderma harzianum／lixii、Hypocrea orientalis -> Trichoderma orientale、Neosartorya fischeri -> Aspergillus fischeri、Eurotium repens -> Aspergillus pseudoglaucus，對應 🟡-1）；9 筆 `Aspergillus genomic`／`clone` 不再當種名（🟢-1）；7 筆標題寫 sp. 但 taxid 有具體種；6 筆 `cf.` 規則恢復生效（🟢-1）。
  - 519 筆 hit 帶有 BLAST 合併的相同序列紀錄，其物種計入候選（🔴-8）：例 62aSCN32a／48xOxTW2b 的 Btub 候選多出 *Cercospora sojina*（Aspergillus 序列被誤標成外屬），`MULTIPLE_GENERA_IN_TOP20` 旗標因此出現；這正是 type 第二輪要濾掉的雜訊。
  - 3 個 marker 的旗標／建議動作變動：51aIPCN2a、66aSCN231a 的 TEF1 因 Hypocrea／Trichoderma 同義名合併，`MULTIPLE_GENERA_IN_TOP20` 消失（建議動作 `CHECK_TAXONOMY_SYNONYMS_OR_ADD_MARKER` -> `SPECIES_NEEDS_MANUAL_TOP_HIT_REVIEW`，🟡-1）；62aSCN32a 的 Btub 反向多出 `MULTIPLE_GENERA_IN_TOP20`（外屬合併紀錄，🔴-8）。
- **開 type 第二輪但離線沒有 type XML**：42 個 query 全部 `TYPE_NOT_RUN`；原本 7 個「兩組 marker 交集為唯一物種」全部降為「交集唯一（待覆核）」，原因欄註明「未做 type 第二輪」。這是預期行為（沒有 type 序列支持就不給綠燈）；這 7 個樣本要跑過 type 第二輪（線上）才可能回到綠燈。
