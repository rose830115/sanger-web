"""accession -> taxid -> NCBI Taxonomy 現行名稱（修正 🔴-5、🟡-1、🟢-1）。

為什麼：學名原本從序列標題前兩個字抓，GenBank 標題會有舊名／同義名／亂寫（Hypocrea lixii、Emericella、
Eurotium、section、complex、Synthetic construct…），字串比對交集會把同一種算成兩種。
做法：每筆 hit 的 accession 先查 taxid（E-utilities esummary，nuccore），再用 taxid 查 NCBI Taxonomy
（E-utilities efetch，taxonomy）得到現行名稱、階層與上層種／屬；舊名、同義名因此在交集前就被正規化。

限速與快取：
- 未帶 API key 每秒最多 3 次請求，這裡保守取 0.4 秒間隔；有環境變數 NCBI_API_KEY 才放寬到 0.12 秒。
- accession->taxid 與 taxid->名稱落檔（JSON）重用；taxid 名稱超過 max_age_days（內部暫定 90 天）才重查，
  因為分類名稱會隨時間修訂。快取檔與信箱無關，信箱只在 HTTP 請求參數裡，不落檔。
- 離線模式（offline=True）絕不連線，只吃快取；查不到的 accession 回傳 None，由呼叫端退回標題猜測。

限制（誠實揭露）：NCBI Taxonomy 的「現行名稱」不等於 MycoBank／Index Fungorum 的正式名稱，
入庫前仍建議對 MycoBank／Index Fungorum（README 已列）。
"""

from __future__ import annotations

import json
import os
import re
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import date
from pathlib import Path
from typing import Callable, Iterable

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
CACHE_VERSION = 1
SUMMARY_BATCH = 200
TAXONOMY_BATCH = 100
_STRICT_BINOMIAL = re.compile(r"^[A-Z][a-z]+(?:-[a-z]+)? [a-z][a-z-]+$")

HttpFn = Callable[[str, dict], bytes]


def default_http(url: str, params: dict) -> bytes:
    """POST 到 E-utilities；暫時性錯誤（429／5xx／連線）重試 3 次。"""
    data = urllib.parse.urlencode(params).encode()
    last: Exception | None = None
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, data=data, headers={"User-Agent": params.get("tool", "sanger-web")})
            with urllib.request.urlopen(req, timeout=60) as resp:  # noqa: S310
                return resp.read()
        except urllib.error.HTTPError as exc:
            last = exc
            if exc.code not in (429, 500, 502, 503, 504):
                break
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            last = exc
        time.sleep(2.0 * (attempt + 1))
    raise RuntimeError(f"E-utilities 請求失敗：{last}")


def classify_taxon(name: str, rank: str, lineage: list[dict]) -> dict:
    """由 taxonomy 節點算出 {species, species_taxid, genus, rank_kind}。

    species：只有「種（含種以下：亞種、變種、品系）」才有值，且必須是乾淨的二名法；
    section／complex／species group／genus／hybrid／環境樣本等一律 species=""（🟢-1）。
    """
    by_rank = {n["rank"]: n for n in lineage if n.get("rank") and n["rank"] != "no rank"}
    genus = by_rank["genus"]["name"] if "genus" in by_rank else (name if rank == "genus" else "")
    species_name = species_taxid = ""
    if rank == "species":
        species_name, species_taxid = name, ""
    elif "species" in by_rank and rank not in ("species group", "species subgroup", "section", "subsection",
                                              "series", "subgenus", "genus", "family", "order", "class"):
        species_name, species_taxid = by_rank["species"]["name"], by_rank["species"]["taxid"]
    lineage_names = {n["name"] for n in lineage}
    non_target = bool(lineage_names & {"environmental samples", "unclassified sequences", "artificial sequences",
                                       "other sequences", "synthetic construct"})
    if species_name:
        if non_target or not _STRICT_BINOMIAL.match(species_name):
            # sp.／hybrid（x）／變體寫法：不當二名法種名，退成「屬 sp.」讓交集排除（與舊行為一致）
            species_name = f"{genus} sp." if genus and not non_target else ""
    return {"species": species_name, "species_taxid": species_taxid, "genus": genus}


def _lineage_from_taxon(taxon: ET.Element) -> list[dict]:
    out = []
    for node in taxon.findall("./LineageEx/Taxon"):
        out.append({
            "taxid": (node.findtext("TaxId") or "").strip(),
            "name": (node.findtext("ScientificName") or "").strip(),
            "rank": (node.findtext("Rank") or "").strip(),
        })
    return out


def parse_taxonomy_xml(payload: bytes) -> dict[str, dict]:
    """efetch db=taxonomy XML -> {taxid: 節點資訊}；已合併的舊 taxid（AkaTaxIds）也各自指向新節點。"""
    root = ET.fromstring(payload)
    out: dict[str, dict] = {}
    for taxon in root.findall("./Taxon"):
        taxid = (taxon.findtext("TaxId") or "").strip()
        name = (taxon.findtext("ScientificName") or "").strip()
        rank = (taxon.findtext("Rank") or "").strip()
        if not taxid or not name:
            continue
        lineage = _lineage_from_taxon(taxon)
        info = {"taxid": taxid, "name": name, "rank": rank, **classify_taxon(name, rank, lineage)}
        if rank == "species":
            info["species_taxid"] = taxid
        out[taxid] = info
        for aka in taxon.findall("./AkaTaxIds/TaxId"):
            old = (aka.text or "").strip()
            if old and old not in out:
                out[old] = dict(info, merged_from=old)
    return out


class TaxonomyResolver:
    def __init__(
        self,
        cache_path: Path | None = None,
        email: str | None = None,
        tool: str = "sanger-web-blast",
        api_key: str | None = None,
        offline: bool = False,
        max_age_days: int = 90,
        http: HttpFn | None = None,
        sleep: Callable[[float], None] = time.sleep,
        today: Callable[[], date] = date.today,
    ) -> None:
        self.cache_path = Path(cache_path) if cache_path else None
        self.email = email
        self.tool = tool
        self.api_key = api_key if api_key is not None else os.environ.get("NCBI_API_KEY")
        self.offline = offline
        self.max_age_days = max_age_days
        self._http = http or default_http
        self._sleep = sleep
        self._today = today
        self.interval = 0.12 if self.api_key else 0.4
        self._last_call = 0.0
        self.cache = self._load()
        self.dirty = False
        self.errors: list[str] = []
        self.stats = {"acc_cached": 0, "acc_fetched": 0, "tax_cached": 0, "tax_fetched": 0, "unresolved": 0}
        self._missing_acc: set[str] = set()

    # ---- 快取 ----
    def _load(self) -> dict:
        if self.cache_path and self.cache_path.is_file():
            try:
                data = json.loads(self.cache_path.read_text(encoding="utf-8"))
                if data.get("version") == CACHE_VERSION:
                    return data
            except (OSError, ValueError):
                pass
        return {"version": CACHE_VERSION, "acc": {}, "tax": {}}

    def save(self) -> None:
        if not self.cache_path or not self.dirty:
            return
        from pipeline.stage_guard import atomic_write_text  # noqa: PLC0415

        self.cache_path.parent.mkdir(parents=True, exist_ok=True)
        atomic_write_text(self.cache_path, json.dumps(self.cache, ensure_ascii=False, sort_keys=True))
        self.dirty = False

    def _fresh(self, entry: dict) -> bool:
        try:
            age = (self._today() - date.fromisoformat(entry.get("fetched", "1970-01-01"))).days
        except ValueError:
            return False
        return age <= self.max_age_days

    # ---- HTTP ----
    def _call(self, endpoint: str, params: dict) -> bytes:
        wait = self.interval - (time.monotonic() - self._last_call)
        if wait > 0:
            self._sleep(wait)
        base = {"tool": self.tool}
        if self.email:
            base["email"] = self.email
        if self.api_key:
            base["api_key"] = self.api_key
        try:
            return self._http(EUTILS + endpoint, {**base, **params})
        finally:
            self._last_call = time.monotonic()

    # ---- 查詢 ----
    def _fetch_taxids(self, accessions: list[str]) -> None:
        for i in range(0, len(accessions), SUMMARY_BATCH):
            chunk = accessions[i : i + SUMMARY_BATCH]
            try:
                raw = self._call("esummary.fcgi", {"db": "nuccore", "id": ",".join(chunk), "retmode": "json"})
                result = json.loads(raw).get("result", {})
            except (RuntimeError, ValueError) as exc:
                self.errors.append(f"esummary：{exc}")
                continue
            got: set[str] = set()
            for uid in result.get("uids", []):
                rec = result.get(str(uid), {})
                taxid = str(rec.get("taxid") or "")
                if not taxid:
                    continue
                for key in (rec.get("accessionversion"), rec.get("caption")):
                    if key:
                        self.cache["acc"][key] = taxid
                        got.add(key)
            self.stats["acc_fetched"] += len(got & set(chunk))
            self._missing_acc |= {a for a in chunk if a not in self.cache["acc"]}
            self.dirty = True

    def _fetch_taxonomy(self, taxids: list[str]) -> None:
        for i in range(0, len(taxids), TAXONOMY_BATCH):
            chunk = taxids[i : i + TAXONOMY_BATCH]
            try:
                raw = self._call("efetch.fcgi", {"db": "taxonomy", "id": ",".join(chunk), "retmode": "xml"})
                nodes = parse_taxonomy_xml(raw)
            except (RuntimeError, ET.ParseError) as exc:
                self.errors.append(f"efetch taxonomy：{exc}")
                continue
            stamp = self._today().isoformat()
            for taxid, info in nodes.items():
                self.cache["tax"][taxid] = {**info, "fetched": stamp}
                self.stats["tax_fetched"] += 1
            self.dirty = True

    def resolve(self, accessions: Iterable[str]) -> dict[str, dict | None]:
        """{accession: 資訊或 None}。資訊：taxid、name（現行名稱）、rank、species、species_taxid、genus。"""
        wanted = list(dict.fromkeys(a for a in accessions if a))
        acc_map = self.cache["acc"]
        need_acc = [a for a in wanted if a not in acc_map and a not in self._missing_acc]
        self.stats["acc_cached"] += len(wanted) - len(need_acc)
        if need_acc and not self.offline:
            self._fetch_taxids(need_acc)
        taxids = list(dict.fromkeys(acc_map[a] for a in wanted if a in acc_map))
        tax = self.cache["tax"]
        need_tax = [t for t in taxids if t not in tax or not self._fresh(tax[t])]
        self.stats["tax_cached"] += len(taxids) - len(need_tax)
        if need_tax and not self.offline:
            self._fetch_taxonomy(need_tax)
        out: dict[str, dict | None] = {}
        for acc in wanted:
            taxid = acc_map.get(acc)
            info = tax.get(taxid) if taxid else None
            if info:
                out[acc] = info
            else:
                out[acc] = None
                self.stats["unresolved"] += 1
        self.save()
        return out
