"""引子名 -> marker 映射，以及正向／反向判斷。

映射從 data/gene_lookup.json 的引子名自動建立（不寫死清單），再加上前綴規則
處理 lookup 沒列、但實驗室常見的變體（如 ITS1／ITS4、EF-1、Bt1a）。
認不得的引子名歸「其他」（UNKNOWN_MARKER），由呼叫端標 WARN，不丟資料。
大小寫不敏感。
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
LOOKUP_PATH = PROJECT_ROOT / "data" / "gene_lookup.json"

UNKNOWN_MARKER = "其他"

# 前綴規則（小寫）。順序無關，取最長前綴；用於 lookup 沒列到的變體。
PREFIX_RULES: dict[str, str] = {
    "its": "ITS",
    "bt": "Btub",
    "ef": "TEF1",
    "cf": "CF",
    "cmd": "CF",
    "cal": "CF",
    "5f": "RPB2",
    "7c": "RPB2",
    "rpb2": "RPB2",
    "gpd": "GAPDH",
    "act": "ACT",
    "nl-": "LSU",
    "lr": "LSU",
    "tub": "Btub",
}

# lookup 只在備註文字提到、沒有完整引子節點的引子（T22 = Nigri 用 Btub 反向引子）
EXTRA_PRIMERS: dict[str, dict[str, str]] = {
    "t22": {"name": "T22", "marker": "Btub", "role": "rev"},
}

# 舊版就有、lookup 不一定列出的正向引子（維持既有排序行為）
LEGACY_FORWARD = {"ITS5", "ITS1", "BT2A", "EF-1", "CF1"}


def _walk_primers(node, marker=None):
    """遞迴找出所有 {name, seq} 引子節點，連同所屬 marker 與方向（fwd/rev）。"""
    if isinstance(node, dict):
        if isinstance(node.get("marker"), str):
            marker = node["marker"]
        primers = node.get("primers")
        if isinstance(primers, dict):
            for role in ("fwd", "rev"):
                p = primers.get(role)
                if isinstance(p, dict) and p.get("name"):
                    yield p["name"], marker, role
        for value in node.values():
            yield from _walk_primers(value, marker)
    elif isinstance(node, list):
        for value in node:
            yield from _walk_primers(value, marker)


@lru_cache(maxsize=4)
def load_primer_table(path: str | None = None) -> dict[str, dict[str, str]]:
    """回傳 {引子名小寫: {"name", "marker", "role"}}。同名多次出現取第一筆。"""
    data = json.loads(Path(path or LOOKUP_PATH).read_text(encoding="utf-8"))
    table: dict[str, dict[str, str]] = {}
    for name, marker, role in _walk_primers(data):
        table.setdefault(name.lower(), {"name": name, "marker": marker or "", "role": role})
    for key, value in EXTRA_PRIMERS.items():
        table.setdefault(key, dict(value))
    return table


def lookup_primer_names(path: str | None = None) -> list[str]:
    return sorted(v["name"] for v in load_primer_table(path).values())


def marker_for_primer(primer: str) -> str:
    """引子名 -> marker；認不得回傳「其他」。"""
    key = primer.strip().lower()
    if not key:
        return UNKNOWN_MARKER
    hit = load_primer_table().get(key)
    if hit and hit["marker"]:
        return hit["marker"]
    for prefix in sorted(PREFIX_RULES, key=len, reverse=True):
        if key.startswith(prefix):
            return PREFIX_RULES[prefix]
    return UNKNOWN_MARKER


def is_forward(primer: str) -> bool:
    """正向引子判斷：lookup 的 fwd 角色、舊版清單；其餘視為反向。"""
    if primer.upper() in LEGACY_FORWARD:
        return True
    hit = load_primer_table().get(primer.strip().lower())
    return bool(hit and hit["role"] == "fwd")
