"""讀取 config.yaml，並解析出各路徑。流程腳本一律從這裡拿設定。"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

import yaml

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CONFIG_PATH = PROJECT_ROOT / "config.yaml"


def load_config(path: str | os.PathLike | None = None) -> dict[str, Any]:
    """載入設定。順序：參數 path > 環境變數 SANGER_CONFIG > 專案根的 config.yaml。

    會把 cfg["_base"] 設為設定檔所在資料夾（相對路徑的基準）。
    環境變數 SANGER_NAS_ROOT / SANGER_NCBI_EMAIL 可覆蓋對應欄位。
    """
    cfg_path = Path(path or os.environ.get("SANGER_CONFIG") or DEFAULT_CONFIG_PATH)
    with cfg_path.open(encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle)
    cfg["_base"] = cfg_path.resolve().parent
    if os.environ.get("SANGER_NAS_ROOT"):
        cfg["paths"]["nas_root"] = os.environ["SANGER_NAS_ROOT"]
    if os.environ.get("SANGER_NCBI_EMAIL"):
        cfg["blast"]["email"] = os.environ["SANGER_NCBI_EMAIL"]
    if os.environ.get("SANGER_ENABLE_LOCAL_MODE"):
        cfg.setdefault("app", {})["enable_local_mode"] = os.environ["SANGER_ENABLE_LOCAL_MODE"] not in ("0", "false", "")
    return cfg


def _resolve(cfg: dict[str, Any], value: str) -> Path:
    p = Path(value).expanduser()
    return p if p.is_absolute() else (cfg["_base"] / p)


def rawdata_dir(cfg: dict[str, Any]) -> Path:
    paths = cfg["paths"]
    if paths.get("rawdata_dir"):
        return _resolve(cfg, paths["rawdata_dir"])
    return _resolve(cfg, paths["nas_root"]) / paths["rawdata_subdir"]


def output_dir(cfg: dict[str, Any]) -> Path:
    paths = cfg["paths"]
    if paths.get("output_dir"):
        return _resolve(cfg, paths["output_dir"])
    return _resolve(cfg, paths["nas_root"]) / paths["output_subdir"]


def gene_lookup_path(cfg: dict[str, Any]) -> Path:
    return _resolve(cfg, cfg["paths"]["gene_lookup"])


def blast_settings(cfg: dict[str, Any]) -> dict[str, Any]:
    """blast 區塊＋解析好的 taxonomy 快取檔絕對路徑（taxid 對照表落檔重用）。"""
    out = dict(cfg["blast"])
    cache = out.get("taxonomy_cache") or ""
    out["taxonomy_cache_file"] = str(_resolve(cfg, cache)) if cache else ""
    return out
