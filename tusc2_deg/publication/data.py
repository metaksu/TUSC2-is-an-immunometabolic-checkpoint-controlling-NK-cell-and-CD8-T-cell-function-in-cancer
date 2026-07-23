"""Loaders for the cached pipeline tables (NK + CD8 main contexts).

One place that knows each pipeline's main context and output paths, so panels
and the ledger read identical data. Reuses each pipeline config's path
constants — no hard-coded output paths.
"""
from __future__ import annotations
import pandas as pd
from tusc2_deg.nk import config as nk_config
from tusc2_deg.cd8 import config as cd8_config

# main context id per pipeline (single-context main figure)
MAIN_CTX = {"nk": "Tang_Tumor", "cd8": "CD8_TIL"}
_CFG = {"nk": nk_config, "cd8": cd8_config}

def _cfg(pipeline: str):
    if pipeline not in _CFG:
        raise ValueError(f"pipeline must be 'nk' or 'cd8', got {pipeline!r}")
    return _CFG[pipeline]

def load_pseudobulk(pipeline: str) -> pd.DataFrame:
    cfg = _cfg(pipeline)
    return pd.read_csv(cfg.DEG_DIR / f"deg_pseudobulk_{MAIN_CTX[pipeline]}.csv")

def load_wilcoxon(pipeline: str) -> pd.DataFrame:
    cfg = _cfg(pipeline)
    return pd.read_csv(cfg.DEG_DIR / f"deg_wilcoxon_{MAIN_CTX[pipeline]}.csv")

def load_gsea(pipeline: str) -> pd.DataFrame:
    cfg = _cfg(pipeline)
    return pd.read_csv(cfg.GSEA_DIR / f"gsea_{MAIN_CTX[pipeline]}.csv")

def load_sig_ladder(pipeline: str) -> pd.DataFrame:
    cfg = _cfg(pipeline)
    return pd.read_csv(cfg.AUDIT_DIR / "deg_sig_ladder.tsv", sep="\t")
