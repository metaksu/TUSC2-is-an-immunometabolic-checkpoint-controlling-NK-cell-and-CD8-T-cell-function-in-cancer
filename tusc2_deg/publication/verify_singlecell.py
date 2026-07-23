"""Single-cell cross-check of the NK cytotoxic-core direction (RAW Tang atlas).

Recomputes the per-cell direction of the NK cytotoxic core (TUSC2-positive vs
TUSC2-negative) and confirms it matches the cached Wilcoxon table. TUSC2-status
is detection-based: a cell is TUSC2-positive iff it has >0 TUSC2 counts.
"""
from __future__ import annotations
import numpy as np
import scanpy as sc
from tusc2_deg.publication import data
from tusc2_deg.publication.ledger import CYTOTOXIC
from tusc2_deg.paths import TANG_H5AD as ATLAS

def recompute_cytotoxic_sc_direction() -> dict[str, bool]:
    a = sc.read_h5ad(ATLAS)
    if "TUSC2" not in a.var_names:
        raise RuntimeError("TUSC2 not in atlas var_names")
    tusc2 = np.asarray(a[:, "TUSC2"].X.todense()).ravel() if hasattr(a[:, "TUSC2"].X, "todense") \
        else np.asarray(a[:, "TUSC2"].X).ravel()
    pos = tusc2 > 0
    cached = data.load_wilcoxon("nk").set_index("gene")
    out = {}
    for g in CYTOTOXIC:
        if g not in a.var_names or g not in cached.index:
            continue
        x = np.asarray(a[:, g].X.todense()).ravel() if hasattr(a[:, g].X, "todense") \
            else np.asarray(a[:, g].X).ravel()
        mean_pos, mean_neg = x[pos].mean(), x[~pos].mean()
        sc_dir_raw = np.sign(mean_neg - mean_pos)         # neg-vs-pos
        sc_dir_cached = np.sign(cached.loc[g, "logfoldchanges"])
        # both should be negative (cytotoxic higher in TUSC2-positive)
        out[g] = bool(sc_dir_raw == sc_dir_cached and sc_dir_raw < 0)
    return out
