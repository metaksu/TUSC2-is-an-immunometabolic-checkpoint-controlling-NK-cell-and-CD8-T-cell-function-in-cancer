"""Panel a — GSEA Hallmark NES bar, NK and CD8 as side-by-side facets.

x = signed NES (right = up in TUSC2-negative). One bar per Hallmark gene set,
ordered by NES, colored by direction. Selects the coordinated programs that
carry the narrative spine (inflammatory TUSC2-neg vs proliferative+OXPHOS
TUSC2-pos). All stats/n/test/source go in the caption, not the panel.
"""
from __future__ import annotations
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tusc2_deg.publication import data, pub_style
from tusc2_deg.publication.ledger import SPINE_TERMS  # canonical spine list (one source)
from tusc2_deg.nk.gsea import _norm_term

def _hallmark(g: pd.DataFrame) -> pd.DataFrame:
    h = g[g["Term"].str.contains("Hallmark", case=False)].copy()
    h["term"] = h["Term"].map(_norm_term)
    h["NES"] = pd.to_numeric(h["NES"], errors="coerce")
    h["FDR q-val"] = pd.to_numeric(h["FDR q-val"], errors="coerce")
    return h.dropna(subset=["NES"])

def select_hallmark_sets() -> dict[str, pd.DataFrame]:
    out = {}
    want = {t.lower() for t in SPINE_TERMS}
    for pipeline in ("nk", "cd8"):
        h = _hallmark(data.load_gsea(pipeline))
        sel = h[h["term"].str.lower().isin(want)].copy()
        out[pipeline] = sel.sort_values("NES")
    return out

def render(fig=None):
    sel = select_hallmark_sets()
    # union of terms present in either lineage, ordered by mean NES for a stable axis
    order = (pd.concat([sel["nk"], sel["cd8"]])
             .groupby("term")["NES"].mean().sort_values().index.tolist())
    if fig is None:
        fig, axes = plt.subplots(1, 2, figsize=(pub_style.COLUMN_FULL_IN, 2.2),
                                 sharey=True, layout="constrained")
    else:
        axes = fig.subplots(1, 2, sharey=True)
    y = np.arange(len(order))
    for ax, pipeline, title in zip(axes, ("nk", "cd8"), ("NK", "CD8")):
        s = sel[pipeline].set_index("term").reindex(order)
        nes = s["NES"].to_numpy()
        colors = [pub_style.UP_COLOR if (v is not None and v > 0) else pub_style.DOWN_COLOR
                  for v in nes]
        ax.barh(y, np.nan_to_num(nes), color=colors, edgecolor="0.2", linewidth=0.4)
        ax.axvline(0, color="0.4", lw=0.6)
        ax.set_xlabel("NES")
        ax.text(0.5, 1.02, title, transform=ax.transAxes, ha="center",
                va="bottom", fontsize=7, fontweight="bold")
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(order, fontsize=6)
    return fig
