"""Panel b — GSEA running-enrichment-score 'mountain' curve.

Recomputes the classic Subramanian running ES from the cached pseudobulk `stat`
ranking against a vendored Hallmark .grp gene set (the cached gsea CSV stores
only the summary, not the per-rank curve). Deterministic, offline, no raw cells.
NES and FDR (read from the cached gsea table) are annotated inside the panel —
the one GSEA idiom where stats legitimately sit on-panel.
"""
from __future__ import annotations
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tusc2_deg.publication import data, pub_style
from tusc2_deg.nk.gsea import build_rank, _norm_term
from tusc2_deg.nk import config as nk_config

# .grp filename per Hallmark set id (vendored under gene_sets/).
GRP = {
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION":
        "HALLMARK_OXIDATIVE_PHOSPHORYLATION.v2025.1.Hs.grp",
    "HALLMARK_E2F_TARGETS": "HALLMARK_E2F_TARGETS.MSigDB_Hallmark_2020.grp",
    "HALLMARK_G2M_CHECKPOINT": "HALLMARK_G2M_CHECKPOINT.MSigDB_Hallmark_2020.grp",
}

def _members(set_id: str) -> set[str]:
    return nk_config._load_grp(GRP[set_id])

def compute(pipeline: str, set_id: str, weight: float = 1.0) -> dict:
    ranked = build_rank(data.load_pseudobulk(pipeline))  # gene -> stat, descending
    genes = ranked.index.to_numpy()
    members = _members(set_id) & set(genes)
    hit = np.array([g in members for g in genes])
    r = np.abs(ranked.to_numpy()) ** weight
    n = len(genes)
    nhit = hit.sum()
    p_hit = np.cumsum(np.where(hit, r, 0.0)) / max((r[hit].sum()), 1e-12)
    p_miss = np.cumsum(np.where(~hit, 1.0, 0.0)) / max((n - nhit), 1)
    es_curve = p_hit - p_miss
    # cached summary numbers for annotation
    g = data.load_gsea(pipeline)
    target = _norm_term(set_id.replace("HALLMARK_", "").replace("_", " ")).lower()
    m = g[g["Term"].map(_norm_term).str.lower() == target]
    # Prefer the Hallmark row (the NES the panel a bar uses) over the GO/Reactome
    # rows of the same pathway name (the NK OXPHOS Hallmark and GO NES differ
    # slightly; live values in output/gsea) so the annotation matches panel a.
    hall = m[m["Term"].str.contains("Hallmark", case=False)]
    row = hall if not hall.empty else m
    nes = float(pd.to_numeric(row["NES"], errors="coerce").iloc[0]) if len(row) else float("nan")
    fdr = float(pd.to_numeric(row["FDR q-val"], errors="coerce").iloc[0]) if len(row) else float("nan")
    stat = ranked.to_numpy()
    zero_cross = int(np.argmax(stat <= 0)) if (stat <= 0).any() else len(stat)
    return {"es_curve": es_curve, "hit_index": np.where(hit)[0],
            "ranked_stat": stat, "nes": nes, "fdr": fdr, "zero_cross": zero_cross}

def _draw(fig, res, *, es_ylim=None, title=None):
    """Draw one GSEA enrichment plot (ES profile + hit barcode + correlation
    gradient + ranked-metric waterfall + zero-cross) into `fig` (Figure/SubFigure)."""
    import matplotlib.colors as mcolors
    gs = fig.add_gridspec(4, 1, height_ratios=[3.0, 0.32, 0.26, 1.25], hspace=0.06)
    ax_es = fig.add_subplot(gs[0])
    ax_hit = fig.add_subplot(gs[1])
    ax_grad = fig.add_subplot(gs[2])
    ax_rank = fig.add_subplot(gs[3])
    x = np.arange(len(res["es_curve"]))

    ax_es.plot(x, res["es_curve"], color="#117733", lw=1.3)
    ax_es.axhline(0, color="0.6", lw=0.5)
    ax_es.grid(True, ls=":", lw=0.3, color="0.85")
    ax_es.set_axisbelow(True)
    ax_es.set_ylabel("Enrichment score (ES)")
    ax_es.set_xticks([])
    ax_es.set_xlim(0, len(x))
    if es_ylim is not None:
        ax_es.set_ylim(*es_ylim)
    if title:
        ax_es.set_title(title, fontsize=8, fontweight="bold")
    nes, fdr = res["nes"], res["fdr"]
    fdr_txt = "= 0.001" if fdr <= 0.001 else f"= {fdr:.3f}"
    ax_es.text(0.5, 0.93, f"NES = {nes:.2f}    FDR q {fdr_txt}",
               transform=ax_es.transAxes, ha="center", va="top", fontsize=6.5)

    ax_hit.vlines(res["hit_index"], 0, 1, color="0.15", lw=0.4)
    ax_hit.set_xlim(0, len(x))
    ax_hit.set_axis_off()

    cmap = mcolors.LinearSegmentedColormap.from_list(
        "tusc2div", [pub_style.UP_COLOR, "#f7f7f7", pub_style.DOWN_COLOR])
    ax_grad.imshow(np.linspace(0, 1, len(x)).reshape(1, -1), aspect="auto",
                   cmap=cmap, extent=(0, len(x), 0, 1))
    ax_grad.set_xlim(0, len(x))
    ax_grad.set_xticks([])
    ax_grad.set_yticks([])
    for s in ax_grad.spines.values():
        s.set_visible(False)

    ax_rank.fill_between(x, res["ranked_stat"], color="0.7", lw=0)
    ax_rank.axhline(0, color="0.6", lw=0.5)
    zc = res["zero_cross"]
    ax_rank.axvline(zc, color="0.5", lw=0.6, ls="--")
    ax_rank.text(zc + len(x) * 0.01, ax_rank.get_ylim()[1] * 0.82,
                 f"zero cross, rank {zc:,}", fontsize=5.5, va="top",
                 ha="left", color="0.4")
    ax_rank.set_xlim(0, len(x))
    ax_rank.set_ylabel("Ranked metric")
    ax_rank.set_xlabel("Gene rank  (TUSC2$^-$ → TUSC2$^+$)")
    return ax_es


def _legend(fig):
    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], color="#117733", lw=1.3, label="Enrichment profile"),
               Line2D([0], [0], color="0.15", lw=0.8, label="Gene-set hits"),
               Line2D([0], [0], color="0.7", lw=4, label="Ranked metric")]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, fontsize=5.5,
               frameon=False, handlelength=1.4, columnspacing=1.4)


def render(pipeline: str, set_id: str, fig=None,
           figsize: tuple[float, float] = (4.6, 3.8)):
    """Single-lineage GSEA enrichment plot."""
    res = compute(pipeline, set_id)
    if fig is None:
        fig = plt.figure(figsize=figsize, layout="constrained")
    _draw(fig, res)
    _legend(fig)
    return fig


def render_pair(set_id, pipelines=("nk", "cd8"), labels=("NK", "CD8"), fig=None,
                figsize: tuple[float, float] = (7.6, 3.9)):
    """Paired enrichment plots (NK and CD8 OXPHOS) sharing the ES y-axis so the
    depth difference — NK shallow/mid-pack vs CD8 deep/discrete — reads honestly
    (the TUSC2 selective-vs-broad metabolic asymmetry). NK is drawn left, CD8
    right, matching panel a. The lineage label is the only on-panel identifier;
    the pathway name goes in the caption."""
    if fig is None:
        fig = plt.figure(figsize=figsize, layout="constrained")
    res = {p: compute(p, set_id) for p in pipelines}
    es_min = min(r["es_curve"].min() for r in res.values())
    ylim = (es_min * 1.08, 0.06)
    subs = fig.subfigures(1, len(pipelines), wspace=0.06)
    for sf, p, lab in zip(subs, pipelines, labels):
        _draw(sf, res[p], es_ylim=ylim, title=lab)
    _legend(fig)
    return fig
