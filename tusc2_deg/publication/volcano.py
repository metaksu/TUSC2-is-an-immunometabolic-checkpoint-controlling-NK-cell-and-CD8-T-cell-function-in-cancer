"""Panels c (CD8) and d (NK) — clean pseudobulk DESeq2 volcanoes.

Reuses each pipeline's significance gate (config.is_significant) and curated
gene-label set (the audit provenance table), rendering a clean 3-category
scatter (up / down / NS) with dashed threshold lines and repelled italic gene
labels. No title, no UP/DOWN counts, no FC-decode formula, no inline legend —
all of that lives in the caption.
"""
from __future__ import annotations
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text
from tusc2_deg.publication import data, pub_style
from tusc2_deg.nk import config as nk_config
from tusc2_deg.cd8 import config as cd8_config

_CFG = {"nk": nk_config, "cd8": cd8_config}


def curated_labels(pipeline: str) -> pd.DataFrame:
    cfg = _CFG[pipeline]
    for name in ("curated_label_provenance.tsv", "volcano_label_provenance.tsv"):
        p = cfg.AUDIT_DIR / name
        if p.exists():
            return pd.read_csv(p, sep="\t")
    raise FileNotFoundError(f"no curated-label provenance for {pipeline}")


# Explicit figure anchors added on top of the curated MUST_LABEL selection.
# HSPA6 is the lone CD8 gene past the 1.5-fold rung (log2FC/padj in output/degs)
# and is called out in the text, so it is anchored here; the curated MUST_LABEL
# set stays the source of truth.
ANCHORS = {"cd8": {"HSPA6"}, "nk": set()}


def label_genes(pipeline: str) -> set[str]:
    """The curated 'genes of interest' for this context
    (config.MUST_LABEL_BY_CONTEXT), plus any explicit figure anchors (HSPA6).
    render() filters to genes present in the table, so an absent gene is silently
    skipped."""
    cfg = _CFG[pipeline]
    must = set(cfg.MUST_LABEL_BY_CONTEXT.get(data.MAIN_CTX[pipeline], {}))
    return must | ANCHORS.get(pipeline, set())


# How many strongest unbiased non-curated hits to disclose in neutral grey
# (anti-cherry-picking device): global top-5 for NK, top-3 per direction for CD8.
_HONEST_TOPN = {"nk": ("global", 5), "cd8": ("per_direction", 3)}

# Display names + a stable order for the functional-family label-colour legend.
_FAMILY_DISPLAY = {
    "effector": "effector", "activation": "activation", "exhaustion": "exhaustion",
    "metabolism": "metabolism", "prolif": "proliferation", "progenitor": "progenitor",
    "stress": "stress/IEG", "other": "other",
}
_FAMILY_ORDER = list(_FAMILY_DISPLAY)


def _panel_legends(fig, pipeline, cfg, curated, extra, shown, family_colors):
    """Two legends beneath the panel (like the running-ES legend): point colour =
    significance/direction; label-text colour = functional family. Only families
    present among this panel's labels are shown."""
    from matplotlib.lines import Line2D
    def dot(c, lab):
        return Line2D([0], [0], marker="o", ls="", mfc=c, mec="none", ms=5, label=lab)
    def sq(c, lab):
        return Line2D([0], [0], marker="s", ls="", mfc=c, mec="none", ms=5, label=lab)
    handles = [dot(pub_style.UP_COLOR, "higher in TUSC2$^-$"),
               dot(pub_style.DOWN_COLOR, "higher in TUSC2$^+$"),
               dot(pub_style.NS_COLOR, "not significant")]
    if family_colors:
        for f in _FAMILY_ORDER:
            if any(g in shown and cfg.gene_family(g) == f for g in curated):
                handles.append(sq(pub_style.FAMILY_COLORS.get(f, "#222222"),
                                   _FAMILY_DISPLAY[f]))
    if extra & shown:
        handles.append(sq(pub_style.HONEST_GREY, "top unbiased"))
    # One legend, row-major with ncol=3: row 1 = the three point colours (circles,
    # significance/direction); rows below = gene-family label colours (squares).
    fig.legend(handles=handles, loc="outside lower center", ncol=3, fontsize=6,
               frameon=False, handletextpad=0.3, columnspacing=1.2, labelspacing=0.5)


def _composite(d: pd.DataFrame) -> pd.Series:
    """|log2FC| * -log10(padj) — the volcano label-ranking composite."""
    return d["log2FoldChange"].abs() * -np.log10(d["padj"].clip(lower=1e-300))


def honest_extra_genes(df: pd.DataFrame, pipeline: str, curated: set[str]) -> set[str]:
    """Strongest SIGNIFICANT non-curated genes by composite — disclosed in grey so
    the curated panel cannot read as a cherry-picked top-hit list."""
    mode, k = _HONEST_TOPN[pipeline]
    sig = df[df["_sig"] & ~df["gene"].isin(curated)].copy()
    if sig.empty:
        return set()
    sig["_cs"] = _composite(sig)
    if mode == "per_direction":
        picks = pd.concat([sig[sig["log2FoldChange"] > 0].nlargest(k, "_cs"),
                           sig[sig["log2FoldChange"] < 0].nlargest(k, "_cs")])
    else:
        picks = sig.nlargest(k, "_cs")
    return set(picks["gene"])


def render(pipeline: str, fig=None, figsize: tuple[float, float] = (5.6, 5.4),
           family_colors: bool = True):
    cfg = _CFG[pipeline]
    df = data.load_pseudobulk(pipeline).dropna(subset=["log2FoldChange", "padj"]).copy()
    df = df[df["gene"] != "TUSC2"]
    df["nlp"] = -np.log10(df["padj"].clip(lower=1e-300))
    up = cfg.is_significant(df, "UP")
    dn = cfg.is_significant(df, "DOWN")
    df["_sig"] = up | dn
    ns = ~df["_sig"]

    # Significance-driven symmetric x-cap: NS junk genes (e.g. NK RGS17, SERPINE1;
    # values in output/degs) must NOT stretch the axis and compress the real biology into
    # the centre. Cap = max(0.6, max|sig LFC| * 1.1); points beyond are clipped to
    # the edge and flagged with an off-scale triangle so the cap is disclosed.
    sig_abs = df.loc[df["_sig"], "log2FoldChange"].abs()
    xcap = max(0.6, float(sig_abs.max()) * 1.1) if len(sig_abs) else 0.6
    df["x_plot"] = df["log2FoldChange"].clip(-xcap, xcap)
    df["_xclip"] = df["log2FoldChange"].abs() > xcap
    strong_cut = cfg.SIG_LFC_LADDER[-1]                       # 0.585 = 1.5x
    strong = df["log2FoldChange"].abs() >= strong_cut

    if fig is None:
        fig, ax = plt.subplots(figsize=figsize, layout="constrained")
    else:
        ax = fig.subplots(1, 1)

    # NS cloud (grey), then significant points by direction; the |LFC|>=0.585
    # strong tail (CD8 only — NK has 0) is outlined + enlarged so the load-bearing
    # large-effect tier is distinguishable from the near-threshold band.
    ax.scatter(df.loc[ns, "x_plot"], df.loc[ns, "nlp"], s=3, c=pub_style.NS_COLOR,
               alpha=0.4, edgecolors="none", rasterized=True)
    for mask, color in ((up, pub_style.UP_COLOR), (dn, pub_style.DOWN_COLOR)):
        band, hi = mask & ~strong, mask & strong
        if band.any():
            ax.scatter(df.loc[band, "x_plot"], df.loc[band, "nlp"], s=5, c=color,
                       alpha=0.85, edgecolors="none", rasterized=True)
        if hi.any():
            ax.scatter(df.loc[hi, "x_plot"], df.loc[hi, "nlp"], s=9, c=color,
                       alpha=0.95, edgecolors="0.25", linewidths=0.3)

    ax.axhline(-np.log10(cfg.PADJ_THRESH), color="0.7", lw=0.5, ls="--")
    for v in (cfg.LFC_THRESH, -cfg.LFC_THRESH):
        ax.axvline(v, color="0.7", lw=0.5, ls="--")
    if (df["_sig"] & strong).any():                          # CD8: 1.5x guide
        for v in (strong_cut, -strong_cut):
            if abs(v) <= xcap:
                ax.axvline(v, color="0.82", lw=0.5, ls=":")

    # off-scale directional markers at the clipped edge (disclose the cap)
    for sub, mk in ((df["_xclip"] & (df["log2FoldChange"] > 0), ">"),
                    (df["_xclip"] & (df["log2FoldChange"] < 0), "<")):
        if sub.any():
            ax.scatter(df.loc[sub, "x_plot"], df.loc[sub, "nlp"], s=24, marker=mk,
                       facecolors="none", edgecolors="0.45", linewidths=0.6, zorder=5)
    ax.set_xlim(-xcap * 1.04, xcap * 1.04)

    # Labels: curated genes coloured by functional FAMILY (hybrid — dots stay
    # direction-coloured); plus the strongest unbiased non-curated hits in neutral
    # grey (honest disclosure). All italic, repelled with a fixed-seed layout
    # (no expand_axes, so the x-cap is not undone by label placement).
    curated = label_genes(pipeline)
    extra = honest_extra_genes(df, pipeline, curated)
    show = df[df["gene"].isin(curated | extra)]
    texts = []
    for r in show.itertuples():
        if r.gene in extra:
            col = pub_style.HONEST_GREY                       # grey honest extra
        elif family_colors:
            col = pub_style.FAMILY_COLORS.get(cfg.gene_family(r.gene), "#222222")
        else:
            col = "0.1"
        texts.append(ax.text(r.x_plot, r.nlp, r.gene, fontsize=6.5,
                             fontstyle="italic", color=col))
    if texts:
        np.random.seed(0)
        adjust_text(texts, ax=ax, expand=(2.4, 2.8), force_text=(1.0, 1.6),
                    force_static=(0.4, 0.7),
                    arrowprops=dict(arrowstyle="-", color="0.5", lw=0.4,
                                    shrinkA=2, shrinkB=2))

    ax.set_xlabel("log$_2$ fold change  (TUSC2$^-$ vs TUSC2$^+$)")
    ax.set_ylabel("-log$_{10}$ adjusted P")
    _panel_legends(fig, pipeline, cfg, curated, extra, set(show["gene"]), family_colors)
    return fig
