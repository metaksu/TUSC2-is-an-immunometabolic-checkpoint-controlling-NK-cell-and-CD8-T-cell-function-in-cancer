"""Gold-standard publication style for the Fig 7 panels.

Modeled on Nature/Cell specs: sans-serif Helvetica/Arial, 5-7 pt text, >=0.5 pt
strokes, colorblind-safe diverging palette (blue<->orange), vector PDF with
embedded fonts (not outlined) + 600-dpi PNG. Panels carry no letter/title.
"""
from __future__ import annotations
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt

# Colorblind-safe diverging pair (Okabe-Ito): orange = up in TUSC2-negative,
# blue = up in TUSC2-positive; neutral grey = NS.
UP_COLOR = "#E69F00"     # orange  (higher in TUSC2-negative / KO analogue)
DOWN_COLOR = "#0072B2"   # blue    (higher in TUSC2-positive / WT analogue)
NS_COLOR = "#BDBDBD"     # neutral grey

# CVD-safe (Paul Tol muted + Okabe-Ito) palette for the gene-FAMILY label-text
# colouring on the volcanoes. Hybrid scheme: dots are DIRECTION-coloured
# (UP/DOWN orange/blue), curated label TEXT is FAMILY-coloured with this set.
# Chosen to stay distinct from the orange/blue direction dots and to remain
# separable under deuteranopia/protanopia. The strongest-unbiased "honest
# disclosure" labels use HONEST_GREY, kept clearly lighter/neutral so they read
# as a distinct (uncurated) class.
FAMILY_COLORS = {
    "effector":   "#CC3311",   # red
    "metabolism": "#999933",   # olive
    "prolif":     "#44AA99",   # teal
    "stress":     "#882255",   # wine
    "exhaustion": "#AA4499",   # purple
    "activation": "#117733",   # green
    "progenitor": "#332288",   # indigo
    "other":      "#CC6677",   # rose
    "tusc2":      "#000000",
}
HONEST_GREY = "#8a8a8a"

# Nature double-column = 183 mm; convert to inches for figsize budgeting.
MM = 1.0 / 25.4
COLUMN_FULL_IN = 183 * MM
COLUMN_SINGLE_IN = 89 * MM

def apply() -> None:
    mpl.rcParams.update({
        "pdf.fonttype": 42,            # embed TrueType (editable text)
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        # Arial is the gold-standard sans (Nature: "preferably Helvetica or
        # Arial"). Helvetica is omitted because it is absent on Windows/Linux and
        # only emits noisy per-glyph "font not found" warnings before falling
        # through to Arial anyway; DejaVu Sans is the metric-compatible fallback.
        "font.family": ["Arial", "DejaVu Sans"],
        "font.size": 7,
        "axes.titlesize": 7,
        "axes.labelsize": 7,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "legend.fontsize": 6,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "lines.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "savefig.dpi": 600,
        "figure.dpi": 600,
        "savefig.bbox": "tight",
        "savefig.facecolor": "white",
    })

def save(fig, base_path) -> None:
    """Write both vector PDF and 600-dpi PNG; no on-figure title is added."""
    base = Path(base_path)
    base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(base.with_suffix(".pdf"))
    fig.savefig(base.with_suffix(".png"), dpi=600)
    plt.close(fig)
