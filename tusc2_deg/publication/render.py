"""Render each Fig 7 panel independently (no collage).

    python -m tusc2_deg.publication.render --panel a
    python -m tusc2_deg.publication.render --all
"""
from __future__ import annotations
import argparse
from pathlib import Path
from tusc2_deg.publication import pub_style, gsea_bar, running_es, volcano

# Default output directory: <repo-root>/figures. render.py lives at
# <repo-root>/tusc2_deg/publication/render.py, so parents[2] is the repo root.
DEFAULT_OUT = Path(__file__).resolve().parents[2] / "figures"

# Figure 7 is four panels.
PANELS = {
    "a": ("fig7_a_gsea_nes_bar", lambda: gsea_bar.render()),
    "b": ("fig7_b_oxphos_runningES",
          lambda: running_es.render_pair("HALLMARK_OXIDATIVE_PHOSPHORYLATION")),
    "c": ("fig7_c_cd8_volcano", lambda: volcano.render("cd8")),
    "d": ("fig7_d_nk_volcano", lambda: volcano.render("nk")),
}

def render_panel(letter: str, out_dir=DEFAULT_OUT) -> str:
    pub_style.apply()
    name, fn = PANELS[letter]
    base = Path(out_dir) / name
    pub_style.save(fn(), base)
    return str(base)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", choices=list(PANELS))
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    a = ap.parse_args()
    letters = list(PANELS) if a.all else [a.panel]
    for L in letters:
        print("wrote", render_panel(L, a.out))

if __name__ == "__main__":
    main()
