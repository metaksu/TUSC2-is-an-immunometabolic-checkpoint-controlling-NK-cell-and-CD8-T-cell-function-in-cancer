"""Offline reproduction of the Figure 7 deliverables.

Regenerates, entirely from the cached pipeline tables under
tusc2_deg/{nk,cd8}/output and with no external data or network access:

  1. the four Figure 7 panels  -> figures/
       a  GSEA NES bar (NK + CD8)
       b  OXPHOS running-enrichment (NK + CD8)
       c  CD8 volcano
       d  NK volcano
  2. the collected Supplementary data files S1-S4 + Table S1 -> supplementary/
  3. the merged workbook Figure7_Supplementary_Data.xlsx     -> supplementary/

Run:  python reproduce.py
Force the Agg backend for headless rendering:  MPLBACKEND=Agg python reproduce.py
"""
from __future__ import annotations
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))               # import this folder's tusc2_deg package

import collect_supplementary
import make_supplementary_xlsx
from tusc2_deg.publication import render

FIGURES = HERE / "figures"


def render_panels() -> None:
    FIGURES.mkdir(exist_ok=True)
    print("[1/3] rendering panels ->", FIGURES)
    for letter in render.PANELS:
        print("  wrote", render.render_panel(letter, FIGURES))


def collect() -> int:
    print("[2/3] collecting Supplementary data files -> supplementary/")
    return collect_supplementary.main()


def build_xlsx() -> int:
    print("[3/3] building Figure7_Supplementary_Data.xlsx")
    return make_supplementary_xlsx.main()


def main() -> int:
    render_panels()
    rc = collect()
    rc |= build_xlsx()
    print("\nreproduce: done" if rc == 0 else "\nreproduce: FAILED (missing inputs)")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
