"""Canonical anchor lists shared by the Fig 7 panel chain.

Two constants with a single definition site so the panels and the single-cell
cross-check read the same lists: CYTOTOXIC (the NK cytotoxic-core genes) and
SPINE_TERMS (the seven pre-specified Hallmark programmes on the GSEA NES bar).
gsea_bar imports SPINE_TERMS; verify_singlecell imports CYTOTOXIC.
"""
from __future__ import annotations

CYTOTOXIC = ["NKG7", "PRF1", "GNLY", "GZMA", "GZMB", "KLRD1"]

# The seven pre-specified Hallmark programmes on the GSEA NES bar (panel a) — the
# canonical spine list; gsea_bar imports it from here so there is one source.
SPINE_TERMS = [
    "TNF-alpha Signaling via NF-kB",
    "Interferon Gamma Response",
    "Interferon Alpha Response",
    "E2F Targets",
    "MYC Targets V1",
    "G2-M Checkpoint",
    "Oxidative Phosphorylation",
]
