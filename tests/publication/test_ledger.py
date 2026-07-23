"""ledger.py contracts: the canonical anchor lists the panel chain imports."""
from tusc2_deg.publication import ledger


def test_spine_terms_are_seven_hallmark_programmes():
    # gsea_bar selects panel a's bars against this list (one source of truth).
    assert len(ledger.SPINE_TERMS) == 7
    for term in ("Oxidative Phosphorylation", "E2F Targets",
                 "Interferon Gamma Response"):
        assert term in ledger.SPINE_TERMS


def test_cytotoxic_core_genes():
    # verify_singlecell reads the NK cytotoxic core from here.
    assert ledger.CYTOTOXIC[:3] == ["NKG7", "PRF1", "GNLY"]
    assert "GZMB" in ledger.CYTOTOXIC
