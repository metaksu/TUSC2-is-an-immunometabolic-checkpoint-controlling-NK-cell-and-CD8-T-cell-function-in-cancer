import math
import pytest
from tusc2_deg.publication import verify_singlecell as vsc

@pytest.mark.slow
def test_singlecell_wilcoxon_matches_cached(has_raw_tang):
    if not has_raw_tang:
        pytest.skip("raw Tang atlas not present")
    # Recompute single-cell logFC for the cytotoxic core from the raw atlas and
    # confirm the cached Wilcoxon table reproduces the same DIRECTION (strongly
    # higher in TUSC2-positive => negative neg-vs-pos logFC) for each gene.
    res = vsc.recompute_cytotoxic_sc_direction()
    for gene, agree in res.items():
        assert agree, f"{gene} single-cell direction disagrees with cached table"
