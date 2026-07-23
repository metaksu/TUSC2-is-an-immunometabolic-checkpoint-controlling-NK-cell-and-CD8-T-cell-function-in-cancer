"""gsea.py contracts: signed-rank construction + offline prerank I/O.

The network prerank call is injected (prerank_fn) so unit tests run offline.
"""
import pandas as pd
import numpy as np

from tusc2_deg.cd8 import gsea


def test_build_rank_dedup_dropna_sorted():
    deg = pd.DataFrame({"gene": ["A", "B", "B", "C"],
                        "stat": [1.0, 3.0, 2.0, np.nan]})
    r = gsea.build_rank(deg)
    # C excluded (NaN); B deduped (keep first=3.0); descending by stat
    assert list(r.index) == ["B", "A"]
    assert r.loc["B"] == 3.0
    assert "C" not in r.index


def test_run_prerank_injected_fn_writes_and_passes_libraries(tmp_path):
    deg = pd.DataFrame({"gene": ["A", "B", "C"], "stat": [3.0, 1.0, -2.0]})
    captured = {}

    class FakeRes:
        res2d = pd.DataFrame({"Term": ["T1"], "NES": [1.5], "FDR q-val": [0.01]})

    def fake_prerank(rnk=None, gene_sets=None, **kw):
        captured["rnk"] = rnk
        captured["gene_sets"] = gene_sets
        return FakeRes()

    out = tmp_path / "gsea.csv"
    df = gsea.run_prerank(deg, "CD8_TIL", out, libraries=["LibA", "LibB"],
                          prerank_fn=fake_prerank)
    assert out.exists()
    assert "NES" in df.columns
    assert captured["gene_sets"] == ["LibA", "LibB"]
    # rnk passed to gseapy as a 2-col (gene, score) frame, descending
    rnk = captured["rnk"]
    assert list(rnk.iloc[:, 0]) == ["A", "B", "C"]
