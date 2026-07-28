import pandas as pd
from pathlib import Path
from tusc2_deg.cd8 import gsea


def _fake_prerank(rnk, gene_sets, **kw):
    class R:
        res2d = pd.DataFrame({
            "Term": ["MSigDB_Hallmark_2020__Oxidative Phosphorylation",
                     "MSigDB_Hallmark_2020__TNF-alpha Signaling via NF-kB"],
            "NES": [-1.94, 3.06],
            "FDR q-val": [0.001, 0.001],
        })
    return R()


def test_run_prerank_writes_csv_offline(tmp_path):
    deg = pd.DataFrame({"gene": [f"g{i}" for i in range(50)],
                        "stat": list(range(50))})
    out = tmp_path / "gsea_toy.csv"
    res = gsea.run_prerank(deg, "toy", out, prerank_fn=_fake_prerank)
    assert out.exists() and len(res) == 2
