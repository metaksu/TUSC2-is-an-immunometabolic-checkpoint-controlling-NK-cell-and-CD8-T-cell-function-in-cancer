import pandas as pd
from tusc2_deg.cd8 import config


def test_is_significant_three_part_gate():
    # 3-part gate: padj-sig & baseMean-sig & |LFC|>LFC_THRESH.
    # Row 0: padj/baseMean OK but |LFC|=0.05 < LFC_THRESH -> NOT sig.
    # Row 1: |LFC| OK but baseMean too low -> NOT sig.
    # Row 2: padj not sig -> NOT sig.
    # Row 3: all three pass -> sig.
    df = pd.DataFrame({
        "padj":           [0.01, 0.01, 0.20, 0.01],
        "baseMean":       [50.0, 2.0,  50.0, 50.0],
        "log2FoldChange": [0.05, 0.90, 0.90, 0.50],
    })
    sig = config.is_significant(df)
    assert list(sig) == [False, False, False, True]


def test_is_significant_direction_uses_lfc_magnitude_and_sign():
    # |LFC| must exceed LFC_THRESH AND have the right sign.
    df = pd.DataFrame({
        "padj":           [0.01, 0.01, 0.01, 0.01],
        "baseMean":       [50.0, 50.0, 50.0, 50.0],
        "log2FoldChange": [0.50, -0.50, 0.05, -0.05],
    })
    # small-|LFC| rows (0.05) fail the magnitude floor in both directions
    assert list(config.is_significant(df, direction="UP"))   == [True, False, False, False]
    assert list(config.is_significant(df, direction="DOWN")) == [False, True, False, False]


import numpy as np, anndata
from pathlib import Path
from tusc2_deg.cd8 import degs


def _toy_pseudobulk():
    rng = np.random.default_rng(0)
    donors = [f"d{i}" for i in range(12)]
    obs = []; rows = []
    for d in donors:
        for st in ("neg", "pos"):
            obs.append({"donor": d, "tusc2_status": st})
            base = rng.poisson(50, 40).astype(float)
            if st == "neg":
                base[0] *= 1.3
            rows.append(base)
    X = np.vstack(rows)
    var = pd.DataFrame(index=[f"g{j}" for j in range(40)])
    ad = anndata.AnnData(X=X, obs=pd.DataFrame(obs), var=var)
    ad.obs_names = [f"s{i}" for i in range(len(obs))]
    return ad


def test_run_deseq2_emits_standard_and_calibrated_padj(tmp_path):
    pb = _toy_pseudobulk()
    res = degs.run_deseq2(pb, "CD8_TIL", out_path=tmp_path / "deg.csv", n_cpus=1)
    assert "stat" in res.columns
    assert "padj" in res.columns and "padj_calibrated" in res.columns
    assert "padj_lfc0" not in res.columns
    s_std = (res["padj"] < 0.05).sum()             # standard (canonical) gate
    s_cal = (res["padj_calibrated"] < 0.05).sum()  # calibrated disclosure
    assert s_cal <= s_std
