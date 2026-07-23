import pandas as pd
from tusc2_deg.publication import data

def test_load_pseudobulk_nk_has_expected_columns():
    df = data.load_pseudobulk("nk")
    for col in ("gene", "baseMean", "log2FoldChange", "stat", "padj"):
        assert col in df.columns
    assert len(df) > 10000

def test_load_pseudobulk_cd8_has_more_degs_than_nk():
    nk = data.load_pseudobulk("nk")
    cd8 = data.load_pseudobulk("cd8")
    # CD8 has a far larger signal than NK (3652 vs 123 at |LFC|>0.1)
    assert (cd8["padj"] < 0.05).sum() > (nk["padj"] < 0.05).sum()

def test_load_gsea_nk_has_nes_and_fdr():
    g = data.load_gsea("nk")
    assert {"Term", "NES", "FDR q-val"} <= set(g.columns)

def test_load_wilcoxon_nk_has_logfoldchanges():
    w = data.load_wilcoxon("nk")
    assert {"gene", "logfoldchanges", "pvals_adj"} <= set(w.columns)

def test_load_sig_ladder_cd8():
    lad = data.load_sig_ladder("cd8")
    assert {"context", "lfc_cutoff", "UP", "DOWN"} <= set(lad.columns)
