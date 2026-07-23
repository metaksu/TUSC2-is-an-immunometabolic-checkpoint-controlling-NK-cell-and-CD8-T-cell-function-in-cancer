"""degs.py contracts: pseudobulk AnnData -> DESeq2 DataFrame."""
import numpy as np
import pandas as pd
import anndata
from scipy.sparse import csr_matrix

from tusc2_deg.cd8 import degs


def _make_toy_pseudobulk(n_donors: int = 6, n_genes: int = 50) -> anndata.AnnData:
    rng = np.random.default_rng(7)
    samples = []
    for d in range(n_donors):
        for s in ("pos", "neg"):
            samples.append(dict(sample=f"d{d}__{s}", donor=f"d{d}",
                                tusc2_status=s, n_cells=50))
    X = rng.poisson(20, size=(len(samples), n_genes)).astype(np.int64)
    var = pd.DataFrame(index=[f"GENE{i}" for i in range(n_genes - 1)] + ["TUSC2"])
    obs = pd.DataFrame(samples).set_index("sample")
    return anndata.AnnData(X=csr_matrix(X), obs=obs, var=var)


def test_run_deseq2_schema(tmp_path):
    pb = _make_toy_pseudobulk()
    out = tmp_path / "deg.csv"
    df = degs.run_deseq2(pb, "CD8_TIL", out_path=out)
    expected_cols = {"gene", "baseMean", "log2FoldChange",
                     "log2FoldChange_unshrunk", "lfcSE", "stat", "pvalue", "padj"}
    assert expected_cols.issubset(df.columns)
    assert out.exists()


def test_run_deseq2_excludes_ribosomal_but_keeps_mt(tmp_path):
    """RPL/RPS are dropped (EXCLUDE_PREFIXES); MT- is intentionally KEPT
    (mtDNA-encoded OXPHOS genes are the most direct respiratory-chain readout).
    Inject one of each and verify."""
    pb = _make_toy_pseudobulk()
    # inject an RPL gene and an MT- gene (replace two non-TUSC2 genes)
    idx = list(pb.var.index)
    idx[0] = "RPL13"
    idx[1] = "MT-ND1"
    pb.var.index = idx
    out = tmp_path / "deg.csv"
    df = degs.run_deseq2(pb, "CD8_TIL", out_path=out)
    genes = set(df["gene"])
    assert not df["gene"].str.startswith(("RPL", "RPS")).any()
    assert "MT-ND1" in genes, "MT- genes must survive (kept in the matrix)"


def test_run_deseq2_excludes_grouping_gene_tusc2(tmp_path):
    """TUSC2 (the grouping variable) must not appear in the DEG table."""
    pb = _make_toy_pseudobulk()
    assert "TUSC2" in pb.var_names
    out = tmp_path / "deg.csv"
    df = degs.run_deseq2(pb, "CD8_TIL", out_path=out)
    assert "TUSC2" not in set(df["gene"])


def test_run_deseq2_stat_sign_matches_wald_lfc(tmp_path):
    """The Wald 'stat' must share sign with the unshrunken Wald log2FC (stat =
    LFC/SE). This guards GSEA prerank rankings, which require 'stat' and the
    unshrunken LFC to agree in sign."""
    pb = _make_toy_pseudobulk()
    df = degs.run_deseq2(pb, "CD8_TIL", out_path=tmp_path / "deg.csv")
    m = df.dropna(subset=["stat", "log2FoldChange_unshrunk"])
    m = m[m["log2FoldChange_unshrunk"].abs() > 1e-9]
    same = np.sign(m["stat"]) == np.sign(m["log2FoldChange_unshrunk"])
    assert same.all(), f"stat sign must match unshrunken Wald LFC; {(~same).sum()}/{len(m)} mismatch"


def _make_toy_pseudobulk_pos_high(target_gene: str = "GENE0",
                                  n_donors: int = 6, n_genes: int = 50) -> anndata.AnnData:
    """Clone of _make_toy_pseudobulk with a real directional signal: `target_gene`
    is driven ~5x high (100) in every 'pos' sample, baseline poisson(20) elsewhere.
    Under the neg-vs-pos contrast a gene genuinely up in 'pos' must read NEGATIVE."""
    pb = _make_toy_pseudobulk(n_donors=n_donors, n_genes=n_genes)
    X = pb.X.toarray()
    gi = list(pb.var_names).index(target_gene)
    pos = (pb.obs["tusc2_status"] == "pos").to_numpy()
    X[pos, gi] = 100  # ~5x the poisson(20) baseline; up in 'pos'
    pb.X = csr_matrix(X.astype(np.int64))
    return pb


def test_run_deseq2_orientation_neg_vs_pos_contrast(tmp_path):
    """Validates the sign ORIENTATION, not trivial sign consistency. A gene
    genuinely higher in the 'pos' group must read NEGATIVE in all three signed
    columns (stat, log2FoldChange, log2FoldChange_unshrunk) under the
    contrast=['tusc2_status','neg','pos'] with the post-shrink log2FoldChange flip.
    This locks the absence-of-TUSC2 framing (+values = higher in TUSC2-negative).
    Unlike test_run_deseq2_stat_sign_matches_wald_lfc (which only checks that stat
    and the LFC agree in sign — automatic since stat=LFC/SE on a signal-free
    fixture), this drives a ~5x signal through the full run_deseq2 path; apeGLM
    shrinkage attenuates magnitude, not sign. Passes with unchanged degs.py."""
    pb = _make_toy_pseudobulk_pos_high("GENE0")
    df = degs.run_deseq2(pb, "CD8_TIL", out_path=tmp_path / "deg.csv")
    row = df[df["gene"] == "GENE0"]
    assert len(row) == 1, "target gene must survive into the DEG table"
    assert row["stat"].iloc[0] < 0, "gene up in 'pos' -> stat NEGATIVE (neg-vs-pos)"
    assert row["log2FoldChange"].iloc[0] < 0, "shrunken LFC NEGATIVE (neg-vs-pos)"
    assert row["log2FoldChange_unshrunk"].iloc[0] < 0, "unshrunk Wald LFC NEGATIVE"


def _make_toy_singlecell(n_cells: int = 200, n_genes: int = 50):
    rng = np.random.default_rng(11)
    X = rng.poisson(3, size=(n_cells, n_genes)).astype(np.float32)
    var = pd.DataFrame(index=[f"GENE{i}" for i in range(n_genes - 1)] + ["TUSC2"])
    obs = pd.DataFrame({
        "tusc2_status": pd.Categorical(["pos"] * (n_cells // 2) + ["neg"] * (n_cells - n_cells // 2)),
    })
    ad = anndata.AnnData(X=X, obs=obs, var=var)
    ad.raw = ad.copy()
    return ad


def test_run_wilcoxon_schema(tmp_path):
    ad = _make_toy_singlecell()
    out = tmp_path / "wilcoxon.csv"
    df = degs.run_wilcoxon(ad, "CD8_TIL", out_path=out)
    expected_cols = {"gene", "scores", "logfoldchanges", "pvals", "pvals_adj",
                     "pct_pos", "pct_neg"}
    assert expected_cols.issubset(df.columns)
    assert out.exists()


def test_concordance_audit_computes_correct_metrics(tmp_path):
    """Hand-crafted toy contract for the concordance audit metrics.
    n_pb_sig=4 (A,B,D,G); co_significance=0.5; sign_agreement=0.75; hi-band=1.0.
    Each metric is derived per-gene from the toy DESeq2 and Wilcoxon tables below."""
    deseq = pd.DataFrame({
        "gene":           ["A",  "B",   "C",  "D",   "E",  "F",   "G"],
        "log2FoldChange": [0.2, -0.15, 0.05, -0.12, 0.03, 0.9,   0.7],
        "padj":           [0.01, 0.04, 0.50,  0.03, 0.80, 0.001, 0.001],
        "baseMean":       [100,  100,  100,   100,  100,  2,     100],
    })
    wilc = pd.DataFrame({
        "gene":           ["A",  "B",   "C",  "D",   "E",  "F",  "G"],
        "logfoldchanges": [0.18, -0.14, 0.04, -0.08, 0.02, 0.5,  0.6],
        "scores":         [3.0, -2.0,   0.5,  1.0,   0.1,  4.0,  5.0],
        "pvals_adj":      [0.02, 0.06,  0.60, 0.04,  0.90, 0.001, 0.001],
    })

    audit = tmp_path / "concordance.tsv"
    degs.write_concordance_audit(deseq, wilc, "test_deterministic", audit_path=audit)

    rows = audit.read_text(encoding="utf-8").splitlines()
    assert len(rows) == 2
    header = rows[0].split("\t")
    assert header == ["timestamp", "context", "spearman_rho", "n_pb_sig",
                      "co_significance_rate", "sign_agreement_rate", "sign_agreement_hi"]
    fields = rows[1].split("\t")
    assert int(fields[3]) == 4
    assert abs(float(fields[4]) - 0.5) < 1e-4   # co_significance
    assert abs(float(fields[5]) - 0.75) < 1e-4  # sign_agreement
    assert abs(float(fields[6]) - 1.0) < 1e-4   # hi-band sign agreement
