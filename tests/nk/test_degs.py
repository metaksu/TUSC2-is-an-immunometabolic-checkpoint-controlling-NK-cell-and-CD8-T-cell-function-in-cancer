"""degs.py contracts: pseudobulk AnnData -> DESeq2 DataFrame."""
import numpy as np
import pandas as pd
import anndata
import scanpy as sc
from scipy.sparse import csr_matrix
from pathlib import Path

from tusc2_deg.nk import degs


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
    df = degs.run_deseq2(pb, "Tang_Tumor", out_path=out)
    expected_cols = {"gene", "baseMean", "log2FoldChange",
                     "log2FoldChange_unshrunk", "lfcSE", "stat", "pvalue", "padj"}
    assert expected_cols.issubset(df.columns)
    assert out.exists()


def test_run_deseq2_excludes_ribosomal_but_keeps_mt(tmp_path):
    """RPL/RPS are excluded (EXCLUDE_PREFIXES); MT- is intentionally KEPT
    (mtDNA-encoded OXPHOS genes are the most direct respiratory-chain readout).
    Inject one of each and verify."""
    pb = _make_toy_pseudobulk()
    # inject an RPL gene and an MT- gene (replace two non-TUSC2 genes)
    idx = list(pb.var.index)
    idx[0] = "RPL13"
    idx[1] = "MT-ND1"
    pb.var.index = idx
    out = tmp_path / "deg.csv"
    df = degs.run_deseq2(pb, "Tang_Tumor", out_path=out)
    genes = set(df["gene"])
    assert not df["gene"].str.startswith(("RPL", "RPS")).any()
    assert "MT-ND1" in genes, "MT- genes must survive (kept in the matrix)"


def test_run_deseq2_excludes_grouping_gene_tusc2(tmp_path):
    """TUSC2 (the grouping variable) must not appear in the DEG table — it is
    quasi-complete-separated and would seed circular enrichment."""
    pb = _make_toy_pseudobulk()
    assert "TUSC2" in pb.var_names  # present in input
    out = tmp_path / "deg.csv"
    df = degs.run_deseq2(pb, "Tang_Tumor", out_path=out)
    assert "TUSC2" not in set(df["gene"])


def test_run_deseq2_stat_sign_matches_wald_lfc(tmp_path):
    """The Wald 'stat' must share sign with the unshrunken Wald log2FC — both are
    the same neg-vs-pos test (stat = LFC/SE, SE>0). Guards against any inconsistency
    where 'stat' and the LFC disagree in sign, which would invert GSEA ranks."""
    pb = _make_toy_pseudobulk()
    df = degs.run_deseq2(pb, "Tang_Tumor", out_path=tmp_path / "deg.csv")
    m = df.dropna(subset=["stat", "log2FoldChange_unshrunk"])
    m = m[m["log2FoldChange_unshrunk"].abs() > 1e-9]
    same = np.sign(m["stat"]) == np.sign(m["log2FoldChange_unshrunk"])
    assert same.all(), f"stat sign must match unshrunken Wald LFC; {(~same).sum()}/{len(m)} mismatch"


def _make_toy_pseudobulk_directional(signal_gene: str = "GENE0", pos_level: int = 100,
                                     n_donors: int = 6, n_genes: int = 50) -> anndata.AnnData:
    """Clone of _make_toy_pseudobulk but with a genuine directional signal: one
    gene is forced HIGH (pos_level, ~5x the ~20 Poisson baseline) in every
    'pos'-group sample and left at baseline in every 'neg'-group sample. Used to
    validate the neg-vs-pos contrast ORIENTATION, which the signal-free fixture
    cannot check."""
    rng = np.random.default_rng(7)
    samples = []
    for d in range(n_donors):
        for s in ("pos", "neg"):
            samples.append(dict(sample=f"d{d}__{s}", donor=f"d{d}",
                                tusc2_status=s, n_cells=50))
    X = rng.poisson(20, size=(len(samples), n_genes)).astype(np.int64)
    var = pd.DataFrame(index=[f"GENE{i}" for i in range(n_genes - 1)] + ["TUSC2"])
    obs = pd.DataFrame(samples).set_index("sample")
    gi = list(var.index).index(signal_gene)
    pos_rows = (obs["tusc2_status"] == "pos").to_numpy()
    X[pos_rows, gi] = pos_level
    return anndata.AnnData(X=csr_matrix(X), obs=obs, var=var)


def test_run_deseq2_orientation_neg_vs_pos_contrast(tmp_path):
    """Orientation (NOT trivial sign consistency): a gene genuinely HIGHER in the
    'pos' group must come out NEGATIVE in the neg-vs-pos contrast — negative
    'stat', shrunken 'log2FoldChange', AND unshrunken 'log2FoldChange_unshrunk'.

    test_run_deseq2_stat_sign_matches_wald_lfc only checks that stat and LFC agree
    in sign with each other (automatic, since stat = LFC/SE, SE>0) on a signal-free
    fixture, so it never verifies the CONTRAST DIRECTION. This locks in that
    contrast=['tusc2_status','neg','pos'] plus the post-shrink sign flip (degs.py)
    yield the absence-of-TUSC2 (neg-vs-pos) framing: a gene up in pos is negative
    in all three columns. It validates the convention against unchanged degs.py;
    it does not fix a bug. Runs real pydeseq2 (integration-style, like the other
    run_deseq2 tests)."""
    signal_gene = "GENE0"
    pb = _make_toy_pseudobulk_directional(signal_gene=signal_gene)
    df = degs.run_deseq2(pb, "Tang_Tumor", out_path=tmp_path / "deg.csv")
    row = df[df["gene"] == signal_gene]
    assert len(row) == 1, f"signal gene {signal_gene} missing from DEG table"
    r = row.iloc[0]
    assert r["stat"] < 0, f"expected stat<0 for a gene up in pos, got {r['stat']}"
    assert r["log2FoldChange"] < 0, \
        f"expected shrunken log2FoldChange<0, got {r['log2FoldChange']}"
    assert r["log2FoldChange_unshrunk"] < 0, \
        f"expected unshrunken log2FoldChange_unshrunk<0, got {r['log2FoldChange_unshrunk']}"


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
    df = degs.run_wilcoxon(ad, "Tang_Tumor", out_path=out)
    expected_cols = {"gene", "scores", "logfoldchanges", "pvals", "pvals_adj",
                     "pct_pos", "pct_neg"}
    assert expected_cols.issubset(df.columns)
    assert out.exists()


def test_concordance_audit_computes_correct_metrics(tmp_path):
    """Deterministic test of the concordance audit metrics.

    Schema: timestamp, context, spearman_rho, n_pb_sig, co_significance_rate,
            sign_agreement_rate, sign_agreement_hi.

    Toy (canonical 3-part gate padj<0.05 & baseMean>5 & |LFC|>0.1):
      A  lfc=+0.2  padj=0.01 bM=100   -> pb-sig
      B  lfc=-0.15 padj=0.04 bM=100   -> pb-sig
      C  lfc=+0.05 padj=0.50 bM=100   -> not sig (lfc & padj)
      D  lfc=-0.12 padj=0.03 bM=100   -> pb-sig
      E  lfc=+0.03 padj=0.80 bM=100   -> not sig
      F  lfc=+0.9  padj=0.001 bM=2    -> NOT sig (baseMean<=5 — tests 3-part gate)
      G  lfc=+0.7  padj=0.001 bM=100  -> pb-sig AND in |LFC|>0.585 hi band
    => n_pb_sig = 4 (A, B, D, G)

    Wilcoxon co-significance (pvals_adj<0.05 & |logfoldchanges|>0.1):
      A True, B False(padj), D False(|lfc|=0.08), G True  => co-sig {A,G} = 2/4 = 0.5
    Sign agreement (sign(pb LFC) == sign(wilcoxon scores)) over pb-sig:
      A +/+ agree, B -/- agree, D -/+ DISAGREE, G +/+ agree => 3/4 = 0.75
    Hi-band sign agreement (|LFC|>0.585): {G} agrees => 1.0
    """
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
    assert len(rows) == 2  # header + one data row
    header = rows[0].split("\t")
    assert header == ["timestamp", "context", "spearman_rho", "n_pb_sig",
                      "co_significance_rate", "sign_agreement_rate", "sign_agreement_hi"]

    fields = rows[1].split("\t")
    n_pb_sig            = int(fields[3])
    co_significance     = float(fields[4])
    sign_agreement      = float(fields[5])
    sign_agreement_hi   = float(fields[6])

    assert n_pb_sig == 4, f"Expected 4 pb-sig genes (A,B,D,G), got {n_pb_sig}"
    assert abs(co_significance - 0.5) < 1e-4, f"Expected co-sig 0.5, got {co_significance}"
    assert abs(sign_agreement - 0.75) < 1e-4, f"Expected sign agreement 0.75, got {sign_agreement}"
    assert abs(sign_agreement_hi - 1.0) < 1e-4, f"Expected hi-band sign agreement 1.0, got {sign_agreement_hi}"


def test_concordance_audit_resets_stale_schema(tmp_path):
    """A pre-existing file with the old 5-column header is replaced, not appended."""
    audit = tmp_path / "concordance.tsv"
    audit.write_text("timestamp\tcontext\tspearman_rho\tn_pb_sig\tconcordance_rate\n"
                     "old\tctx\t0.1\t5\t0.5\n", encoding="utf-8")
    deseq = pd.DataFrame({"gene": ["A"], "log2FoldChange": [0.2], "padj": [0.01], "baseMean": [100]})
    wilc = pd.DataFrame({"gene": ["A"], "logfoldchanges": [0.18], "scores": [3.0], "pvals_adj": [0.02]})
    degs.write_concordance_audit(deseq, wilc, "ctx", audit_path=audit)
    rows = audit.read_text(encoding="utf-8").splitlines()
    assert rows[0].split("\t")[-1] == "sign_agreement_hi"  # new schema
    assert len(rows) == 2  # stale row gone, fresh header + one new row
