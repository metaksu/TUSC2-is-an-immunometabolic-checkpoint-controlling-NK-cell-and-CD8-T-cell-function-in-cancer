"""DESeq2 (primary) and Wilcoxon (supplementary) DEGs per context.

run_deseq2 takes a pseudobulk AnnData (built by data.build_pseudobulk) rather
than building the pseudobulk internally — data and degs are separate concerns.
Thresholds are imported from config.

DESeq2 column convention:
    log2FoldChange         = apeGLM-shrunken LFC (post-lfc_shrink)
    log2FoldChange_unshrunk = unshrunken Wald LFC (saved before lfc_shrink)
    padj                   = BH-adjusted p-value from Wald test
"""
# Parallel to cd8/degs.py; the two differ only by their config import.
from __future__ import annotations
from pathlib import Path
import datetime
import numpy as np
import pandas as pd
import anndata
import scanpy as sc
from scipy.sparse import issparse
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds  import DeseqStats
from pydeseq2.default_inference import DefaultInference

from tusc2_deg.nk.config import (
    EXCLUDE_PREFIXES, GROUPING_GENES, LFC_THRESH, PADJ_THRESH, SIG_LFC_LADDER,
    is_significant,
)


def _is_excluded(g: str) -> bool:
    return g in GROUPING_GENES or any(g.startswith(p) for p in EXCLUDE_PREFIXES)


def _find_tusc2_coeff(ds: DeseqStats) -> str:
    """pydeseq2 0.5.x names the coeff 'tusc2_status[T.pos]' (alphabetical
    default; ref_level is deprecated/no-op). LFC reported under that coeff
    is pos-vs-neg; we negate post-shrinkage in run_deseq2 to flip to
    neg-vs-pos (absence-of-TUSC2 framing)."""
    candidates = ["tusc2_status_pos_vs_neg", "tusc2_status[T.pos]"]
    for c in candidates:
        if c in ds.LFC.columns:
            return c
    raise KeyError(
        f"Cannot find tusc2_status coeff in LFC columns. "
        f"Available: {list(ds.LFC.columns)}. Tried: {candidates}"
    )


def run_deseq2(pb: anndata.AnnData, ctx: str, out_path: Path | str,
               n_cpus: int = 8) -> pd.DataFrame:
    """Run DESeq2 with formula '~donor + tusc2_status' on pseudobulk AnnData."""
    X = pb.X.toarray() if issparse(pb.X) else np.asarray(pb.X)
    counts_df = pd.DataFrame(X.astype(np.int64), index=pb.obs_names,
                             columns=pb.var_names)
    meta_df = pb.obs[["donor", "tusc2_status"]].copy()

    print(f"[{ctx}] DESeq2 input: samples={counts_df.shape[0]} genes={counts_df.shape[1]} "
          f"donors={meta_df['donor'].nunique()}", flush=True)

    inference = DefaultInference(n_cpus=n_cpus)
    dds = DeseqDataSet(
        counts=counts_df, metadata=meta_df,
        design="~donor + tusc2_status",
        refit_cooks=True, inference=inference, quiet=True,
    )
    dds.deseq2()

    from tusc2_deg.nk.config import LFC_NULL, ALT_HYPOTHESIS

    # --- standard test (H0: LFC=0): signed Wald `stat` for GSEA + unshrunk LFC ---
    ds = DeseqStats(dds, contrast=["tusc2_status", "neg", "pos"], alpha=0.05)
    ds.summary()
    wald_lfc  = ds.results_df["log2FoldChange"].copy()
    std_stat  = ds.results_df["stat"].copy()
    std_padj  = ds.results_df["padj"].copy()
    coeff = _find_tusc2_coeff(ds)
    ds.lfc_shrink(coeff=coeff)
    res = ds.results_df.copy()
    # lfc_shrink writes pos-vs-neg shrunken values back into log2FoldChange
    # (the model's only fitted coeff is 'tusc2_status[T.pos]'). Flip its sign to
    # match the contrast orientation (neg-vs-pos = absence-of-TUSC2 framing).
    # Invariant: 'stat'/'pvalue'/'padj' are NOT overwritten by lfc_shrink — they
    # remain the neg-vs-pos Wald-test values from
    # DeseqStats(contrast=[.., 'neg', 'pos']), so 'stat' is already neg-vs-pos
    # and must NOT be flipped (do not negate it). We restore the standard `stat`
    # explicitly below to keep it aligned with log2FoldChange and GSEA prerank.
    res["log2FoldChange"] = -res["log2FoldChange"]
    res["log2FoldChange_unshrunk"] = wald_lfc.reindex(res.index)
    res["stat"]      = std_stat.reindex(res.index)
    res["padj"]      = std_padj.reindex(res.index)   # STANDARD (H0:LFC=0) BH padj = canonical gate

    # --- calibrated lfcThreshold test (H0: |LFC|<=LFC_NULL): disclosure column ---
    # The Wald test against a non-zero null with alt_hypothesis 'greaterAbs' bakes
    # the fold-change condition into padj (Love 2014). It is not the significance
    # gate; it is exposed as a per-gene reference column `padj_calibrated` so
    # Methods can report how many genes exceed the calibrated lfcThreshold null.
    # The canonical `padj` above stays the standard Wald padj; `stat` stays the
    # standard signed Wald statistic (GSEA rank metric).
    ds_cal = DeseqStats(dds, contrast=["tusc2_status", "neg", "pos"], alpha=0.05,
                        lfc_null=LFC_NULL, alt_hypothesis=ALT_HYPOTHESIS)
    ds_cal.summary()
    res["padj_calibrated"] = ds_cal.results_df["padj"].reindex(res.index)

    res["gene"] = res.index.astype(str)
    res = res[~res["gene"].apply(_is_excluded)].reset_index(drop=True)

    cols_ordered = ["gene", "baseMean", "log2FoldChange", "log2FoldChange_unshrunk",
                    "lfcSE", "stat", "pvalue", "padj", "padj_calibrated"]
    cols_present = [c for c in cols_ordered if c in res.columns]
    cols_extra   = [c for c in res.columns if c not in cols_present]
    res = res[cols_present + cols_extra]

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(out_path, index=False)
    print(f"[{ctx}] wrote {out_path}", flush=True)
    return res


def run_wilcoxon(adata: anndata.AnnData, ctx: str,
                 out_path: Path | str) -> pd.DataFrame:
    """scanpy Wilcoxon DEGs on cell-level data (supplementary).

    Schema: gene, scores, logfoldchanges, pvals, pvals_adj, pct_pos, pct_neg.
    """
    sc.tl.rank_genes_groups(
        adata, groupby="tusc2_status", groups=["neg"], reference="pos",
        method="wilcoxon", use_raw=True, corr_method="benjamini-hochberg",
        pts=True,
        n_genes=adata.raw.shape[1] if adata.raw is not None else adata.shape[1],
        key_added="degs_t",
    )
    df = sc.get.rank_genes_groups_df(adata, group="neg", key="degs_t").rename(
        columns={"names": "gene", "pct_nz_group": "pct_neg"}
    )
    pts = adata.uns["degs_t"]["pts"]
    if "pos" not in pts.columns:
        raise RuntimeError(
            f"[{ctx}] Expected 'pos' column in adata.uns['degs_t']['pts']; "
            f"found columns={list(pts.columns)}"
        )
    pct_pos = pts["pos"].rename_axis("gene").reset_index().rename(columns={"pos": "pct_pos"})
    df = df.merge(pct_pos, on="gene", how="left")
    df = df[~df["gene"].apply(_is_excluded)].reset_index(drop=True)

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    print(f"[{ctx}] wrote {out_path}: rows={len(df)}", flush=True)
    return df


CONCORDANCE_HEADER = ("timestamp\tcontext\tspearman_rho\tn_pb_sig\t"
                      "co_significance_rate\tsign_agreement_rate\tsign_agreement_hi\n")


def write_concordance_audit(deseq_df: pd.DataFrame, wilc_df: pd.DataFrame,
                            ctx: str, audit_path: Path | str) -> None:
    """Append one row per context to concordance_pb_vs_wilcoxon.tsv.

    Columns:
      spearman_rho          Spearman corr of pb vs Wilcoxon LFC (all overlapping genes)
      n_pb_sig              # pb-significant genes (canonical 3-part is_significant gate)
      co_significance_rate  fraction of pb-sig genes ALSO Wilcoxon-significant.
                            This is co-occurrence of significance, NOT a direction
                            check.
      sign_agreement_rate   fraction of pb-sig genes whose pb LFC sign matches the
                            Wilcoxon rank-score sign — the real directional metric.
      sign_agreement_hi     sign_agreement restricted to pb-sig genes with
                            |pb LFC| > 0.585 (~1.5x); NaN if none.

    Sign uses the Wilcoxon rank 'scores' (z-statistic), not 'logfoldchanges'
    (which carries a reference-relative log1p offset). Pseudobulk is primary;
    Wilcoxon is a confirmatory directional check (Squair 2021).
    """
    audit_path = Path(audit_path)
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    merged = deseq_df.merge(wilc_df, on="gene", suffixes=("_pb", "_wc"))
    if merged.empty:
        rho = co_sig = sign_rate = sign_hi = float("nan")
        n_pb_sig = 0
    else:
        from scipy.stats import spearmanr
        rho = spearmanr(merged["log2FoldChange"], merged["logfoldchanges"]).correlation
        pb_sig = is_significant(merged)
        wc_sig = (merged["pvals_adj"] < PADJ_THRESH) & (merged["logfoldchanges"].abs() > LFC_THRESH)
        n_pb_sig = int(pb_sig.sum())
        co_sig = float((pb_sig & wc_sig).sum() / max(1, n_pb_sig))
        sig = merged[pb_sig]
        if len(sig):
            agree = np.sign(sig["log2FoldChange"]) == np.sign(sig["scores"])
            sign_rate = float(agree.mean())
            hi = sig["log2FoldChange"].abs() > SIG_LFC_LADDER[-1]  # high-effect tail (~1.5x)
            sign_hi = float(agree[hi].mean()) if hi.any() else float("nan")
        else:
            sign_rate = sign_hi = float("nan")

    # Schema-aware: if an existing file has a header that does not match the
    # current schema, start fresh so the audit is internally consistent.
    if audit_path.exists():
        first_line = audit_path.read_text(encoding="utf-8").splitlines()[:1]
        if not first_line or (first_line[0] + "\n") != CONCORDANCE_HEADER:
            audit_path.unlink()
    # De-dup-on-write: keep the file bounded to the latest row per context.
    # If a schema-current file already exists, rewrite header + only the data rows
    # for OTHER contexts (context is column index 1); the new row for this ctx is
    # appended below. This is a human QC log with no reader, so it trades run
    # history for a bounded, one-row-per-context file — the per-row timestamp still
    # records the surviving row's provenance.
    if audit_path.exists():
        existing = audit_path.read_text(encoding="utf-8").splitlines()
        kept = [ln for ln in existing[1:] if ln and ln.split("\t")[1:2] != [ctx]]
        audit_path.write_text(
            CONCORDANCE_HEADER + ("\n".join(kept) + "\n" if kept else ""),
            encoding="utf-8",
        )
    else:
        audit_path.write_text(CONCORDANCE_HEADER, encoding="utf-8")
    ts = datetime.datetime.now(datetime.timezone.utc).isoformat()
    with audit_path.open("a", encoding="utf-8") as f:
        f.write(f"{ts}\t{ctx}\t{rho:.4f}\t{n_pb_sig}\t{co_sig:.4f}\t"
                f"{sign_rate:.4f}\t{sign_hi:.4f}\n")
