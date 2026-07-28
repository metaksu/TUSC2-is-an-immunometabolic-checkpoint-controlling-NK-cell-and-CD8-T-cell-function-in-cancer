"""GSEA prerank (threshold-free pathway enrichment).

Ranks the full transcriptome by the signed Wald statistic (column `stat`,
neg-vs-pos so positive = up in TUSC2-negative) and runs gseapy.prerank against
the configured gene-set libraries. Being threshold-free, GSEA sidesteps the
|LFC| gate and foreground-size choices required by list-based enrichment
(Subramanian 2005 PNAS; sc-best-practices Ch.18).

The gp.prerank call is injectable (prerank_fn) so the ranking and I/O logic is
unit-testable offline without a network call.
"""
from __future__ import annotations
import re
from pathlib import Path
import pandas as pd

from tusc2_deg.cd8.config import ENRICHR_LIBS

RANK_METRIC = "stat"  # signed Wald z (neg-vs-pos): positive = up in TUSC2-neg
_GO_ID_RE = re.compile(r"\s*\(GO:\d+\)\s*$")


def _norm_term(t: str) -> str:
    s = str(t)
    if "__" in s:                       # gseapy prefixes terms with "Library__"
        s = s.split("__", 1)[1]
    return _GO_ID_RE.sub("", s).strip()  # drop trailing "(GO:nnnn)"


def build_rank(deg_df: pd.DataFrame) -> pd.Series:
    """gene -> signed Wald stat: NaN-dropped, deduped (keep first), descending."""
    df = deg_df[["gene", RANK_METRIC]].dropna(subset=[RANK_METRIC])
    df = df.drop_duplicates("gene", keep="first")
    s = df.set_index("gene")[RANK_METRIC].astype(float)
    return s.sort_values(ascending=False)


def run_prerank(deg_df: pd.DataFrame, ctx: str, out_path: Path | str,
                libraries: list[str] | None = None, *,
                permutation_num: int = 1000, min_size: int = 15,
                max_size: int = 500, seed: int = 0, prerank_fn=None) -> pd.DataFrame:
    """Run gseapy.prerank on the signed-stat ranking; write res2d; return it."""
    libraries = list(libraries) if libraries is not None else list(ENRICHR_LIBS)
    rnk = build_rank(deg_df).reset_index()  # 2-col (gene, stat) for gseapy
    if prerank_fn is None:
        import gseapy as gp
        prerank_fn = gp.prerank
    res = prerank_fn(rnk=rnk, gene_sets=libraries, min_size=min_size,
                     max_size=max_size, permutation_num=permutation_num,
                     seed=seed, no_plot=True, outdir=None, verbose=False)
    df = res.res2d.copy() if hasattr(res, "res2d") else pd.DataFrame(res)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    print(f"[{ctx}] GSEA prerank wrote {out_path}: rows={len(df)}", flush=True)
    return df
