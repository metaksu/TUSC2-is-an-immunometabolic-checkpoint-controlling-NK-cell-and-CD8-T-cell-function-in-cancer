"""huARdb v2 CD8 atlas (Xue 2025 Nat Methods, Zenodo 12542577) loader
+ pseudobulk builder.

Tumor-only build: the CD8 context is CD8_TIL (solid-tumor TIL).

The huARdb h5ad ships .X as integer raw counts (float32-typed but
fractionally integer) and has no .raw layer. After subsetting per context
we mirror .X into .raw so the downstream code path (which reads TUSC2 from
adata.raw) works unchanged.

build_pseudobulk: AnnData -> AnnData. Sums raw counts per (donor,
tusc2_status) into pseudosamples; drops donors with < MIN_CELLS_PER_GROUP
cells in either group. TUSC2 is asserted present after gene filtering —
if it isn't, the pseudobulk is too sparse for the analysis.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from scipy.sparse import issparse, csr_matrix

from tusc2_deg.paths import HUARDB_CD8_ALLGENES_H5AD
from tusc2_deg.cd8.config import MIN_CELLS_PER_GROUP, CONTEXTS, GENE_SUM_FLOOR

TUSC2_THRESHOLD = 0
DONOR_COL = "individual_id"  # 526 unique IDs, present on every cell

# Per-context boolean masks over obs. Each value is a callable
# obs -> pd.Series[bool] so we can re-evaluate against any subset.
CONTEXT_FILTERS = {
    "CD8_TIL": lambda obs: (obs["disease_type"] == "Solid tumor")
                            & (obs["meta_tissue_type"] == "TIL"),
}


def _read_tusc2(adata: anndata.AnnData) -> np.ndarray:
    """Read TUSC2 expression from .raw (preferred) or .X."""
    if adata.raw is not None and "TUSC2" in adata.raw.var_names:
        x = adata.raw[:, "TUSC2"].X
    elif "TUSC2" in adata.var_names:
        x = adata[:, "TUSC2"].X
    else:
        raise RuntimeError("TUSC2 missing from adata — cannot annotate status.")
    expr = np.asarray(x.todense() if hasattr(x, "todense") else x).flatten()
    # NaN would be silently classed neg (NaN > 0 is False); fail loud instead.
    assert not np.isnan(expr).any(), "TUSC2 expression contains NaN — cannot classify tusc2_status"
    return expr


def _annotate(adata: anndata.AnnData, donor_col: str) -> anndata.AnnData:
    """Tag each cell with tusc2_status and remember donor_col on adata.uns."""
    adata = adata.copy()
    expr = _read_tusc2(adata)
    adata.obs["tusc2_status"] = np.where(expr > TUSC2_THRESHOLD, "pos", "neg")
    adata.obs["tusc2_status"] = adata.obs["tusc2_status"].astype("category")
    adata.uns["donor_col"] = donor_col
    return adata


def load_all_contexts() -> dict[str, anndata.AnnData]:
    """Return dict ctx -> AnnData (subset, raw mirrored, status annotated).

    Reads the full huARdb h5ad once into memory, then subsets per context.
    Memory: ~12 GB peak (1.13M cells × 19957 genes float32 sparse). The
    pipeline runs once per context; loading once and subsetting is faster
    than reloading per call.
    """
    print(f"Loading huARdb v2 CD8 atlas: {HUARDB_CD8_ALLGENES_H5AD}", flush=True)
    full = sc.read_h5ad(HUARDB_CD8_ALLGENES_H5AD)
    print(f"Loaded {full.n_obs:,} cells × {full.n_vars:,} genes", flush=True)
    if "TUSC2" not in full.var_names:
        raise RuntimeError("TUSC2 missing from huARdb var_names — wrong file?")

    # Mirror .X -> .raw so the downstream code path (which reads TUSC2
    # from adata.raw) works unchanged. .X is integer raw counts already.
    full.raw = full

    out: dict[str, anndata.AnnData] = {}
    for ctx in CONTEXTS:
        if ctx not in CONTEXT_FILTERS:
            raise KeyError(f"No CONTEXT_FILTERS entry for {ctx}")
        raw_mask = CONTEXT_FILTERS[ctx](full.obs)
        n_dropped = int(raw_mask.isna().sum())
        mask = raw_mask.fillna(False).values
        if n_dropped:
            print(f"  [{ctx}] {n_dropped:,} cells dropped by NaN context mask "
                  f"(fillna(False))", flush=True)
        n_cells = int(mask.sum())
        if n_cells == 0:
            raise RuntimeError(f"[{ctx}] zero cells after filtering")
        sub = full[mask, :].copy()
        n_donors = sub.obs[DONOR_COL].nunique()
        print(f"  [{ctx}] {n_cells:,} cells, {n_donors} donors", flush=True)
        out[ctx] = _annotate(sub, donor_col=DONOR_COL)
    return out


def build_pseudobulk(adata: anndata.AnnData, ctx: str) -> anndata.AnnData:
    """Sum raw counts per (donor, tusc2_status) into pseudosamples.

    Drops donors that don't have >= MIN_CELLS_PER_GROUP cells in BOTH
    tusc2_status groups. Returns AnnData with obs = (donor, tusc2_status,
    n_cells), var indexed by gene name, .X = integer count matrix (CSR).
    Asserts TUSC2 survived the gene-sum filter — if it didn't, the
    pseudobulk is too sparse to interpret.
    """
    donor_col = adata.uns["donor_col"]
    counts = (
        adata.obs[[donor_col, "tusc2_status"]]
        .groupby([donor_col, "tusc2_status"], observed=True)
        .size().unstack(fill_value=0)
    )
    keep_donors = counts.index[
        (counts.get("pos", 0) >= MIN_CELLS_PER_GROUP)
        & (counts.get("neg", 0) >= MIN_CELLS_PER_GROUP)
    ].tolist()
    if not keep_donors:
        raise RuntimeError(
            f"[{ctx}] no donors satisfy MIN_CELLS_PER_GROUP={MIN_CELLS_PER_GROUP} "
            f"in both tusc2_status groups"
        )

    var_names = adata.raw.var_names if adata.raw is not None else adata.var_names
    X = adata.raw.X if adata.raw is not None else adata.X
    cell_donor = adata.obs[donor_col].values
    cell_status = adata.obs["tusc2_status"].values

    rows, meta_rows = [], []
    for d in keep_donors:
        for s in ("pos", "neg"):
            mask = (cell_donor == d) & (cell_status == s)
            # Defensive/unreachable: keep_donors already guarantees both
            # arms clear MIN_CELLS_PER_GROUP and nothing mutates obs in
            # between, so this never fires on current data. Kept as a guard
            # in case that upstream invariant ever changes.
            if mask.sum() < MIN_CELLS_PER_GROUP:
                continue
            sub = X[mask, :]
            arr = sub.toarray() if issparse(sub) else np.asarray(sub)
            # The source .raw may carry decontaminated/fractional per-cell counts
            # (the Tang atlas stores non-integer values for a subset of donors);
            # round per cell to the integer counts DESeq2 expects.
            arr = np.rint(arr).astype(np.int64)
            rows.append(arr.sum(axis=0))
            meta_rows.append(dict(sample=f"{d}__{s}", donor=str(d),
                                  tusc2_status=s, n_cells=int(mask.sum())))

    assert var_names.is_unique, (
        f"[{ctx}] var_names are not unique — cannot build counts_df safely"
    )
    counts_df = pd.DataFrame(np.vstack(rows),
                             index=[m["sample"] for m in meta_rows],
                             columns=var_names)
    keep_genes = counts_df.sum(axis=0) >= GENE_SUM_FLOOR
    counts_df = counts_df.loc[:, keep_genes]
    if "TUSC2" not in counts_df.columns:
        raise RuntimeError(
            f"[{ctx}] TUSC2 dropped by sum>={GENE_SUM_FLOOR} gene filter — pseudobulk too sparse"
        )

    pb = anndata.AnnData(
        X=csr_matrix(counts_df.values.astype(np.int64)),
        obs=pd.DataFrame(meta_rows).set_index("sample"),
        var=pd.DataFrame(index=counts_df.columns),
    )
    return pb
