"""Tang NK atlas loader (Tang_Tumor context) + pseudobulk builder.

build_pseudobulk: AnnData -> AnnData. Sums raw counts per (donor, tusc2_status)
into pseudosamples; drops donors with < MIN_CELLS_PER_GROUP cells in either
group. TUSC2 is asserted present after gene filtering — if it isn't, the
pseudobulk is too sparse for the analysis.

load_all_contexts: loads each configured context AnnData with tusc2_status +
donor_col annotated.

TUSC2 status is read from adata.raw (see _read_tusc2); the Tang h5ad is assumed
to already provide .raw. build_pseudobulk integerizes the summed counts via
np.rint before handing them to DESeq2.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from scipy.sparse import issparse, csr_matrix

from tusc2_deg.paths import TANG_H5AD
from tusc2_deg.nk.config import CONTEXTS, MIN_CELLS_PER_GROUP, GENE_SUM_FLOOR

TANG_TISSUE_COL = "meta_tissue_in_paper"
TUSC2_THRESHOLD = 0
DONOR_COLS = {
    "Tang_Tumor":     "meta_patientID",
}
# tissue label (meta_tissue_in_paper value) per context
CONTEXT_TISSUE = {
    "Tang_Tumor":     "Tumor",
}


def _read_tusc2(adata: anndata.AnnData) -> np.ndarray:
    if adata.raw is None or "TUSC2" not in adata.raw.var_names:
        raise RuntimeError("TUSC2 missing from adata.raw — cannot annotate status.")
    x = adata.raw[:, "TUSC2"].X
    expr = np.asarray(x.todense() if hasattr(x, "todense") else x).flatten()
    # NaN would be silently classed neg (NaN > 0 is False); fail loud instead.
    assert not np.isnan(expr).any(), "TUSC2 expression contains NaN — cannot classify tusc2_status"
    return expr


def _annotate(adata: anndata.AnnData, donor_col: str) -> anndata.AnnData:
    adata = adata.copy()
    expr = _read_tusc2(adata)
    adata.obs["tusc2_status"] = np.where(expr > TUSC2_THRESHOLD, "pos", "neg")
    adata.obs["tusc2_status"] = adata.obs["tusc2_status"].astype("category")
    adata.uns["donor_col"] = donor_col
    return adata


def load_all_contexts() -> dict[str, anndata.AnnData]:
    """Return dict ctx -> AnnData with tusc2_status + donor_col annotated."""
    out: dict[str, anndata.AnnData] = {}
    print(f"Loading Tang h5ad: {TANG_H5AD}", flush=True)
    tang = sc.read_h5ad(TANG_H5AD)
    if TANG_TISSUE_COL not in tang.obs.columns:
        raise RuntimeError(f"Tang missing tissue column {TANG_TISSUE_COL}")
    for ctx in CONTEXTS:
        label = CONTEXT_TISSUE[ctx]
        m = tang.obs[TANG_TISSUE_COL] == label
        if m.sum() == 0:
            raise RuntimeError(f"Tang has no cells with tissue={label}")
        out[ctx] = _annotate(tang[m, :], donor_col=DONOR_COLS[ctx])
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
            # Defensive: unreachable given keep_donors already enforces
            # >= MIN_CELLS_PER_GROUP in both arms and obs is not mutated here.
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
        f"[{ctx}] duplicate gene names in var_names — pseudobulk columns would be ambiguous"
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
