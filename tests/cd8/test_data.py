"""data.py contracts: AnnData -> AnnData pseudobulk on CD8 huARdb."""
import numpy as np
import pandas as pd
import anndata
import pytest

from tusc2_deg.cd8 import data, config


def _make_toy_adata(n_cells_per_group: int = 12, n_genes: int = 20) -> anndata.AnnData:
    rng = np.random.default_rng(42)
    n_cells = 4 * n_cells_per_group
    X = rng.poisson(5, size=(n_cells, n_genes)).astype(np.int64)
    var = pd.DataFrame(index=[f"GENE{i}" for i in range(n_genes - 1)] + ["TUSC2"])
    obs = pd.DataFrame({
        "individual_id": (["d1"] * (2 * n_cells_per_group) + ["d2"] * (2 * n_cells_per_group)),
        "tusc2_status": (["pos"] * n_cells_per_group + ["neg"] * n_cells_per_group) * 2,
    })
    adata = anndata.AnnData(X=X, obs=obs, var=var)
    adata.raw = adata.copy()
    adata.uns["donor_col"] = "individual_id"
    return adata


def test_build_pseudobulk_returns_anndata():
    adata = _make_toy_adata()
    pb = data.build_pseudobulk(adata, "CD8_TIL")
    assert isinstance(pb, anndata.AnnData)


def test_build_pseudobulk_obs_schema():
    adata = _make_toy_adata()
    pb = data.build_pseudobulk(adata, "CD8_TIL")
    for col in ("donor", "tusc2_status", "n_cells"):
        assert col in pb.obs.columns
    assert pb.n_obs == 4


def test_build_pseudobulk_drops_underpowered_donors():
    adata = _make_toy_adata(n_cells_per_group=config.MIN_CELLS_PER_GROUP - 1)
    with pytest.raises(RuntimeError, match="MIN_CELLS_PER_GROUP"):
        data.build_pseudobulk(adata, "CD8_TIL")


def test_build_pseudobulk_counts_are_integer():
    adata = _make_toy_adata()
    pb = data.build_pseudobulk(adata, "CD8_TIL")
    X = pb.X.toarray() if hasattr(pb.X, "toarray") else pb.X
    assert np.issubdtype(X.dtype, np.integer)


def test_build_pseudobulk_preserves_tusc2():
    adata = _make_toy_adata()
    pb = data.build_pseudobulk(adata, "CD8_TIL")
    assert "TUSC2" in pb.var_names


def test_context_filters_use_huardb_columns():
    """The CD8_TIL context filter must reference the disease_type/meta_tissue_type
    columns that exist on huARdb v2 obs (solid-tumor TIL cells only)."""
    obs = pd.DataFrame({
        "disease_type":     ["Solid tumor", "Solid tumor", "Healthy"],
        "meta_tissue_type": ["TIL",         "PBMC",        "PBMC"],
    })
    assert data.CONTEXT_FILTERS["CD8_TIL"](obs).tolist() == [True, False, False]
    assert set(data.CONTEXT_FILTERS) == {"CD8_TIL"}
