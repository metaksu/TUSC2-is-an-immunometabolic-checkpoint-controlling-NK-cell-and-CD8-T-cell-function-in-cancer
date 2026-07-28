"""data.py contracts: AnnData -> AnnData pseudobulk."""
import numpy as np
import pandas as pd
import anndata
import pytest

from tusc2_deg.nk import data, config


def _make_toy_adata(n_cells_per_group: int = 12, n_genes: int = 20) -> anndata.AnnData:
    rng = np.random.default_rng(42)
    n_cells = 4 * n_cells_per_group  # 2 donors x 2 groups x cells_per_group
    X = rng.poisson(5, size=(n_cells, n_genes)).astype(np.int64)
    var = pd.DataFrame(index=[f"GENE{i}" for i in range(n_genes - 1)] + ["TUSC2"])
    obs = pd.DataFrame({
        "donor": (["d1"] * (2 * n_cells_per_group) + ["d2"] * (2 * n_cells_per_group)),
        "tusc2_status": (["pos"] * n_cells_per_group + ["neg"] * n_cells_per_group) * 2,
    })
    adata = anndata.AnnData(X=X, obs=obs, var=var)
    adata.raw = adata.copy()
    adata.uns["donor_col"] = "donor"
    return adata


def test_build_pseudobulk_returns_anndata():
    adata = _make_toy_adata()
    pb = data.build_pseudobulk(adata, "Tang_Tumor")
    assert isinstance(pb, anndata.AnnData)


def test_build_pseudobulk_obs_schema():
    adata = _make_toy_adata()
    pb = data.build_pseudobulk(adata, "Tang_Tumor")
    for col in ("donor", "tusc2_status", "n_cells"):
        assert col in pb.obs.columns, f"missing {col}"
    # 2 donors x 2 groups = 4 pseudosamples
    assert pb.n_obs == 4


def test_build_pseudobulk_drops_underpowered_donors():
    """Donor with <MIN_CELLS_PER_GROUP cells in either group must be dropped."""
    adata = _make_toy_adata(n_cells_per_group=config.MIN_CELLS_PER_GROUP - 1)
    with pytest.raises(RuntimeError, match="MIN_CELLS_PER_GROUP"):
        data.build_pseudobulk(adata, "Tang_Tumor")


def test_build_pseudobulk_counts_are_integer():
    adata = _make_toy_adata()
    pb = data.build_pseudobulk(adata, "Tang_Tumor")
    X = pb.X.toarray() if hasattr(pb.X, "toarray") else pb.X
    assert np.issubdtype(X.dtype, np.integer)


def test_build_pseudobulk_preserves_tusc2():
    adata = _make_toy_adata()
    pb = data.build_pseudobulk(adata, "Tang_Tumor")
    assert "TUSC2" in pb.var_names
