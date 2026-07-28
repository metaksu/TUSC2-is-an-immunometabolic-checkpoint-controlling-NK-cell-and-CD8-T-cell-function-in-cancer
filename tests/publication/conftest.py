import os
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import pytest

REPO = Path(__file__).resolve().parents[2]            # repository root
PKG = REPO / "tusc2_deg"

@pytest.fixture(scope="session")
def tmp_out(tmp_path_factory):
    return tmp_path_factory.mktemp("fig7_panels")

@pytest.fixture(scope="session")
def has_raw_tang():
    return (REPO.parent / "data" / "tang" / "comb_CD56_CD16_NK.h5ad").exists()
