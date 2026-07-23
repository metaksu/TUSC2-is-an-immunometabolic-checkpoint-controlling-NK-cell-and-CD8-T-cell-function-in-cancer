"""Input-data paths for the DEG pipelines.

All inputs resolve under a single data root, overridable with the
TUSC2_DATA_DIR environment variable (default: a ``data/`` directory beside the
package). The README documents the expected layout and the atlas accessions.
These paths are only needed for the optional from-raw pipeline; the offline
reproduce path renders every panel from the cached tables under
``tusc2_deg/{nk,cd8}/output`` and touches none of them.
"""
from __future__ import annotations
import os
from pathlib import Path

DATA_DIR = Path(os.environ.get("TUSC2_DATA_DIR", Path(__file__).resolve().parents[2] / "data"))

# NK atlas (Tang 2023).
TANG_H5AD = DATA_DIR / "tang" / "comb_CD56_CD16_NK.h5ad"

# CD8 atlas (huARdb v2, Xue 2025).
HUARDB_CD8_ALLGENES_H5AD = DATA_DIR / "huardb" / "huARdb_v2_GEX.CD8.all_genes.h5ad"
