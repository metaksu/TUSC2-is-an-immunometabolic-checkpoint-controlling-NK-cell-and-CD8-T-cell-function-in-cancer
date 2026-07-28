"""pipeline.py contracts: stage caching + --force."""
import pandas as pd
import json
import pytest

from tusc2_deg.nk import pipeline, config


def test_stage_inputs_present(tmp_path, monkeypatch):
    """If stage outputs already exist on disk and --force not set, the
    stage runner returns 'cached' without rerunning."""
    monkeypatch.setattr(config, "DEG_DIR", tmp_path)
    expected = tmp_path / "deg_pseudobulk_Tang_Tumor.csv"
    pd.DataFrame({"gene": ["X"]}).to_csv(expected, index=False)

    result = pipeline._cached_or_run(
        out_path=expected,
        runner=lambda: (_ for _ in ()).throw(AssertionError("should not run")),
        force=False,
    )
    assert result == "cached"


def test_force_busts_cache(tmp_path):
    expected = tmp_path / "x.csv"
    pd.DataFrame({"a": [1]}).to_csv(expected, index=False)

    called = {"n": 0}
    def _runner():
        called["n"] += 1
        pd.DataFrame({"a": [2]}).to_csv(expected, index=False)

    pipeline._cached_or_run(out_path=expected, runner=_runner, force=True)
    assert called["n"] == 1


def test_curated_label_provenance_writer(tmp_path, monkeypatch):
    """The provenance writer emits one row per curated gene per context with
    live LFC/padj/baseMean and a [sig]/[relaxed]/[apriori-disclose] tag, so the
    curated-label set is auditable against the live DEG table."""
    monkeypatch.setattr(config, "AUDIT_DIR", tmp_path)
    ctx = "Tang_Tumor"
    curated = config.MUST_LABEL_BY_CONTEXT[ctx]
    # Build a DEG frame containing the curated genes plus extras; mix gate states.
    import numpy as np
    genes = list(curated) + ["EXTRA1", "EXTRA2"]
    deg = pd.DataFrame({
        "gene": genes,
        "log2FoldChange": [0.3] * len(genes),
        "padj": [0.001] * len(genes),
        "baseMean": [100.0] * len(genes),
    })
    out = pipeline._write_curated_label_provenance({ctx: deg}, [ctx])
    assert out.exists()
    prov = pd.read_csv(out, sep="\t")
    prov_ctx = prov[prov["context"] == ctx]
    # one row per curated gene for the context
    assert set(prov_ctx["gene"]) == set(curated)
    assert "tag" in prov_ctx.columns
    assert {"log2FoldChange", "padj", "baseMean"}.issubset(prov_ctx.columns)


def test_manifest_writes(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "AUDIT_DIR", tmp_path)
    pipeline._write_manifest({"stage_a": "cached", "stage_b": "ran"})
    manifest = tmp_path / "pipeline_manifest.json"
    assert manifest.exists()
    data = json.loads(manifest.read_text(encoding="utf-8"))
    assert "stages" in data
    assert "config_snapshot" in data
