"""pipeline.py contracts: stage caching + manifest."""
import json
import pandas as pd

from tusc2_deg.cd8 import pipeline, config


def test_cached_when_output_exists(tmp_path):
    expected = tmp_path / "x.csv"
    pd.DataFrame({"a": [1]}).to_csv(expected, index=False)

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


def test_force_from_propagates_downstream():
    """`--force-from degs` must rerun degs and gsea (downstream), not pseudobulk."""
    assert pipeline.STAGES == ["pseudobulk", "degs", "gsea"]
    assert pipeline._force_from("degs", "degs",       False) is True
    assert pipeline._force_from("degs", "gsea",       False) is True
    assert pipeline._force_from("degs", "pseudobulk", False) is False
    # `--force-from gsea` reruns only gsea
    assert pipeline._force_from("gsea", "gsea",       False) is True
    assert pipeline._force_from("gsea", "degs",       False) is False


def test_volcano_label_provenance_writer(tmp_path, monkeypatch):
    """The provenance writer emits one row per curated MUST_LABEL gene with
    its live LFC/padj/baseMean + a [sig]/[relaxed]/[apriori-disclose] tag."""
    monkeypatch.setattr(config, "AUDIT_DIR", tmp_path)
    ctx = config.MAIN_CONTEXTS[0]  # CD8_TIL
    curated = list(config.MUST_LABEL_BY_CONTEXT[ctx])
    # toy DEG frame: first curated gene fully-sig, second relaxed, third absent.
    deg = pd.DataFrame({
        "gene":           [curated[0], curated[1], "OTHER"],
        "log2FoldChange": [0.5,        0.03,       0.0],
        "padj":           [1e-5,       1e-3,       0.9],
        "baseMean":       [100.0,      100.0,      100.0],
    })
    out = pipeline._write_curated_label_provenance({ctx: deg}, [ctx])
    assert out.exists()
    prov = pd.read_csv(out, sep="\t")
    # one row per curated gene (present-or-absent disclosed)
    assert set(prov["gene"]) == set(curated)
    assert "tag" in prov.columns and "context" in prov.columns
    row0 = prov[prov["gene"] == curated[0]].iloc[0]
    assert row0["tag"] == "[sig]"


def test_manifest_writes(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "AUDIT_DIR", tmp_path)
    pipeline._write_manifest({"stage_a": "cached", "stage_b": "ran"})
    manifest = tmp_path / "pipeline_manifest.json"
    assert manifest.exists()
    data = json.loads(manifest.read_text(encoding="utf-8"))
    assert "stages" in data
    assert "config_snapshot" in data
    assert data["config_snapshot"]["CONTEXTS"] == config.CONTEXTS
