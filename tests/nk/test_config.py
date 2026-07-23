"""config.py contracts: single source of truth for thresholds and paths."""
import pandas as pd
from tusc2_deg.nk import config


def test_contexts_locked():
    # Tumor-only reproducibility build: the NK context set is Tang_Tumor only.
    assert config.CONTEXTS == ["Tang_Tumor"]
    assert config.MAIN_CONTEXTS == ["Tang_Tumor"]


def test_directions_locked():
    assert config.DIRECTIONS == ["UP", "DOWN"]


def test_lfc_calibrated_at_0_1():
    """The LFC threshold is a single shared constant so the volcano plot and
    the enrichment input use the same cutoff — must equal 0.1."""
    assert config.LFC_THRESH == 0.1


def test_padj_threshold():
    assert config.PADJ_THRESH == 0.05


def test_paths_resolve_under_project():
    assert config.ROOT.name == "nk"
    assert config.OUT.name == "output"


def test_enrichr_libs_present():
    assert "GO_Biological_Process_2025" in config.ENRICHR_LIBS or \
           "GO_Biological_Process_2023" in config.ENRICHR_LIBS
    assert "Reactome_Pathways_2024" in config.ENRICHR_LIBS
    assert "MSigDB_Hallmark_2020" in config.ENRICHR_LIBS


def test_grouping_genes_contains_tusc2():
    """TUSC2 is the experimental axis (quasi-complete separation) and must be
    excluded from DEG tables / enrichment foreground / background."""
    assert "TUSC2" in config.GROUPING_GENES


def test_gene_sum_floor():
    """The pseudobulk gene-sum floor is a single named constant so every module
    filters low-count genes at the same cutoff (>=10)."""
    assert config.GENE_SUM_FLOOR == 10


def test_sig_lfc_ladder():
    """Effect-size sensitivity ladder; primary gate first, then stricter cuts."""
    assert config.SIG_LFC_LADDER[0] == config.LFC_THRESH
    assert config.SIG_LFC_LADDER == [0.1, 0.25, 0.585]


def test_is_significant_three_part_gate():
    """Significance is a 3-part gate: padj<PADJ AND baseMean>BASEMEAN AND
    |log2FC|>LFC_THRESH. A small-|LFC| gene (|LFC|<LFC_THRESH) is NOT significant
    even when padj and baseMean pass."""
    df = pd.DataFrame({
        "gene":           ["sig", "small_lfc_fail", "fail_padj", "fail_basemean"],
        "padj":           [0.01,  0.01,             0.20,        0.01],
        "log2FoldChange": [0.2,   0.05,             0.5,         0.2],
        "baseMean":       [100.0, 100.0,            100.0,       2.0],
    })
    assert list(config.is_significant(df)) == [True, False, False, False]


def test_is_significant_direction():
    df = pd.DataFrame({
        "padj":           [0.01,  0.01],
        "log2FoldChange": [0.2,  -0.2],
        "baseMean":       [100.0, 100.0],
    })
    assert list(config.is_significant(df, direction="UP"))   == [True, False]
    assert list(config.is_significant(df, direction="DOWN")) == [False, True]


def test_eomes_unified_out_of_exhaustion():
    """Cross-pipeline parity guard: EOMES is a memory/transitional maturation
    rheostat, NOT terminal exhaustion, in both pipelines (CD8 -> progenitor;
    NK -> stress/IEG via STRESS_IEG.grp). Mirrors the CD8
    test_eomes_not_in_exhaustion_genes so the NK side cannot silently classify
    EOMES into EXHAUSTION_GENES. Collins-Zhang 2021 (NK);
    Intlekofer 2005 / Zheng 2021 (CD8)."""
    assert "EOMES" not in config.EXHAUSTION_GENES
    assert config.gene_family("EOMES") != "exhaustion"
    # the genuinely-dysfunction checkpoints stay in exhaustion
    assert "PDCD1" in config.EXHAUSTION_GENES
