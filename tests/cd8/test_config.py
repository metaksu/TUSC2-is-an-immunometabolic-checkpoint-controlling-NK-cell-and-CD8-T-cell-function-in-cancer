"""config.py contracts: thresholds match NK invariants; CD8-specific extensions present."""
import pandas as pd
from tusc2_deg.cd8 import config


def test_contexts_locked():
    # Tumor-only reproducibility build: the CD8 context set is CD8_TIL only,
    # and MAIN_CONTEXTS must be a subset.
    assert set(config.CONTEXTS) == {"CD8_TIL"}
    assert set(config.MAIN_CONTEXTS).issubset(set(config.CONTEXTS))


def test_directions_locked():
    assert config.DIRECTIONS == ["UP", "DOWN"]


def test_lfc_calibrated_at_0_1():
    """LFC threshold is the single source of truth shared with the NK pipeline."""
    assert config.LFC_THRESH == 0.1


def test_padj_threshold():
    assert config.PADJ_THRESH == 0.05


def test_paths_resolve_under_project():
    assert config.ROOT.name == "cd8"
    assert config.OUT.name == "output"


def test_enrichr_libs_present():
    assert "GO_Biological_Process_2025" in config.ENRICHR_LIBS
    assert "Reactome_Pathways_2024" in config.ENRICHR_LIBS
    assert "MSigDB_Hallmark_2020" in config.ENRICHR_LIBS


def test_grouping_genes_contains_tusc2():
    """TUSC2 is the experimental axis (quasi-complete separation) and must be
    excluded from DEG tables / enrichment foreground / background."""
    assert "TUSC2" in config.GROUPING_GENES


def test_gene_sum_floor():
    """Pseudobulk gene-sum floor is exposed as a single named constant so the
    same value is used everywhere it gates low-count genes."""
    assert config.GENE_SUM_FLOOR == 10


def test_sig_lfc_ladder():
    assert config.SIG_LFC_LADDER[0] == config.LFC_THRESH
    assert config.SIG_LFC_LADDER == [0.1, 0.25, 0.585]


def test_context_n_donors_verified():
    """A verified per-context donor-count dict sourced from the README context
    table, so figures can self-disclose sample size. Donor n is displayed in the
    panel title/caption. Verified count: CD8_TIL=106."""
    assert hasattr(config, "CONTEXT_N_DONORS"), "CONTEXT_N_DONORS dict must exist"
    assert config.CONTEXT_N_DONORS["CD8_TIL"] == 106
    # every MAIN context must have a verified donor count for the disclosure
    for ctx in config.MAIN_CONTEXTS:
        assert ctx in config.CONTEXT_N_DONORS


def test_context_n_donors_modeled_cd8_til():
    """The paired-entry DESeq2 model donor n for CD8_TIL is 99 (of the 106
    context-filtered donors; 7 fail the >=10-cells-in-both-arms floor). This is
    the single source of truth for the "of 106" qualifier in the volcano title
    and composite footer. It must differ from CONTEXT_N_DONORS (the filter count)
    so the formatter renders the split."""
    assert hasattr(config, "CONTEXT_N_DONORS_MODELED"), "modeled-n dict must exist"
    assert config.CONTEXT_N_DONORS_MODELED["CD8_TIL"] == 99
    assert config.CONTEXT_N_DONORS_MODELED["CD8_TIL"] != config.CONTEXT_N_DONORS["CD8_TIL"]


def test_is_significant_three_part_gate():
    """Canonical significance = padj<PADJ AND baseMean>BASEMEAN AND |LFC|>LFC."""
    df = pd.DataFrame({
        "gene":           ["sig", "fail_lfc", "fail_padj", "fail_basemean"],
        "padj":           [0.01,  0.01,       0.20,        0.01],
        "log2FoldChange": [0.2,   0.05,       0.5,         0.2],
        "baseMean":       [100.0, 100.0,      100.0,       2.0],
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


def test_exhaustion_genes_present():
    """CD8-specific 5th family. Sources cited inline in config.py
    (Zheng 2021 Science; van der Leun 2020 Nat Rev Cancer; Wherry & Kurachi 2015)."""
    core = {"PDCD1", "HAVCR2", "LAG3", "TIGIT", "TOX", "ENTPD1", "CTLA4"}
    assert core.issubset(config.EXHAUSTION_GENES)


def test_gene_family_precedence_progenitor_restored():
    """The progenitor/memory family is its own category, and EOMES classifies
    into it rather than terminal exhaustion (matching NK, where EOMES is a
    maturation rheostat, not terminal). Precedence:
    progenitor > exhaustion > effector > activation > metabolism > prolif >
    stress > other. So the stem/Tpex genes (TCF7, SLAMF6, LEF1, IL7R, EOMES)
    resolve to 'progenitor'; terminal checkpoints (PDCD1) stay 'exhaustion'."""
    assert config.gene_family("PDCD1") == "exhaustion"
    assert config.gene_family("EOMES") == "progenitor"   # classified as progenitor, not exhaustion
    assert config.gene_family("GZMB")  == "effector"
    assert config.gene_family("IFNG")  == "effector"
    assert config.gene_family("CXCR4") == "effector"
    assert config.gene_family("TCF7")  == "progenitor"   # Tpex/stem master TF
    assert config.gene_family("SLAMF6") == "progenitor"  # Tpex surface marker
    assert config.gene_family("LEF1")  == "progenitor"   # stemness TF
    assert config.gene_family("IL7R")  == "progenitor"   # memory-precursor (precedence over activation)
    assert config.gene_family("ITK")   == "activation"
    assert config.gene_family("FOO")   == "other"


def test_eomes_not_in_exhaustion_genes():
    """EOMES must not be in EXHAUSTION_GENES (it is a memory/transitional
    rheostat, not a terminal-exhaustion TF), consistent with the NK taxonomy."""
    assert "EOMES" not in config.EXHAUSTION_GENES
    # the genuinely-terminal TFs stay
    assert "PRDM1" in config.EXHAUSTION_GENES
    assert "ID2" in config.EXHAUSTION_GENES


def test_progenitor_genes_used_by_gene_family():
    """PROGENITOR_GENES is the gene set that gene_family() uses to assign the
    progenitor family."""
    assert hasattr(config, "PROGENITOR_GENES")
    core = {"TCF7", "SLAMF6", "LEF1", "IL7R", "EOMES"}
    assert core.issubset(config.PROGENITOR_GENES)


def test_gene_family_palette_extended():
    """The colour palette must define an entry for every gene family."""
    for fam in ("progenitor", "exhaustion", "effector", "activation",
                "metabolism", "prolif", "stress", "other"):
        assert fam in config.GENE_FAMILY_COLORS, f"missing palette entry for {fam}"
    # progenitor colour must be distinct from the others
    prog = config.GENE_FAMILY_COLORS["progenitor"]
    others = [v for k, v in config.GENE_FAMILY_COLORS.items() if k != "progenitor"]
    assert prog not in others, "progenitor colour must be distinct"
