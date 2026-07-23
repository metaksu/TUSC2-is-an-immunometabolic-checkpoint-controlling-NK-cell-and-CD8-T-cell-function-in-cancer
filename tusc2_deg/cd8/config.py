"""Single source of truth for thresholds, paths, and palette.

Every module imports its thresholds from this file so that the volcano gate
and the enrichment input gate cannot drift apart. No hardcoded 0.1, 0.5,
0.05, or 5.0 should appear in any non-config module.
"""
from __future__ import annotations
import math
from pathlib import Path

# Tumor-only reproducibility build: the CD8 context is CD8_TIL on huARdb v2
# (Xue 2025 Nat Methods) — solid-tumor TIL (intratumoral, most disease-relevant;
# ~106 donors, 245k cells), paired with NK Tang_Tumor in the combined NK+CD8
# figure. It exceeds the >=15-donor cutoff (Squair 2021 NatComm) for donor-level
# pseudobulk DESeq2.
CONTEXTS = ["CD8_TIL"]
MAIN_CONTEXTS = ["CD8_TIL"]
DIRECTIONS = ["UP", "DOWN"]

# Per-context donor count from the huARdb v2 context filter
# (disease_type x meta_tissue_type): the donors passing the context filter. The
# pseudobulk DESeq2 design is ~donor + tusc2_status over the subset that also
# clears MIN_CELLS_PER_GROUP in both groups.
CONTEXT_N_DONORS = {
    "CD8_TIL": 106,
}

# Paired-entry DESeq2 MODEL donor n. CONTEXT_N_DONORS above is the context-FILTER
# count (disease_type x meta_tissue_type); this is the subset that additionally
# clears MIN_CELLS_PER_GROUP in BOTH TUSC2 arms and therefore actually enters the
# paired ~donor + tusc2_status pseudobulk model. For CD8_TIL, 7 of the 106 context
# donors fail the >=10-cells-in-both-arms floor, so the model n is 99. NK has no
# entry here (its filter and paired-entry counts are equal — no attrition).
CONTEXT_N_DONORS_MODELED = {
    "CD8_TIL": 99,
}

# DEG thresholds — single source of truth
LFC_THRESH          = 0.1   # |log2FC| gate on apeGLM-shrunken LFC (Zhu 2019).
                            # Discovery-permissive (~7% FC); not tuned to a target
                            # DEG count. Robustness shown via SIG_LFC_LADDER.
PADJ_THRESH         = 0.05
BASEMEAN_THRESH     = 5
MIN_CELLS_PER_GROUP = 10
GENE_SUM_FLOOR      = 10    # drop pseudobulk genes with summed counts below this

# Experimental-axis genes excluded from the DEG table, the enrichment
# foreground, AND the background. TUSC2 is the grouping variable and is
# quasi-complete-separated (TUSC2- cells are zero by construction); it is the
# only |LFC|>5 gene (built-in positive control) and would otherwise seed
# circular enrichment. Dropped in degs.run_deseq2/run_wilcoxon via _is_excluded
# so it never reaches a CSV. TUSC2 nonetheless remains in the DESeq2 counts input
# (and therefore the BH multiple-testing denominator); it is removed only from the
# output CSV via _is_excluded — the effect on padj is negligible (one gene of
# ~15k, factor ~(n+1)/n).
GROUPING_GENES = {"TUSC2"}

# Effect-size sensitivity ladder. DEG counts are reported at each |log2FC|
# cutoff (0.1 primary, 0.25 ~= 1.19x, 0.585 ~= 1.5x = sc-best-practices) so the
# biology is shown to survive a stricter gate rather than being calibrated to
# one count.
SIG_LFC_LADDER = [0.1, 0.25, 0.585]

# The calibrated DESeq2 lfcThreshold Wald test (Love 2014; null
# H0:|log2FC| <= LFC_NULL, alt 'greaterAbs') is NOT the significance definition.
# It is kept only to compute a per-gene calibrated padj as a reference/disclosure
# column (degs.run_deseq2 -> padj_calibrated). Significance is the
# |LFC|>LFC_THRESH 3-part gate (discovery floor + SIG_LFC_LADDER robustness +
# GSEA), which is shared with the NK pipeline for cross-pipeline comparability.
LFC_NULL        = math.log2(1.1)   # ~= 0.1375 — calibrated-test null (disclosure only)
ALT_HYPOTHESIS  = "greaterAbs"


def is_significant(df, direction=None):
    """Canonical significance gate (3-part): padj < PADJ_THRESH AND baseMean >
    BASEMEAN_THRESH AND |log2FoldChange| > LFC_THRESH. `direction='UP'/'DOWN'`
    restricts the fold-change sign. The calibrated lfcThreshold variant is not
    used here; it is retained only as the `padj_calibrated` disclosure column.
    Returns a boolean Series aligned to df.index.
    """
    base = (df["padj"] < PADJ_THRESH) & (df["baseMean"] > BASEMEAN_THRESH)
    lfc = df["log2FoldChange"]
    if direction == "UP":
        return base & (lfc > LFC_THRESH)
    if direction == "DOWN":
        return base & (lfc < -LFC_THRESH)
    return base & (lfc.abs() > LFC_THRESH)

# GSEA gene-set libraries.
# KEGG is excluded from the main panel: KEGG disease pathways (Parkinson Disease,
# Alzheimer Disease, Toll-Like Receptor Signaling Pathway, IL-17 Signaling
# Pathway) surface in immune contexts because they share OXPHOS / NF-kB / cytokine
# genes with NK biology, and the ALL-CAPS KEGG term casing visually outweighs
# Hallmark/Reactome at equal x-height. Hallmark + GO-BP + Reactome only.
ENRICHR_LIBS = [
    "GO_Biological_Process_2025",
    "Reactome_Pathways_2024",
    "MSigDB_Hallmark_2020",
]

# Visualisation palette
# Volcano dot colours — muted/desaturated off-red (UP) and off-blue (DOWN) so
# the dots recede and do NOT compete with the saturated gene-family LABEL colours
# (effector dark-red, prolif sky-blue, metab chocolate, etc.).
COLOR_UP   = "#d98880"   # dusty rose (off-red)
COLOR_DOWN = "#85a9c9"   # muted steel (off-blue)
COLOR_NS   = "#cfcfcf"   # light grey

# Volcano label gene families (highlight biology in the volcano panels).
#
# Sources (all citable, no heuristic prefix matching):
#   cytotox    — curated cytotoxicity-machinery gene set
#   cytokine   — curated cytokine/chemokine gene set, expanded with the
#                chemokine receptors that anchor NK trafficking biology
#                (CXCR3/4, CCR5/7) — these otherwise surface as METAB false
#                positives because Hallmark Glycolysis lists CXCR4 via HIF1A
#                regulation.
#   metabolism — union of HALLMARK_OXIDATIVE_PHOSPHORYLATION,
#                HALLMARK_GLYCOLYSIS, HALLMARK_FATTY_ACID_METABOLISM
#                (MSigDB v2025.1)
#   prolif     — union of HALLMARK_E2F_TARGETS, HALLMARK_G2M_CHECKPOINT,
#                HALLMARK_MYC_TARGETS_V1, HALLMARK_MITOTIC_SPINDLE
#                (MSigDB_Hallmark_2020)
#
# Precedence on overlap: exhaustion > cytotox > cytokine > metab > prolif > other.
# Rationale: CD8 dysfunction (exhaustion) is the central biology in tumor
# TILs and is the most informative readout for TUSC2+ vs TUSC2- contrasts.
# Cytotox/cytokine sit beneath because effector and exhausted CD8 share
# many granzyme/IFNG transcripts at intermediate states. METAB > PROLIF
# because TUSC2 → OXPHOS is the central finding; overlap genes (CYC1, COX5A,
# GOT2) are primarily mitochondrial — Myc Targets includes them only
# because Myc regulates mitochondrial biogenesis.
# Gene-set definitions are vendored with the package (MSigDB + curated).
GENE_SET_DIR = Path(__file__).resolve().parents[1] / "gene_sets"


def _load_grp(filename: str) -> set[str]:
    """Read MSigDB .grp file. Skip comment lines (#), URL footer, header."""
    path = GENE_SET_DIR / filename
    out: set[str] = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        s = line.strip()
        if not s or s.startswith("#") or s.startswith(">") or s.startswith("https"):
            continue
        if s.startswith("HALLMARK_"):
            continue
        out.add(s)
    return out


# Exhaustion (CD8-specific 5th family).
# Sources:
#   Zheng et al. 2021 *Science* 374:abe6474 — pan-cancer T atlas exhaustion module
#                (PDCD1, HAVCR2, LAG3, TIGIT, TOX, ENTPD1, CTLA4)
#   van der Leun et al. 2020 *Nat Rev Cancer* 20:218-232, Table 1 — TIL
#                dysfunction core markers
#   Wherry & Kurachi 2015 *Nat Rev Immunol* 15:486-499 — chronic-stimulation
#                transcription factors (TOX, BATF, NR4A1/2; PRDM1/ID2 terminal).
#                EOMES is NOT here — it is a memory/transitional rheostat
#                (PROGENITOR family; Intlekofer 2005; Zheng 2021).
#   TNFRSF9 (4-1BB) is included as a TIL-activation/exhaustion-coupled marker
#   (Bassez 2021 *Nat Med*; Caushi 2021 *Nature*).
# Stress / IEG / antioxidant / UPR family — cross-pipeline-conserved program.
# NR4A1/NR4A2/NR4A3 live in EXHAUSTION (Wherry & Kurachi 2015 — chronic-
# stimulation TFs; Chen 2019; Liu 2019), not here. TSC22D3 (GILZ; Yang 2022) is
# included. CD69/REL belong to activation. The antioxidant/redox subset stays.
STRESS_GENES = {
    "FOS", "FOSB", "JUN", "JUNB", "JUND", "EGR1", "EGR2", "EGR3",
    "IER2", "IER3", "ATF3", "BTG1", "BTG2",
    "HSPA1A", "HSPA1B", "HSPA6", "HSPA8", "HSPB1",
    "DNAJB1", "DNAJB4", "DNAJA1", "HSP90AA1",
    "TSC22D3",   # GILZ, glucocorticoid-induced AP-1/NF-kB inhibitor (Yang 2022)
    "TXNIP", "PRDX1", "PRDX2", "PRDX3", "GSTP1", "SOD1", "GPX1", "TXNRD1",
    "DDIT3", "GADD45A", "GADD45B", "DUSP1", "DUSP2", "PPP1R15A",
    "NFKBIA", "ZFP36", "ZFP36L1", "ZFP36L2", "TNFAIP3",
}

# EXHAUSTION (CD8-expanded): terminal-exhaustion TFs + TIL markers.
# NR4A3 is included (Chen 2019; Liu 2019). TNFRSF18 = GITR (Lowery 2024).
# EOMES, TCF7, and SLAMF6 are PROGENITOR, not exhaustion: EOMES is a
# memory/transitional rheostat (Intlekofer 2005; Zheng 2021) and TCF7 is the
# Tpex stem pole (Im 2016; Siddiqui 2019), so neither is a terminal-exhaustion
# TF. The genuinely-terminal TFs PRDM1/ID2 stay here (Wherry & Kurachi 2015;
# Paley 2012).
EXHAUSTION_GENES = {
    "PDCD1", "HAVCR2", "LAG3", "TIGIT", "CD96", "TOX", "TOX2", "ENTPD1", "CTLA4",
    "TNFRSF9", "TNFRSF18",
    "BATF", "NR4A1", "NR4A2", "NR4A3",
    "PRDM1", "ID2",        # terminal TFs (Wherry & Kurachi 2015; Paley 2012)
    "KLRG1", "ITGAE", "LAYN", "CXCL13", "CD200",
}

# EFFECTOR: cytotoxic machinery and effector cytokines on one co-axial effector
# axis (Zheng 2021 Science / Szabo 2019). Identical to NK EFFECTOR_GENES.
EFFECTOR_GENES = {
    # --- lytic-granule cargo / death ligands (cytotox core) ---
    "PRF1", "GZMA", "GZMB", "GZMH", "GZMK", "GZMM",      # granzymes/perforin (Turiello 2025)
    "NKG7",                                               # granule-effector (Ng 2020 Nat Immunol)
    "GNLY", "FGFBP2", "SERPINB9", "SRGN",                 # CD56dim cytotoxic sig (Crinier 2018)
    "FASLG", "TNFSF10",                                   # death-effector ligands
    "LGALS1", "EMP3",                                     # CD56dim/NK1 sig (Crinier 2018)
    # --- granule secretory / docking machinery ---
    "CTSC", "CTSW", "CST7", "LAMP1",                      # granule-resident proteases (Magister 2015)
    "SYTL1", "SYTL2", "SYTL3",                            # Slp granule delivery (Kurowska 2012 SYTL3)
    "RAB27A", "UNC13D", "STX11", "STXBP2",                # Rab27a/SNARE docking (Stinchcombe 2001)
    "ADGRG1",                                             # GPR56 pan-cytotoxic marker (Peng 2011)
    # --- secreted effector cytokines ---
    "IFNG", "TNF", "LTA", "LTB", "CSF2", "IL10", "IL32",  # T/NK cytokine module (Szabo 2019)
    # --- recruitment chemokines ---
    "XCL1", "XCL2",                                       # cDC1 recruitment (Bottcher 2018)
    "CCL3", "CCL4", "CCL4L2", "CCL3L1", "CCL5",           # effector chemokines (Tang 2023)
    "CXCL8", "CXCL10",                                    # inflammatory chemokines (Szabo 2019)
    # --- trafficking receptors (CXCR4 overrides Hallmark-Glycolysis tag) ---
    "CXCR3", "CXCR4", "CCR5", "CXCR6",
    # --- TNF-superfamily recruitment ligands ---
    "TNFSF8", "TNFSF14",
    # --- effector-lineage transcription factor ---
    "TBX21",                                             # T-bet, master effector TF (AACR CIR 2015)
}
# CCR7 is NOT here — it is a PROGENITOR (Tcm) marker.

# ACTIVATION (shared NK/CD8). Identical to NK ACTIVATION_GENES. CD69/REL/MCL1
# have their single home here. In CD8, CD2 resolves to activation (no nk_identity
# family); IL7R resolves to progenitor (precedence).
ACTIVATION_GENES = {
    "CD69",                          # earliest-activation IEG-surface (Cibrian 2017)
    "FYN", "LCK", "ITK", "TXK",      # Src/Tec proximal kinases (Salmond 2009)
    "ZAP70", "CD247", "PLCG1", "VAV1",  # TCR/CD3-zeta signalosome (Salmond 2009)
    "PIK3R1",                        # PI3K regulatory subunit
    "TESPA1",                        # TCR-signalling adaptor
    "REL", "NFKB1", "NFKB2", "RELB", # activation-induced NF-kB arm
    "TNFAIP3",                       # NF-kB feedback (activation-induced)
    "MCL1", "BCL2A1",                # activation-induced survival (Dunkle 2013)
    "IL7R",                          # homeostatic survival receptor (-> progenitor in CD8)
    "CD2", "CD27", "CD28", "CD44",   # costim / activation surface markers
    "KLRK1",                         # NKG2D, costimulatory TCR receptor (PubMed 15814668)
}

# PROGENITOR (stem/Tpex/Tcm/memory pole). Polar opposite of terminal exhaustion;
# up-in-TUSC2-NEG (positive LFC). Occupies the slot-3 legend position (NK's
# nk_identity analogue). Owns TCF7/SLAMF6/LEF1/IL7R + EOMES (the memory/
# transitional rheostat — Intlekofer 2005; Zheng 2021).
PROGENITOR_GENES = {
    "TCF7", "LEF1", "BACH2",        # stemness TFs (Im 2016 Nature; Yao 2021)
    "IL7R", "SLAMF6", "SATB1",      # progenitor survival/surface (Miller & Sen 2019 Nat Immunol)
    "EOMES",                        # memory/transitional rheostat (Intlekofer 2005; Zheng 2021) — NOT terminal
    "CCR7", "KLF2", "KLF3",         # Tcm/circulation (Zheng 2021)
    "NELL2", "FOXP1", "ZEB1",       # supporting memory TFs
}

# METABOLISM: union of Hallmark sets with mtDNA OXPHOS + lipid/glycolysis
# orphans. NOT direction-coherent (ACOT7/DBI negative, MT-/SLC2A3 positive) —
# read as "metabolic rewiring". Identical to NK METAB_EXTRA.
METAB_EXTRA = {
    "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6",
    "MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB", "MT-ATP6", "MT-ATP8",  # 13 mtDNA OXPHOS subunits
    "ACOT7",    # mito acyl-CoA thioesterase, lipid (Hunt 2002)
    "DBI",      # acyl-CoA-binding protein, lipid
    "SLC2A3",   # GLUT3 glucose transporter (Schafer 2019)
    "CMC1",     # Complex-IV/COX assembly factor
    "GAPDH",    # glycolysis enzyme absent from the Hallmark-Glycolysis .grp
}
METAB_GENES = (
    _load_grp("HALLMARK_OXIDATIVE_PHOSPHORYLATION.v2025.1.Hs.grp")
    | _load_grp("HALLMARK_GLYCOLYSIS.v2025.1.Hs.grp")
    | _load_grp("HALLMARK_FATTY_ACID_METABOLISM.v2025.1.Hs.grp")
    | METAB_EXTRA
)
# PROLIF: union of 4 Hallmark sets with the canonical Tirosh 2016 / Seurat
# cc.genes.updated.2019 S+G2M cycling signature. Capturing the cycling signature
# explicitly prevents ~3000 CD8 cell-cycle genes spilling into "other".
# DIRECTION: UP IN TUSC2-POSITIVE (CDT1/MKI67 negative neg-vs-pos log2FC; see output/degs), i.e. a
# proliferation surge in TUSC2+. Identical to NK PROLIF_EXTRA.
PROLIF_EXTRA = {
    # S-phase: replication licensing / fork / nucleotide supply
    "MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","MCM10","CDT1","CDC6","ORC1",
    "GINS2","CLSPN","PCNA","TYMS","RRM1","RRM2","DHFR","FEN1","UHRF1",
    "RAD51","RAD51AP1","DTL","GMNN","CHAF1B","POLE2","PCLAF","CDCA7",
    "FANCI","FANCD2","RBBP8","WDHD1","ATAD2","EXO1","HELLS","SLBP","CCNE2",
    # G2/M: kinetochore / centromere / condensin / spindle / mitotic kinases
    "MKI67","TOP2A","CCNB1","CCNB2","CCNA2","CDK1","PKMYT1","CDC25C","FOXM1",
    "BUB1","BUB1B","AURKA","AURKB","PLK1","BIRC5","UBE2C","CKS1B","CKS2",
    "PTTG1","STMN1","TUBA1B","TUBB","TUBB6","MYBL2","NUSAP1","TPX2","DLGAP5",
    "CENPA","CENPE","CENPF","CENPH","CENPN","CENPU","CENPW","CENPX","CENPO","CENPP",
    "ZWINT","NUF2","NDC80","SKA1","SKA3","HJURP","KIF11","KIF14","KIF23","KIFC1",
    "ASPM","NCAPG","NCAPG2","NCAPH","NCAPD2","SMC2","SMC4","SGO1","SGO2","ESCO2",
    "CDCA2","CDCA5","CDCA8","GTSE1","TTK","NEK2","ANLN","ECT2","HMMR","CEP55",
    # replication-dependent histones (S-phase-restricted transcription)
    "HIST1H1A","HIST1H1B","HIST1H2AJ","HIST1H3B","HIST1H3C","HIST1H3F",
    "HIST1H3G","HIST1H4C","HIST2H2AB",
}
PROLIF_GENES = (
    _load_grp("HALLMARK_E2F_TARGETS.MSigDB_Hallmark_2020.grp")
    | _load_grp("HALLMARK_G2M_CHECKPOINT.MSigDB_Hallmark_2020.grp")
    | _load_grp("HALLMARK_MYC_TARGETS_V1.MSigDB_Hallmark_2020.grp")
    | _load_grp("HALLMARK_MITOTIC_SPINDLE.MSigDB_Hallmark_2020.grp")
    | PROLIF_EXTRA
)

# Force-label genes. CD8 axes: checkpoint markers + TIL proliferation surge +
# cross-pipeline stress program. Forced labels do not inflate top_n; they replace
# the lowest-composite non-family pick. The selector silently drops any force gene
# that fails the gate in a panel, so a gene renders only where it is significant.
FORCE_LABEL_GENES = {
    # effector / cytotoxic / cytokine
    "GZMH", "GZMK", "GZMA",            # effector granzymes
    "GZMM", "IFNG", "TNF", "PRF1",     # complete cytotoxic story
    # proliferation surge (UP in TUSC2-POS)
    "MKI67", "TYMS", "TOP2A",
    # cross-pipeline conserved with NK
    "HSPA1A", "HSPA6", "CXCR4", "CCL5", "FOS", "JUN", "JUNB",
    # high-baseMean lymphotoxin
    "LTB",
    # progenitor / exhaustion axis
    "TCF7", "IL7R",                    # progenitor
    "SLAMF6", "NR4A1", "NR4A2", "NR4A3", "PRDM1", "VSIR",  # exhaustion/progenitor axis
    # activation / stress
    "CD69", "ITK",
    "TSC22D3",
    # terminal-exhaustion checkpoints
    "PDCD1", "TOX2", "TNFRSF9", "TNFRSF18", "KLRG1",
}

# Curated per-context volcano labels. Literature-tiered set spanning all four CD8
# axes; preserves the directional split (effector/progenitor/activation UP-in-NEG
# vs proliferation/OXPHOS/TOX2 UP-in-POS). Live per-gene significance status
# (full gate vs relaxed padj-sig |LFC|<0.1) is recorded in
# output/audit/curated_label_provenance.tsv, verified against
# output/degs/deg_pseudobulk_CD8_TIL.csv with is_significant().
# Sign: POSITIVE log2FC = higher in TUSC2-NEGATIVE cells.
# Honoured by the volcano label selector's MUST_LABEL/relaxed-gate branch;
# non-sig curated genes are silently dropped.
MUST_LABEL_BY_CONTEXT = {
    "CD8_TIL": {
        # effector / cytotoxic + cytokine (UP in TUSC2-NEG)
        "GZMK",   # pre-exhausted/effector-memory hallmark (Zheng 2021 Science)
        "PRF1",   # perforin lytic core
        "IFNG",   # effector cytokine (Grasso 2020 Cancer Cell)
        "TNF",    # effector cytokine
        "CCL5",   # effector chemokine
        "TBX21",  # T-bet, effector-lineage TF — non-IEG counter to the stress caveat
        # progenitor / stem-like (UP in TUSC2-NEG)
        "TCF7",   # stem/Tpex master TF (Im 2016 Nature; Siddiqui 2019 Immunity)
        "IL7R",   # memory-precursor survival (Joshi/Kaech 2007)
        "SLAMF6", # Tpex surface marker (Miller/Sen 2019 Nat Immunol)
        # exhaustion / checkpoint (mixed direction is the finding)
        "PDCD1",  # [UP in NEG] PD-1 (Zheng 2021; Wherry 2015)
        "CD96",   # [UP in NEG] co-inhibitory checkpoint, CD226 competitor
        "TOX2",   # [UP in POS] exhaustion TF (Khan 2019 Nature)
        "KLRG1",  # terminal-differentiation marker (Joshi/Kaech 2007)
        # activation / costimulation (UP in TUSC2-NEG)
        "KLRK1",  # NKG2D costimulatory receptor — non-IEG
        "CD69",   # activation surface
        "FOS",    # AP-1
        # proliferation surge (UP in TUSC2-POS; dominant axis) — 3 reps so the DOWN
        # side is not flooded by cell-cycle at the expense of the programs below.
        "MKI67",  # universal proliferation marker (Tirosh 2016)
        "TOP2A",  # G2-M
        "MYBL2",  # cell-cycle TF
        # other downregulated programs (UP in TUSC2-POS) — varied beyond prolif:
        "COX8A",  # OXPHOS Complex IV (metabolism)
        "LDHA",   # glycolysis terminal enzyme (metabolism, high baseMean)
        "GSTP1",  # glutathione-S-transferase — antioxidant axis (stress)
        "PRDX1",  # peroxiredoxin-1 — antioxidant axis (stress)
        "LGALS1", # galectin-1, immunosuppressive effector (high baseMean)
    },
}

GENE_FAMILY_COLORS = {
    "progenitor":  "#1F78B4",   # blue (distinct from prolif sky-blue) — stem/Tpex/memory: TCF7/SLAMF6/LEF1/IL7R/EOMES
    "exhaustion":  "#7B1FA2",   # deep magenta — PDCD1/HAVCR2/LAG3/TIGIT/TOX + terminal TFs (PRDM1/ID2)
    "effector":    "#8B0000",   # dark red — granzymes/perforin + effector cytokines (merged)
    "activation":  "#117733",   # green (Wong/Okabe-Ito) — proximal signalosome + CD69/REL/MCL1 (shared)
    "metabolism":  "#D2691E",   # chocolate — OXPHOS / glycolysis / FA oxidation / mtDNA
    "prolif":      "#56B4E9",   # sky blue (Okabe-Ito) — E2F/G2M/MYC/spindle + Tirosh cycling
    "stress":      "#0B7A5F",   # deep teal — AP-1 / HSP / antioxidant / TSC22D3
    "other":       "#2c3e50",   # dark slate (default)
    "tusc2":       "#000000",   # black bold (reserved)
}


def gene_family(gene: str) -> str:
    """Return one of {'progenitor', 'exhaustion', 'effector', 'activation',
    'metabolism', 'prolif', 'stress', 'other'}.

    EOMES belongs to 'progenitor', not terminal exhaustion (memory/transitional
    rheostat — matching NK). The stem/Tpex genes (TCF7/SLAMF6/LEF1/IL7R/EOMES)
    resolve to 'progenitor'.

    Precedence: progenitor > exhaustion > effector > activation > metabolism >
    prolif > stress > other. progenitor FIRST so the stem/memory pole wins over
    the dual-membership genes (IL7R also in ACTIVATION; the stem TFs are the polar
    opposite of terminal exhaustion and must not be coloured as exhausted).
    exhaustion next (CD8 dysfunction axis, owns terminal checkpoints/TFs).
    effector, then activation. metabolism > prolif (TUSC2->OXPHOS thesis).
    stress last (artifact-confounded).
    """
    g = str(gene).upper()
    if g in PROGENITOR_GENES:
        return "progenitor"
    if g in EXHAUSTION_GENES:
        return "exhaustion"
    if g in EFFECTOR_GENES:
        return "effector"
    if g in ACTIVATION_GENES:
        return "activation"
    if g in METAB_GENES:
        return "metabolism"
    if g in PROLIF_GENES:
        return "prolif"
    if g in STRESS_GENES:
        return "stress"
    return "other"

# Paths (relative to this file)
ROOT         = Path(__file__).resolve().parent
OUT          = ROOT / "output"
DEG_DIR      = OUT / "degs"
GSEA_DIR     = OUT / "gsea"
AUDIT_DIR    = OUT / "audit"

# Gene-name exclusions for DEG output. Only ribosomal RPL/RPS prefixes are
# dropped; mtDNA-encoded "MT-" genes are kept. Best-practice (decoupler-py,
# muscat, sc-best-practices Ch.18, Love DESeq2 vignette) keeps mtDNA-encoded
# genes in the matrix — they are the most direct readout of mitochondrial
# respiratory chain function. Mirrors the NK config so both pipelines are
# consistent.
EXCLUDE_PREFIXES = ("RPL", "RPS")
