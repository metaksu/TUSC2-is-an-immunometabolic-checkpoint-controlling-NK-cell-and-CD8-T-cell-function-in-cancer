"""Single source of truth for thresholds, paths, and palette.

Every module imports its thresholds from this file so that the volcano and
enrichment stages share one set of cutoffs. Hardcoded numeric thresholds
(0.1, 0.5, 0.05, 5.0) should not appear in non-config modules.
"""
from __future__ import annotations
import math
from pathlib import Path

# Tumor-only reproducibility build: the NK main context is Tang_Tumor,
# paired with CD8_TIL in the combined NK+CD8 figure.
CONTEXTS = ["Tang_Tumor"]
MAIN_CONTEXTS = ["Tang_Tumor"]
DIRECTIONS = ["UP", "DOWN"]

# Verified per-context donor count: donors that clear the MIN_CELLS_PER_GROUP
# gate in BOTH tusc2_status arms and therefore enter the paired pseudobulk DESeq2
# model. Held as a constant so it need not be threaded from the live pseudobulk;
# re-sync if the donor gate or input data change.
CONTEXT_N_DONORS = {
    "Tang_Tumor":     110,
}

# DEG thresholds — single source of truth
LFC_THRESH          = 0.1   # |log2FC| gate on apeGLM-shrunken LFC (Zhu 2019).
                            # Discovery-permissive (~7% FC), not tuned to a target
                            # count. Robustness is shown via SIG_LFC_LADDER below
                            # rather than asserted at a single cutoff.
PADJ_THRESH         = 0.05
BASEMEAN_THRESH     = 5
MIN_CELLS_PER_GROUP = 10
GENE_SUM_FLOOR      = 10    # drop pseudobulk genes with summed counts below this

# Experimental-axis genes excluded from the DEG table, the enrichment
# foreground, AND the background. TUSC2 is the grouping variable and is
# quasi-complete-separated (TUSC2- cells are zero by construction), so apeGLM
# inflates its |LFC| (it is the single strongest DOWN "hit") and it would seed
# its own enrichment terms — circular. Dropped in degs.run_deseq2/run_wilcoxon
# via _is_excluded so it never reaches a CSV.
# BH-denominator note: TUSC2 nonetheless remains in the DESeq2 counts input (and
# therefore in the Wald tests and the BH multiple-testing denominator); it is
# removed only from the output CSV via _is_excluded. The effect on padj is
# negligible (one extra gene among ~15k).
GROUPING_GENES = {"TUSC2"}

# Effect-size sensitivity ladder. DEG counts are reported at each |log2FC|
# cutoff (0.1 primary, 0.25 ~= 1.19x, 0.585 ~= 1.5x = sc-best-practices), so the
# biology is shown to survive a stricter gate rather than being calibrated to one
# count.
SIG_LFC_LADDER = [0.1, 0.25, 0.585]

# The calibrated DESeq2 lfcThreshold Wald test (Love 2014; null H0:|log2FC| <=
# LFC_NULL, alt 'greaterAbs') is not the significance definition: near-threshold,
# composition-mediated effects need not exceed a 1.1x null given their SE.
# Significance stays the |LFC|>LFC_THRESH 3-part gate (discovery floor +
# SIG_LFC_LADDER robustness + GSEA). LFC_NULL/ALT_HYPOTHESIS exist only to compute
# the calibrated padj as a per-gene disclosure/reference column
# (degs.run_deseq2 -> padj_calibrated) so Methods can report how many genes exceed
# the calibrated lfcThreshold null. This is a reference computation, not the
# significance definition.
LFC_NULL        = math.log2(1.1)   # ~= 0.1375 — calibrated-test null (disclosure only)
ALT_HYPOTHESIS  = "greaterAbs"


def is_significant(df, direction=None):
    """Canonical significance gate (3-part): padj < PADJ_THRESH AND baseMean >
    BASEMEAN_THRESH AND |log2FoldChange| > LFC_THRESH. `direction='UP'/'DOWN'`
    restricts the fold-change sign. The calibrated lfcThreshold variant is not
    used as the gate; its result is kept only as the `padj_calibrated` disclosure
    column. Returns a boolean Series aligned to df.index.
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
# the dots recede and do not compete with the saturated gene-family LABEL colours
# (effector dark-red, prolif sky-blue, metab chocolate, etc.).
COLOR_UP   = "#d98880"   # dusty rose (off-red)
COLOR_DOWN = "#85a9c9"   # muted steel (off-blue)
COLOR_NS   = "#cfcfcf"   # light grey

# Volcano label gene families (highlight biology in the volcano panels).
#
# Sources (all citable, no heuristic prefix matching):
#   effector   — curated cytotoxicity-machinery and cytokine/chemokine gene
#                sets, expanded with the chemokine receptors that anchor NK
#                trafficking biology (CXCR3/4, CCR5/7) — these surface
#                otherwise as METAB false positives because Hallmark
#                Glycolysis lists CXCR4 via HIF1A regulation.
#   metabolism — union of HALLMARK_OXIDATIVE_PHOSPHORYLATION,
#                HALLMARK_GLYCOLYSIS, HALLMARK_FATTY_ACID_METABOLISM
#                (MSigDB v2025.1, vendored under gene_sets/)
#   prolif     — union of HALLMARK_E2F_TARGETS, HALLMARK_G2M_CHECKPOINT,
#                HALLMARK_MYC_TARGETS_V1, HALLMARK_MITOTIC_SPINDLE
#                (MSigDB_Hallmark_2020, vendored under gene_sets/)
#
# Family precedence on overlap is defined in gene_family() below
# (exhaustion > effector > activation > metabolism > prolif > stress > other);
# metabolism wins over prolif because TUSC2 → OXPHOS is the central finding and
# the shared genes (CYC1, COX5A, GOT2) are primarily mitochondrial (Myc Targets
# includes them only because Myc regulates mitochondrial biogenesis).
# Gene-set definitions are vendored with the package (MSigDB + curated).
GENE_SET_DIR = Path(__file__).resolve().parents[1] / "gene_sets"


def _load_grp(filename: str) -> set[str]:
    """Read MSigDB-format .grp file. Skip comment lines (#), URL footer, header."""
    path = GENE_SET_DIR / filename
    out: set[str] = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        s = line.strip()
        if not s or s.startswith("#") or s.startswith(">") or s.startswith("https"):
            continue
        if s.startswith("HALLMARK_") or s.startswith("NK_TUSC2_"):
            continue
        out.add(s)
    return out


# Stress / IEG / antioxidant / UPR family — loaded from disk; provenance lives in
# the .grp file header. Loaded the same way as metab/prolif so methods can cite a
# single set of .grp files rather than a mix of literals and disk reads.
# STRESS_EXTRA: TSC22D3 (GILZ, glucocorticoid-induced AP-1/NF-kB inhibitor; Yang
# 2022) is a real orphan in both pipelines. CD69/REL/MCL1 do not go here (->
# activation).
STRESS_EXTRA = {"TSC22D3"}
STRESS_GENES = _load_grp("NK_TUSC2_STRESS_IEG.curated.grp") | STRESS_EXTRA

# NK surface-receptor family — loaded from disk. Captures the KLR / NCR / NCAM1 /
# SELL / FCGR3A / KIR axis + inhibitory checkpoints (PDCD1, LAG3, HAVCR2, TIGIT).
NK_RECEPTOR_GENES = _load_grp("NK_TUSC2_NK_RECEPTORS.curated.grp")

# Exhaustion / dysfunction axis — shared family taxonomy with the CD8 volcanos.
# PDCD1/HAVCR2/LAG3/TIGIT/TOX/CTLA4 also live in NK_RECEPTOR_GENES (the .grp file
# lists inhibitory checkpoints); exhaustion takes precedence in gene_family().
# NR4A1/NR4A2 sit in exhaustion (Wherry & Kurachi 2015 — chronic-stimulation TFs);
# NR4A3 stays in stress.
# EOMES is intentionally excluded from NK exhaustion: it is the NK maturation
# rheostat (Collins/Zhang 2021), not a dysfunction marker, and sits in
# STRESS_IEG.grp. EOMES is a memory/transitional rheostat in both pipelines
# (never exhaustion; Intlekofer 2005 / Zheng 2021).
EXHAUSTION_GENES = {
    "PDCD1", "HAVCR2", "LAG3", "TIGIT", "TOX", "ENTPD1", "CTLA4",
    "TNFRSF9", "BATF", "NR4A1", "NR4A2",
}
# EFFECTOR: lytic-effector and effector-cytokine output share one colour. Zheng
# 2021 Science / Szabo 2019 Nat Commun treat them as one co-axial functional axis.
# Per-context label curation (MUST_LABEL_BY_CONTEXT) can still pick granzyme vs
# IFNG reps. Members: lytic-granule cargo + death ligands + secretory/docking
# machinery + secreted effector cytokines/chemokines + trafficking receptors.
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
    "XCL1", "XCL2",                                       # NK-cDC1 recruitment (Bottcher 2018)
    "CCL3", "CCL4", "CCL4L2", "CCL3L1", "CCL5",           # IFNG-hi NK co-express (Tang 2023)
    "CXCL8", "CXCL10",                                    # inflammatory chemokines (Szabo 2019)
    # --- trafficking receptors (CXCR4 overrides Hallmark-Glycolysis tag) ---
    "CXCR3", "CXCR4", "CCR5", "CXCR6",
    # --- TNF-superfamily recruitment ligands (CD8-relevant, kept for parity) ---
    "TNFSF8", "TNFSF14",
}
# CCR7 is NOT here — it belongs to nk_identity (Tcm/circulation marker).

# ACTIVATION (shared NK/CD8). Proximal TCR/NK-receptor signalosome
# kinases/adaptors + early-activation surface markers + activation-induced
# NF-kB/survival arm. Home for ITK/FYN/CD69/REL/MCL1. CD69/REL/MCL1 live here,
# not in stress.
ACTIVATION_GENES = {
    "CD69",                          # earliest-activation IEG-surface (Cibrian 2017)
    "FYN", "LCK", "ITK", "TXK",      # Src/Tec proximal kinases (Salmond 2009)
    "ZAP70", "CD247", "PLCG1", "VAV1",  # TCR/CD3-zeta signalosome (Salmond 2009)
    "PIK3R1",                        # PI3K regulatory subunit
    "TESPA1",                        # TCR-signalling adaptor
    "REL", "NFKB1", "NFKB2", "RELB", # activation-induced NF-kB arm
    "TNFAIP3",                       # NF-kB feedback (activation-induced)
    "MCL1", "BCL2A1",                # activation-induced survival (Dunkle 2013)
    "IL7R",                          # homeostatic survival receptor
    "CD2", "CD27", "CD28", "CD44",   # costim / activation surface markers
}
# nk_identity. Base = KLR/NCR/NCAM1/SELL/FCGR3A/KIR + inhibitory checkpoints from
# the .grp, plus documented orphans inline:
NK_IDENTITY_EXTRA = {
    "CLEC2D",   # LLT1, cognate ligand of KLRB1/CD161 (Germain 2011; Braud 2018)
    "CD2",      # CD2-CD58 NK synapse / adaptive-NK receptor (Rolle 2016)
    "CCR7",     # Tcm/circulation marker
}
# LAIR2 NOT added (soluble decoy, negative direction, 1 ctx). Stays 'other'.
NK_IDENTITY_GENES = NK_RECEPTOR_GENES | NK_IDENTITY_EXTRA

# METABOLISM: union of Hallmark sets with documented orphans. Lane is not
# direction-coherent (ACOT7/DBI negative, MT-/SLC2A3 positive) — frame as
# "metabolic rewiring", show per-gene sign.
METAB_EXTRA = {
    "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6",
    "MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB", "MT-ATP6", "MT-ATP8",  # 13 mtDNA OXPHOS subunits
    "ACOT7",    # mito acyl-CoA thioesterase, lipid (Hunt 2002)
    "DBI",      # acyl-CoA-binding protein, lipid (admitted for direction-rule consistency)
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
# PROLIF: union of the 4 Hallmark sets with the canonical Tirosh 2016 / Seurat
# cc.genes.updated.2019 S+G2M cycling signature. DIRECTION: this module is UP IN
# TUSC2-POSITIVE (CDT1/MKI67 negative in neg-vs-pos; values in output/degs/) — caption must say
# "proliferation surge in TUSC2+". (NK has ~0 cycling orphans; applied for config
# parity with CD8, harmless.)
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

# Force-label genes.
# Guaranteed volcano labels when the gene is a significant DEG in that panel.
# Picked from the top-effect-size and sanity-check lists because the standard
# composite (|LFC| × log10(baseMean+1)) ranking undersold them. Forced labels do
# not inflate top_n; they replace the lowest-composite non-family pick.
# No per-gene LFC/padj/baseMean values are hard-coded here — they drift every
# re-run; live values live in output/degs/. Cross-pipeline anchors are listed for
# parity with the CD8 panel; the selector silently drops any gene that fails the
# gate in a panel, so a gene only renders where it is genuinely significant.
FORCE_LABEL_GENES = {
    "CXCR4",     # NK trafficking receptor
    "HSPA6",     # cross-pipeline stress
    "CCL4L2",    # NK effector chemokine
    "KLRF1",     # NKp80 activating receptor
    "KLRB1",     # CD161 NK identity
    "KLRC1",     # NKG2A inhibitory receptor
    "NCAM1",     # CD56 NK marker
    "SELL",      # CD62L homing marker
    "GZMH",      # cross-pipeline effector
    "GZMK",      # effector
    "NKG7",      # pan-cytotoxicity anchor (Turiello 2025; Ng 2020)
    "FOS", "JUN", "JUNB",  # AP-1 axis anchors
    "DNAJB1", "HSPA1A",    # HSP axis anchors
    # in coloured families:
    "CD69", "MCL1", "IL7R",  # activation
    "CLEC2D",                # nk_identity
    "CMC1",                  # metabolism
    "SYTL3",                 # effector (granule delivery)
    "TSC22D3",               # stress
}

# Curated per-context volcano labels.
#
# For MAIN_CONTEXTS this overrides the family-quota selector so the volcanos
# surface biologically representative DEGs from each family instead of the
# highest-composite housekeeping/spliceosome hits that crowd the top_n picker.
# Selector contract: every gene listed here that is significant (padj <
# PADJ_THRESH, |LFC| > LFC_THRESH, baseMean > BASEMEAN_THRESH) is labeled;
# non-significant entries are silently dropped (no off-axis labels).
#
# Per-context rationale, with sig DEGs from output/degs/ at the pipeline-config
# thresholds (LFC>0.1, padj<0.05, baseMean>5):
#
# Tang_Tumor (20 labels): 3 cytokines + 2 NK receptors + 4 stress/IEG (AP-1 +
# NFkB + nuclear-receptor) + 4 metabolism (glycolysis + OXPHOS reps) + 3 prolif
# (E2F/replication forks, not spliceosome) + 4 "other" (NK biology:
# IL7R/CD69/BCL2A1/CCL4L2). Favors canonical COX5A/ALDOA over generic stats
# winners (SNRPB, SNRPD3, RANBP1, EIF5A, PPIA, TALDO1) and low-baseMean NDUFA6.
# Literature-tiered curated label set. Each gene's live
# significance status (full-gate hit / relaxed sub-threshold |LFC| / NS a-priori
# canonical anchor) is recorded per-run in
# output/audit/curated_label_provenance.tsv rather than as a hard-coded tag here;
# the tags would drift every re-run. The selector renders each gene as a full
# hit, via the MUST_LABEL relaxed-gate branch, or (for NS canonical anchors) only
# where it becomes significant.
# Sign: POSITIVE log2FC = higher in TUSC2-NEGATIVE cells.
MUST_LABEL_BY_CONTEXT = {
    "Tang_Tumor": {          # MAIN; receptor/AP-1/chemokine UP in TUSC2-NEG;
                             # glycolysis/OXPHOS reps DOWN in TUSC2-NEG (= up in TUSC2+).
        # effector / cytotoxic
        "CCL3",   # NK-cDC1 recruitment chemokine (Bottcher 2018 Cell)
        "CCL4",   # NK effector chemokine (Bottcher 2018 Cell)
        "SRGN",   # serglycin granule proteoglycan (Crinier 2018 Immunity)
        "NKG7",   # pan-cytotoxicity anchor, more stable than GZMB/PRF1 (Turiello 2025 EJI; Ng 2020 Nat Immunol)
        # NK receptors ('other' colour). KLRF1/KLRB1 stay; KLRD1 + TYROBP are
        # excluded as sub-threshold relaxed anchors with no family colour.
        "KLRF1",  # NKp80 activating receptor
        "KLRB1",  # CD161 NK identity
        # cytokine / trafficking
        "CXCR4",  # NK trafficking receptor
        # activation / stress (artifact-robust AP-1/NF-kB arm)
        "FOS",    # AP-1 IEG (van den Brink 2017 — disclose dissociation caveat)
        "JUNB",   # AP-1 IEG
        "NFKBIA", # NF-kB feedback (TUSC2 suppresses Stressed program)
        "DUSP2",  # strongest AP-1 MAPK phosphatase
        "CD69",   # earliest-activation surface (Cibrian 2017)
        # metabolism (glycolysis + OXPHOS reps; central finding, DOWN in TUSC2-neg)
        "ENO1",   # glycolysis (O'Brien & Finlay 2019 NRI; Assmann 2017 Nat Immunol)
        "PKM",    # glycolysis
        "ALDOA",  # glycolysis
        "COX5A",  # OXPHOS subunit rep (disclose composition-mediated)
        "NDUFA6", # OXPHOS Complex I subunit (DOWN = up in TUSC2+)
        "COX6A1", # OXPHOS Complex IV subunit
        "GAPDH",  # glycolysis, canonical (high baseMean)
        # proliferation / DNA replication (DOWN in TUSC2-neg = up in TUSC2+)
        "MCM7",   # replication licensing, E2F target
        "MCM5",   # replication licensing
        "TFDP1",  # E2F dimerization partner
        # other NK biology (genuinely sig, not stats winners)
        "IL7R",   # NK homeostasis (Lopez-Verges 2010 Blood)
    },
}

# ONE Okabe-Ito family palette, identical dict literal in both pipeline configs
# so the figure legend + any shared gene is identical across the NK and CD8
# panels. NK genes never resolve to 'progenitor'/'tusc2' (those families are
# CD8-only / reserved), but the keys are present for parity so the shared
# volcano legend logic is uniform. Shared families are hex-identical to CD8;
# progenitor is a distinct blue (#1F78B4) from prolif sky-blue (#56B4E9).
GENE_FAMILY_COLORS = {
    "progenitor":   "#1F78B4",   # blue (distinct from prolif sky-blue) — CD8 stem/Tpex/memory; unused by NK genes
    "exhaustion":   "#7B1FA2",   # deep magenta — PDCD1/HAVCR2/LAG3/TIGIT/TOX (shared with CD8)
    "effector":     "#8B0000",   # dark red — granzymes/perforin + effector cytokines (merged)
    "activation":   "#117733",   # green (Wong/Okabe-Ito) — proximal signalosome + CD69/REL/MCL1 (shared)
    "metabolism":   "#D2691E",   # chocolate — OXPHOS / glycolysis / FA oxidation / mtDNA
    "prolif":       "#56B4E9",   # sky blue (Okabe-Ito) — E2F/G2M/MYC/spindle + Tirosh cycling
    "stress":       "#0B7A5F",   # deep teal — AP-1 / HSP / NR4A3 / antioxidant / TSC22D3
    "other":        "#2c3e50",   # dark slate (default)
    "tusc2":        "#000000",   # black bold (reserved; parity with CD8)
}


def gene_family(gene: str) -> str:
    """Return one of {'exhaustion', 'effector', 'activation', 'metabolism',
    'prolif', 'stress', 'other'}. KLR/NCR/NCAM1/SELL/CLEC2D genes resolve to
    'other'; the dual-membership CD2 resolves to 'activation'.

    Precedence: exhaustion > effector > activation > metabolism > prolif >
    stress > other. exhaustion first (dysfunction signal, owns shared
    checkpoints). effector second. activation third. metabolism > prolif (TUSC2 ->
    OXPHOS finding). stress last (artifact-confounded; only colours when not
    claimed elsewhere).
    """
    g = str(gene).upper()
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
AUDIT_DIR    = OUT / "audit"
GSEA_DIR     = OUT / "gsea"

# Gene-name exclusions for DEG output.
# "MT-" is intentionally NOT excluded: mtDNA-encoded protein-coding OXPHOS genes
# (MT-ND1-6, MT-CO1-3, MT-ATP6/8, MT-CYB) are the most direct readout of
# mitochondrial respiratory-chain function and should reach the volcano.
# Best-practice (decoupler-py, muscat, sc-best-practices Ch.18, Love DESeq2
# vignette) keeps MT- in the matrix and masks it only at visualization if needed.
# RPL/RPS stay excluded (translational housekeeping noise dominates library-size
# signal).
EXCLUDE_PREFIXES = ("RPL", "RPS")
