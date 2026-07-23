# TUSC2 differential expression in tumor-infiltrating NK and CD8 T cells (Figure 7)

Self-contained reproducibility package for Figure 7: TUSC2-negative vs
TUSC2-positive differential expression and rank-based pathway enrichment in
tumour-infiltrating **NK cells** (Tang et al. 2023 pan-cancer NK atlas) and
**CD8 T cells** (huARdb v2, Xue et al. 2025).

The pipeline is **tumour-only** (one context per lineage: NK `Tang_Tumor`, CD8
`CD8_TIL`) and **enrichment is threshold-free GSEA-prerank only**.

The package ships the cached differential-expression and GSEA tables, so every
Figure 7 panel and the Supplementary workbook regenerate **offline, from cache,
with no external data and no network access**.

## Package layout

```
reproduce.py                 # offline: renders the 4 panels + builds the Supplementary xlsx
collect_supplementary.py     # gathers S1-S4 + Table S1 into supplementary/
make_supplementary_xlsx.py   # merges the collected files into one workbook
tusc2_deg/
├── paths.py                 # input-data paths (only needed for the from-raw path)
├── gene_sets/               # vendored MSigDB Hallmark .grp sets (panel B + GSEA)
├── nk/                      # NK pipeline  (Tang 2023)
├── cd8/                     # CD8 pipeline (huARdb v2, Xue 2025)
│   └── output/{degs,gsea,audit}/   # cached tables that drive every panel
└── publication/            # panel renderers (volcano, GSEA bar, running-ES)
tests/                       # offline test suite
figures/                     # written by reproduce.py — the 4 Figure 7 panels
supplementary/               # written by reproduce.py — S1-S4 + Table S1 + xlsx
```

`nk/` and `cd8/` are parallel pipelines with the same structure
(`config · data · degs · gsea · pipeline`). Each pipeline's `config.py` is the
single source of truth for its thresholds, gene-family taxonomy, palette and
paths. Both render their panels through `tusc2_deg.publication` so the NK and CD8
figures are visually uniform.

## Method

Per lineage, in its tumour context:

1. **pseudobulk** — sum raw counts per `(donor, TUSC2 status)`; keep groups with
   ≥ 10 cells and genes with summed counts ≥ 10. A donor enters only if it
   contributes ≥ 10 cells to both the TUSC2-positive and TUSC2-negative arms, so
   each retained donor supplies a matched pair of pseudosamples.
2. **degs** — donor-paired pyDESeq2 (`~ donor + TUSC2 status`), median-of-ratios
   size factors, empirical-Bayes dispersions (Love et al. 2014), apeglm LFC
   shrinkage (Zhu et al. 2019). Contrast: TUSC2-negative vs TUSC2-positive. A
   per-cell Wilcoxon table is written as a confirmatory cross-check.
3. **gsea** — GSEA-prerank (gseapy) on the signed Wald statistic, against
   GO Biological Process 2025 + Reactome Pathways 2024 + MSigDB Hallmark 2020
   (size filter 15–500, 1,000 permutations, fixed seed). This is the sole,
   threshold-free enrichment readout.

Alongside the stages, the pipeline writes two audit tables from the live DEG
frames: the effect-size significance ladder (`deg_sig_ladder.tsv`) and the
curated volcano-label provenance (`curated_label_provenance.tsv`).

**Direction convention.** A positive log₂ fold change or NES denotes higher
expression in TUSC2-negative cells; a negative value denotes higher expression
in TUSC2-positive cells.

**Significance gate** (three-part): `padj < 0.05`, `baseMean > 5`,
`|log2FC| > 0.1`; robustness is reported across an effect-size ladder
(`0.1 / 0.25 / 0.585 = log₂1.5`) and against a calibrated 1.1-fold null.

The full Methods and panel legends appear in the accompanying manuscript (Figure 7).

## Install

Python **3.13**. From this folder:

```bash
pip install -e .
# or, for pinned versions used to produce the cached tables:
pip install -r requirements.txt
```

## Offline reproduce (default path — no external data)

```bash
MPLBACKEND=Agg python reproduce.py
```

This renders the four Figure 7 panels into `figures/` and builds the
Supplementary data files + workbook into `supplementary/`, entirely from the
cached tables under `tusc2_deg/{nk,cd8}/output`. No atlas download, no network.

Outputs:

| Path | Contents |
|------|----------|
| `figures/fig7_a_gsea_nes_bar.{png,pdf}` | GSEA Hallmark NES bars (NK + CD8) |
| `figures/fig7_b_oxphos_runningES.{png,pdf}` | OXPHOS running-enrichment (NK + CD8) |
| `figures/fig7_c_cd8_volcano.{png,pdf}` | CD8 TIL volcano |
| `figures/fig7_d_nk_volcano.{png,pdf}` | NK Tang tumour volcano |
| `supplementary/SupplData_S1..S4_*` | per-gene DE, effect-size ladder, GSEA, label provenance |
| `supplementary/SupplTable_S1_*` | data and code availability |
| `supplementary/Figure7_Supplementary_Data.xlsx` | all of the above, one sheet per file |

Panel B is recomputed deterministically from the cached ranking. The GSEA bar and
volcanos read the cached DE/GSEA tables. Volcano gene-label placement uses
`adjustText`, whose iterative solver introduces sub-pixel jitter in label
positions between runs; the data, axes and colours are fully deterministic.

## Optional: rebuild the cache from raw atlases

Only needed to regenerate the cached tables from scratch. Download the processed
atlases, point `TUSC2_DATA_DIR` at them, then run each pipeline:

```bash
export TUSC2_DATA_DIR=/path/to/atlases
python -m tusc2_deg.nk.pipeline
python -m tusc2_deg.cd8.pipeline
```

Flags: `--force` recomputes every stage; `--force-from <stage>` recomputes from a
named stage onward.

Accessions:

| Dataset | Reference | Accession |
|---------|-----------|-----------|
| Tumour-infiltrating NK-cell atlas | Tang et al. 2023 | Zenodo `10.5281/zenodo.8275845` (raw: NCBI BioProject `PRJNA877828`) |
| CD8 T-cell atlas, huARdb v2 | Xue et al. 2025, *Nat Methods* | Zenodo `10.5281/zenodo.12542577` |

## Tests

```bash
MPLBACKEND=Agg python -m pytest -q
```

The suite runs offline against the cached tables and injectable renderers. A
single `slow`-marked single-cell cross-check loads a raw atlas if one is present
and otherwise skips.

## Citation

This analysis code is archived at Zenodo: **10.5281/zenodo.21516067** (concept
DOI — always resolves to the latest version). Please cite it alongside the
accompanying publication.

## License

MIT — see `LICENSE`.
