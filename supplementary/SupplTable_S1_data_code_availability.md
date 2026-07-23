# Suppl. Table S1 — Data and code availability (Figure 7)

Renderable content for the Supplementary data-and-code-availability table. All three data
accessions were verified to the correct dataset; the code/repository DOI must be minted
before submission (see note).

| Resource | Source | Accession / DOI | Verified |
|---|---|---|---|
| Tumor-infiltrating NK-cell atlas (Tang et al. 2023), processed | reference 1 | Zenodo **doi:10.5281/zenodo.8275845** | ✓ |
| Tang et al. 2023 atlas, raw sequencing | reference 1 | NCBI BioProject **PRJNA877828** (provided by the original authors upon request) | ✓ |
| CD8⁺ T-cell atlas, huARdb v2 (Xue et al. 2025), processed | reference 2 | Zenodo **doi:10.5281/zenodo.12542577** | ✓ |
| Analysis code + cached differential-expression / GSEA tables | this study | **repository archive DOI — to mint** | pending |

## Notes

- **Raw-data provenance:** the Tang raw sequencing data used here were provided by the original authors upon request; the CD8 huARdb v2 atlas (Xue et al. 2025) was obtained from its public Zenodo deposit.
- **Code DOI:** the analysis repository (pipeline code + the cached DEG, GSEA and audit
  tables that regenerate every panel) should be archived (e.g. a Zenodo snapshot of the tagged
  release) and the minted DOI substituted for the Methods "(DOI at submission)" placeholder.
- **huARdb DOI reconciliation:** the CD8 atlas used here is **huARdb v2**, described in **Xue et
  al., *Nat Methods* 2025;22(2):435-445 (doi:10.1038/s41592-024-02530-0)** with the processed
  data at **Zenodo 12542577** — that is the pairing to cite. The separate *Nucleic Acids Res*
  huARdb database paper (doi:10.1093/nar/gkae1003) describes the broader resource and is **not**
  the v2 atlas analyzed here; do not substitute it for reference 2.
- Direction convention and all per-context provenance are recorded in the other Supplementary
  items (Suppl. Data S1-S4).
