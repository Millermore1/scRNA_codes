# PRMT8 in inhibitory neurons — single-nucleus RNA-seq analysis

Analysis code for the study of **PRMT8 in inhibitory neurons across
neurological disease and in a mouse conditional-knockout model with
chemogenetic manipulation**.

The repository contains the complete, re-runnable pipeline behind the
figures: from Cell Ranger output to the publication panels. Raw data are not
included (see [`docs/data_availability.md`](docs/data_availability.md)).

---

## What is analysed

| Dataset | Design | Role in the study |
|---|---|---|
| Human multi-cohort brain snRNA-seq | AD, PD, MS, EPI and controls (54 samples, 9 GEO series) | PRMT8 expression and co-expression across disease |
| Mouse *Prmt8* cKO + DREADD snRNA-seq | `WT`, `KOF`, `KOF_vector_CNO`, `KOF_Gq_CNO` | causal test of *Prmt8* loss and Gq activation |
| Human cortical organoids (GSE203201) | ASD vs control, D0–D100 | developmental trajectory of PRMT8 |

Analyses: QC and Harmony integration, cell-type annotation, composition,
per-cell-type differential expression, GO/KEGG enrichment with term networks,
CellChat communication, hdWGCNA co-expression modules, module scoring and a
large set of publication figures.

---

## Repository layout

```
.
├── README.md                  ← you are here
├── README.zh-CN.md            中文说明
├── run_all.R                  run every module, or a subset
├── config/
│   ├── config.R               paths, thresholds, contrasts, gene panels
│   └── palettes.R             colour schemes and level orders
├── R/
│   ├── 00_setup.R             entry point: root, config, packages
│   ├── 01_prepare_datasets.R  read + QC + integrate (human, mouse)
│   ├── 02_annotation.R        cell types, marker plots, UMAP panels
│   ├── 03_cell_composition.R  proportions + sample-level tests
│   ├── 04_differential_expression.R
│   ├── 05_inhibitory_neurons.R
│   ├── 06_excitatory_neurons.R
│   ├── 07_enrichment_networks.R
│   ├── 08_cellchat.R
│   ├── 09_wgcna.R
│   ├── 10_figures_heatmaps.R
│   ├── 11_asd_organoid.R
│   └── utils/                 io, analysis and plotting helpers
├── data/
│   ├── metadata/              cohort mapping, marker panel, gene sets
│   ├── raw/                   your input data (git-ignored)
│   └── processed/             Seurat objects written by the pipeline (git-ignored)
├── results/
│   ├── figures/               PDF + PNG for every panel (git-ignored)
│   ├── tables/                CSV tables (git-ignored)
│   └── objects/               intermediate .rds files (git-ignored)
├── env/                       package installation + tested versions
├── docs/
│   ├── methods.md             every parameter, ready for the Methods section
│   ├── figure_mapping.md      figure ↔ script ↔ output
│   ├── data_availability.md   accessions and data statements
│   └── migration_notes.md     how the original 43 scripts map onto this repo
├── tests/
│   ├── smoke_test.R           data-free self-check (runs in CI)
│   └── integration_test.R     exercises the helpers on a synthetic object
├── .github/workflows/ci.yml   runs the smoke test on push
└── archive/original_scripts/  the 43 original scripts, kept for provenance
```

---

## Quick start

### 1. Install the packages

```r
source("env/install_packages.R")
```

Tested with R 4.3.0, Seurat 5.3.0, CellChat 1.6.1, clusterProfiler 4.10.0
(full list in [`env/package_versions.tsv`](env/package_versions.tsv)).

### 2. Put the data in place

Either copy the inputs into `data/raw/` (see the default paths in
`config/config.R`), or point the pipeline at them with environment variables —
this is the recommended way if the data live on another drive:

```r
# add to ~/.Renviron
PRMT8_INPUT_DIR   = "D:/scrnadata"
# or set every location individually
PRMT8_HUMAN_10X_DIR = "D:/scrnadata/human_10x"
PRMT8_MOUSE_ATLAS_RDS = "D:/scrnadata/scRNA_HKU.rds"
PRMT8_ASD_10X_DIR   = "D:/scrnadata/GSE203201"
PRMT8_BRAIN_COR_XLSX = "D:/scrnadata/correlations_PRMT8_brain_region.xlsx"
```

### 3. Run

```bash
Rscript run_all.R          # everything, in order
Rscript run_all.R 01 02 05 # only these modules
Rscript run_all.R --list   # what each module does
```

or from RStudio, run the modules one by one (`R/01_...`, `R/02_...`). Every
module is self-contained and reads the `.rds` written by the previous one, so
you can stop and resume anywhere.

Before the first run, check that everything is wired correctly:

```bash
Rscript tests/smoke_test.R        # no data, no Seurat needed (~2 s)
Rscript tests/integration_test.R  # needs Seurat + ComplexHeatmap (~20 s)
```

The smoke test parses every script, validates the configuration and palettes,
and exercises the helper functions on toy data. The integration test builds a
small synthetic Seurat object and runs the shared analysis and plotting
functions end to end, so interface mistakes are caught without the real data.

---

## Preparing the input objects

Two objects must exist before module 02:

- **`data/processed/human_atlas_integrated.rds`** — built by module 01 from
  Cell Ranger directories, or supplied by you if the atlas was integrated
  elsewhere. Required metadata: `orig.ident`, `group`, `seurat_clusters`, a
  `umap` reduction.
- **`data/processed/mouse_ko_integrated.rds`** — the mouse cKO atlas. Export it
  from your existing object (the old `scRNA_HKU`):

```r
library(Seurat)
scRNA_HKU$group <- factor(scRNA_HKU$group,
                          levels = c("WT", "KOF", "KOF_vector_CNO", "KOF_Gq_CNO"))
saveRDS(scRNA_HKU, "data/processed/mouse_ko_integrated.rds")
```

Required metadata: `orig.ident`, `group`, `celltype`, `seurat_clusters`, a
`umap` reduction. If `mt_percent` is missing it is computed by module 01.

The organoid cohort (GSE203201) is read from `data/raw/GSE203201/`; if the
files still carry GEO prefixes, module 11 renames them into per-sample
Cell Ranger folders automatically.

---

## Outputs

Every figure is written twice — a `cairo_pdf` for editing and submission, and a
300 dpi PNG for quick inspection:

```
results/figures/<module>_<name>.pdf
results/figures/<module>_<name>.png
results/tables/<module>_<name>.csv
```

The file-name prefix is the module number, so any output can be traced back to
its script. [`docs/figure_mapping.md`](docs/figure_mapping.md) lists every
output with the corresponding manuscript panel (to be filled in).

---

## Reproducibility

- All thresholds, contrasts and gene panels are in `config/config.R`; scripts
  never hard-code a path or a cut-off.
- Cluster → cell-type mappings are validated against the number of clusters and
  fail loudly instead of mislabelling cells.
- `save_session_info()` writes the exact package versions of every run to
  `results/session_info.txt`.
- `set.seed(cfg$params$seed)` (42) is applied in `R/00_setup.R`.
- The `archive/` folder keeps the original scripts untouched, so the history of
  every analysis step is preserved.

---

## Publishing to GitHub

```bash
cd PRMT8-inhibitory-neuron-scRNAseq
git init
git add .
git commit -m "Initial commit: PRMT8 interneuron snRNA-seq analysis pipeline"

# create an empty repository on GitHub first, then:
git branch -M main
git remote add origin https://github.com/MillerMore1/scRNA_codes.git
git push -u origin main

git tag -a v1.0.0 -m "Version accompanying the manuscript"
git push origin v1.0.0
```

Before pushing, replace the placeholders:

| Where | What |
|---|---|
| `CITATION.cff` | full name, ORCID, affiliation, email |
| `LICENSE` | copyright holder |
| `README.md`, `README.zh-CN.md` | contact email |
| `docs/data_availability.md` | accession for the mouse data, repository URL |

After the first push, archive the tagged release on Zenodo to obtain a DOI and
add it to `docs/data_availability.md`. The smoke test runs automatically on
every push (`.github/workflows/ci.yml`).

## Citation

If you use this code, please cite the paper (details in
[`CITATION.cff`](CITATION.cff)) and this repository. Please also cite the
original datasets listed in
[`docs/data_availability.md`](docs/data_availability.md).

## Licence

Code: [MIT](LICENSE). Figures and tables derived from controlled-access data
remain subject to the terms of the original data providers.

## Contact

William Huang — MillerMore1@users.noreply.github.com — issues and pull requests are welcome.
