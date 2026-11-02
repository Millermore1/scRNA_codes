# data/

| Folder | Holds | In git? |
|---|---|---|
| `raw/` | downloaded matrices (Cell Ranger folders, `.h5`, MatrixMarket triplets, the correlation workbook) | no — too large, and public data should be fetched from the source |
| `processed/` | Seurat objects written by `R/01`–`R/05` | no — regenerate by running the pipeline |
| `metadata/` | small, human-readable tables that belong with the paper | **yes** |

## `metadata/`

- `cohort_mapping.csv` — every sample of the human atlas with its cohort
  (`NC`, `AD`, `PD`, `MS`, `EPI`) and accession. 54 samples. Read by
  `assign_cohort_groups()` in `R/01_prepare_datasets.R`.
- `mouse_groups.csv` — the four mouse experimental groups and what each
  manipulation was.
- `celltype_markers.csv` — the canonical marker panel behind the annotation,
  species by species. Suitable as a supplementary table.
- `hp_autism_genes.txt` — the 204-gene autism risk list used for the module
  score, one symbol per line.

## What goes in `raw/`

```
data/raw/
├── human_10x/<sample>/          # barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz
├── GSE203201/<sample>/          # same layout; the pipeline can rename GEO files
├── correlations_PRMT8_brain_region.xlsx
└── MSCtr_snRNA_ExpressionMatrix_R.txt
```

Set `PRMT8_INPUT_DIR` to a different root if these live elsewhere; see the
"Quick start" section of `README.md`.

The mouse atlas is a single `.rds` and belongs in `processed/`, because it is
the starting point of module 02 rather than a raw matrix. Instructions for
exporting it are in `README.md`.
