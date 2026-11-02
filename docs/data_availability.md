# Data availability

## Datasets analysed

| Dataset | Cohort | Accession / source | Availability in this repository |
|---|---|---|---|
| Human brain snRNA-seq | Alzheimer's disease | GSE157827, GSE161045 | Public (GEO); accessions only |
| Human brain snRNA-seq | Parkinson's disease | GSE140231, GSE161045 | Public (GEO); accessions only |
| Human brain snRNA-seq | Multiple sclerosis | MS cohort (sample IDs `MS121`, `MS122`, `MS176`, `MS242`) | Public; accessions only |
| Human brain snRNA-seq | Epilepsy | GSE302285 | Public (GEO); accessions only |
| Human brain snRNA-seq | Neurological controls | GSE157827, plus control samples (`CO14`, `CO25`, `CO28`, `CO39`, `SD48.16`) | Public (GEO); accessions only |
| Human cortical organoids | ASD vs control, D0–D100 | GSE203201 | Public (GEO); accessions only |
| Human brain snRNA-seq | Additional cohort, h5 matrix | GSE174367 | Public (GEO); reader in `01_prepare_datasets.R` |
| Human brain snRNA-seq | Additional cohort, MatrixMarket files | GSE186537 | Public (GEO); reader in `01_prepare_datasets.R` |
| Mouse brain snRNA-seq | *Prmt8* cKO and DREADD groups | This study | **To be deposited**; see below |
| PRMT8 correlation matrix | Brain-region correlation of PRMT8 and partner genes | `correlations_PRMT8_brain_region.xlsx` | Values are from a public resource — cite the source |

The exact sample → cohort mapping used in the analysis is in
[`data/metadata/cohort_mapping.csv`](../data/metadata/cohort_mapping.csv)
(54 samples). Original scripts referenced these cohorts by accession prefix;
the CSV is the reviewed version of that mapping.

## What is not in this repository

Raw sequencing data are not stored here — they are too large for Git and are
already public for the human cohorts. The repository contains code,
configuration, metadata and generated tables/figures only:

```
data/raw/          <- put downloaded matrices here (git-ignored)
data/processed/    <- Seurat objects written by the pipeline (git-ignored)
results/           <- figures and tables (git-ignored except the manifest)
```

## Before submission

1. **Deposit the mouse data.** Upload the raw FASTQ or Cell Ranger
   `filtered_feature_bc_matrix` files for all four groups to GEO/SRA and insert
   the accession below. Reviewers routinely check this for a cKO + DREADD
   study; a "data available on request" statement is often flagged.

   > Mouse snRNA-seq data have been deposited in GEO under accession
   > **GSE______**.

2. **Confirm the MS cohort.** The four `MS*` samples are cited without an
   accession in the analysis scripts. Add the source publication/accession, or
   state that they were generated for this study.

3. **Add the processed object.** If the journal requires it, deposit the
   annotated `.rds` files (human atlas, mouse atlas, organoids) in a repository
   that supports large files (Zenodo, Figshare) and cite the DOI. The `.rds`
   files produced here are `data/processed/*.rds`.

4. **Gene list provenance.** `data/metadata/hp_autism_genes.txt` holds the
   autism risk-gene list used for the module score (265 genes, human symbols).
   Record its source (database and release, or the publication it was compiled
   from) in the Methods.

5. **Correlation matrix provenance.** State which dataset produced
   `correlations_PRMT8_brain_region.xlsx` (brain region, platform, number of
   samples), and whether the reported values are Pearson or Spearman
   correlations.

## Suggested data availability statement

> Single-nucleus RNA-seq data for the human cohorts were obtained from GEO
> under accessions GSE157827, GSE161045, GSE140231, GSE302285, GSE203201,
> GSE174367 and GSE186537. Mouse *Prmt8* cKO and chemogenetic snRNA-seq data
> generated in this study have been deposited in GEO under accession
> GSE______. All analysis code, configuration files and the metadata needed to
> reproduce every figure are available at
> https://github.com/<user>/PRMT8-inhibitory-neuron-scRNAseq (tag
> `v1.0.0`).

## Suggested code availability statement

> Analysis code is available at
> https://github.com/<user>/PRMT8-inhibitory-neuron-scRNAseq under the MIT
> licence. The repository is archived at Zenodo
> (DOI: 10.5281/zenodo.______). R 4.3.0 was used with the package versions
> listed in `env/package_versions.tsv`.
