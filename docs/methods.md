# Analysis methods

This document records every processing step and parameter that produced the
figures. It is written so that the "Single-cell RNA-seq analysis" subsection of
the Methods can be assembled from it directly, and so that a reviewer can see
the complete parameter set without reading the code.

Parameter values live in [`config/config.R`](../config/config.R); the section
below names the corresponding `cfg` entry.

---

## 1. Datasets

### 1.1 Human multi-cohort brain snRNA-seq atlas

Single-nucleus RNA-seq from human brain tissue, assembled from public cohorts
and analysed as one object. Cohorts and sample identifiers are listed in
[`data/metadata/cohort_mapping.csv`](../data/metadata/cohort_mapping.csv):

| Cohort | Accession(s) | Samples |
|--------|--------------|---------|
| Alzheimer's disease (AD) | GSE157827, GSE161045 | 16 |
| Parkinson's disease (PD) | GSE140231, GSE161045 | 14 |
| Multiple sclerosis (MS) | MS cohort | 4 |
| Epilepsy (EPI) | GSE302285 | 6 |
| Neurological controls (NC) | GSE157827, control cohort | 14 |

The number of cells per cohort used for the composition figure was AD 109,037;
EPI 34,083; MS 11,208; NC 100,777; PD 31,919 (taken from the assembled object;
`cfg$design$human_group_order`).

`TODO (author)`: add the sample-level table (sex, age, post-mortem interval,
brain region, library protocol) as a supplementary table, and confirm the
number of donors per cohort after QC.

### 1.2 Mouse Prmt8 conditional-knockout atlas

snRNA-seq of mouse brain from four experimental groups
([`data/metadata/mouse_groups.csv`](../data/metadata/mouse_groups.csv)):

| Group | Genotype / manipulation |
|-------|-------------------------|
| `WT` | wild-type littermate control |
| `KOF` | *Prmt8* conditional knockout |
| `KOF_vector_CNO` | cKO + control virus + CNO (chemogenetic control) |
| `KOF_Gq_CNO` | cKO + hM3Dq + CNO (Gq-DREADD activation) |

### 1.3 Human cortical organoids (GSE203201)

ASD patient-derived and control organoids profiled at D0, D6, D12, D30, D60
and D100. Group labels (`ASD`, `Control`) are parsed from the sample names
(`cfg$params$organoid`).

---

## 2. Quality control

Per sample, before merging (`apply_qc_filter()`, `cfg$params$qc`):

| Metric | Threshold |
|--------|-----------|
| Detected genes (`nFeature_RNA`) | 300 – 5,000 |
| UMI counts (`nCount_RNA`) | > 1,000 and below the 97th percentile of the sample |
| Mitochondrial reads | < 10 % (< 5 % for the mouse object) |
| Haemoglobin reads | < 3 % (human samples) |

The upper UMI bound is computed within each sample so that libraries with
different sequencing depth are not penalised globally. Per-sample QC statistics
are written to `results/tables/01_qc_summary_*.csv`.

`TODO (author)`: state whether doublets were removed before or after QC, and
which tool was used (DoubletFinder was installed but its call is not part of
these scripts).

---

## 3. Normalisation, integration and clustering

- `NormalizeData()` — LogNormalize, scale factor 10,000.
- `FindVariableFeatures()` — vst, 3,000 genes (2,000 for the sub-clustering
  steps).
- `ScaleData()`, `RunPCA(npcs = 30)`.
- Batch correction: `harmony::RunHarmony(group.by.vars = "orig.ident")`.
- Graph: `FindNeighbors(dims = 1:10)` for the atlases; `1:20` for the
  interneuron/excitatory sub-clustering; `1:30` for the neuron compartment.
- Clustering: Louvain (`FindClusters`), resolution 0.05 for the atlases,
  0.1 for the inhibitory and excitatory sub-clustering.
- Embedding: UMAP on the Harmony space.

*(`integrate_atlas()`, `run_seurat_standard()`, `cfg$params$integration`,*
*`cfg$params$inhibitory`, `cfg$params$excitatory`,*
*`cfg$params$neuron_subclustering`)*

---

## 4. Cell-type annotation

Clusters were annotated manually from canonical markers
([`data/metadata/celltype_markers.csv`](../data/metadata/celltype_markers.csv)),
using dot plots of the marker panels
(`cfg$params$annotation$human_markers`, `$mouse_markers`).

Major classes: excitatory neurons, inhibitory neurons, oligodendrocytes
(+ OPC in the mouse), astrocytes, microglia, endothelial cells. Inhibitory
neurons were further split into PVALB, SST and VIP subtypes; excitatory neurons
into cortical layers L2/3, L4, L5 and L6 using layer markers
(`cfg$genes$excit_layers`).

Cluster→label mappings are explicit vectors in `cfg$design`
(`human_cluster_mapping`, `mouse_inhibitory_mapping`,
`mouse_excitatory_mapping`, `mouse_neuron_mapping`). Each mapping is validated
against the number of clusters produced by the run; a mismatch stops the
pipeline with an actionable message rather than mislabelling cells.

---

## 5. Cell-type composition

Proportions are computed per sample (donor) and then averaged within groups,
because cells from one donor are not independent observations
(`proportion_by_sample()` in `03_cell_composition.R`). Group differences are
tested with a two-sided Wilcoxon rank-sum test on the per-sample proportions,
with Benjamini–Hochberg correction across cell types
(`test_proportion_by_group()`).

---

## 6. Differential expression

- Test: Wilcoxon rank-sum (`test.use = "wilcox"`), implemented in
  `Seurat::FindMarkers()`.
- Minimum expression: `min.pct = 0.1`; log fold-change threshold 0.25 for the
  test, reporting threshold `|log2FC| > 0.25`.
- Significance: Benjamini–Hochberg adjusted p < 0.05.
- Cell types with fewer than 5 cells in either group are skipped and logged.
- Contrasts: KOF vs WT, and KOF_Gq_CNO vs KOF_vector_CNO (mouse); each disease
  cohort vs NC (human); ASD vs control at D30 (organoids).
  `ident1` is the numerator, so positive log2FC means higher in `ident1`.

*(`run_deg_by_celltype()`, `filter_deg()`; `cfg$params$de`,
`cfg$design$*_contrasts`)*

---

## 7. Functional enrichment

- Over-representation analysis with `clusterProfiler::enrichGO()`
  (Biological Process), background = all genes detected in the object,
  `pAdjustMethod = "BH"`, `pvalueCutoff = 0.05`, `qvalueCutoff = 0.2`.
- Redundancy reduction with `simplify(cutoff = 0.6, by = "p.adjust")`.
- Organism-specific annotation: `org.Mm.eg.db` (mouse), `org.Hs.eg.db` (human).
- KEGG: `enrichKEGG(organism = "mmu")` with the same thresholds.
- Term networks: the 10 most significant terms are nodes; two terms are
  connected when they share at least one gene and Jaccard similarity ≥ 0.05;
  node pies show the cell-type composition of the term's genes.
- The "GO bar with genes" panel uses fold enrichment
  (GeneRatio / BgRatio) as the bar length, with the term and up to eight member
  genes printed inside the bar.

*(`07_enrichment_networks.R`, `run_enrich_go()`, `cfg$params$enrichment`)*

`TODO (author)`: if the published panels used Metascape Z-scores rather than
fold enrichment, pass that table to `plot_go_bar_with_genes(df, score = "Z.score")`
— the layout function is independent of how the score was computed.

---

## 8. Cell–cell communication (CellChat)

- Object: neuron compartment (PN, PVALB, SST, VIP) of the mouse atlas.
- Database: `CellChatDB.mouse` (full, not restricted to secreted signalling).
- Inference per group: `computeCommunProb(type = "triMean")`,
  `filterCommunication(min.cells = 10)`, `computeCommunProbPathway()`,
  `aggregateNet()`, `netAnalysis_computeCentrality()`.
- Group comparison: `mergeCellChat()`, then overall interaction counts and
  weights (`compareInteractions`), differential interactions
  (`netVisual_diffInteraction`), information flow (`rankNet`), signalling roles
  (`netAnalysis_signalingRole_scatter`) and a synapse-pathway bubble plot.
- Circle plots are drawn with a common edge-weight maximum across groups so
  that line thickness is comparable between panels.

*(`08_cellchat.R`, `cfg$params$cellchat`)*

---

## 9. Co-expression network (hdWGCNA)

- Object: inhibitory neurons of the human atlas.
- Gene selection: expressed in ≥ 5 % of the metacells (`fraction = 0.05`).
- Metacells: `k = 25`, `max_shared = 10`, built within group × donor.
- Network: signed, soft threshold chosen from the scale-free topology fit
  (`TestSoftPowers`; the value used is `cfg$params$wgcna$soft_power`).
- Module detection with the default dynamic tree cut; module eigengenes and
  kME computed per module.
- Module–trait correlation against a binary disease indicator
  (`PD_Status`), with the trait definition in `cfg$params$wgcna$trait`.

*(`09_wgcna.R`)*

---

## 10. Module score (HP_AUTISM)

`Seurat::AddModuleScore()` with the autism risk-gene list in
[`data/metadata/hp_autism_genes.txt`](../data/metadata/hp_autism_genes.txt)
(265 genes, human symbols). The function matches gene-name style automatically
(human upper-case vs mouse title-case) and reports how many genes were found;
fewer than 10 matched genes raises a warning.

*(`add_module_score_auto()` in `R/utils/analysis.R`)*

---

## 11. Software

| Component | Version |
|-----------|---------|
| R | 4.3.0 |
| Seurat | 5.3.0 |
| harmony | 1.2.0 |
| CellChat | 1.6.1 |
| clusterProfiler | 4.10.0 |
| ComplexHeatmap | 2.18.0 |
| hdWGCNA | 0.4.12 |
| WGCNA | 1.74 |
| DoubletFinder | 2.0.6 |

The full list is in [`env/package_versions.tsv`](../env/package_versions.tsv);
`save_session_info()` writes the exact environment of any run to
`results/session_info.txt`.

---

## 12. Reproducibility notes

- `cfg$params$seed = 42` is set in `R/00_setup.R`; it controls the layout
  seeds of the force-directed plots and the heat-map sampling.
- Every figure is written as a `cairo_pdf` (vector, editable in Illustrator)
  plus a 300 dpi PNG preview.
- Modules are numbered in execution order; each writes an `.rds` that the next
  one reads, so the pipeline can be resumed from any step.
