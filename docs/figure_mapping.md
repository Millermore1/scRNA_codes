# Figure and table inventory

Every output the pipeline produces, what it shows, and which module and
original script it comes from. Fill in the "manuscript panel" column once the
figures are fixed — that column is then the authoritative figure ↔ file map for
the submission checklist.

All figures are written twice: `.pdf` (vector, `cairo_pdf`) and `.png`
(300 dpi preview). Tables are UTF-8 CSV unless noted.

## Module 01 — data preparation

| Output (stem) | Shows | Original script(s) | Manuscript panel |
|---|---|---|---|
| `01_qc_violin_human` | QC metric distributions per cohort | `single cell-20240607.R`, `QC.R` | |
| `01_qc_violin_mouse` | QC metric distributions per group | `QC.R` | |
| `01_qc_summary_*.csv` | Cells and median QC metrics per sample | `single cell-20240607.R` | suppl. table |

## Module 02 — annotation

| Output (stem) | Shows | Original script(s) | Manuscript panel |
|---|---|---|---|
| `02_umap_celltype_human` | Cell classes of the human atlas | `正式code908.R`, `full_codes.R` | |
| `02_umap_celltype_by_group_human` | Same embedding, one facet per cohort | `full_codes.R`, `正式code908-2.R` | |
| `02_dotplot_markers_human` | Canonical markers per class | `marker gene.R`, `full_codes.R` | |
| `02_umap_prmt8_human` | PRMT8 across the atlas | `Heatmap-last-2.R` | |
| `02_marker_heatmap_human` | Top markers per class, Z-scored | `QC.R` | |
| `02_umap_celltype_mouse` | Cell classes of the mouse cKO atlas | `inhibit_6_29.1.R`, `doubleFinder_28.1.R` | |
| `02_umap_celltype_by_group_mouse` | Mouse atlas by group | `QC.R` | |
| `02_dotplot_markers_mouse` | Mouse marker dot plot | `doubleFinder_28.1.R` | |
| `02_markers_*.csv` | Cluster marker tables | `single cell-20240607.R` | suppl. table |

## Module 03 — composition

| Output (stem) | Shows | Original script(s) | Manuscript panel |
|---|---|---|---|
| `03_composition_human` | Stacked cell-class proportions per cohort | `细胞比例.R`, `inhibit_6_29.1.R` | |
| `03_composition_mouse` | Stacked proportions per mouse group | `inhibit_6_29.1.R` | |
| `03_celltype_counts_*` | Absolute cells per class | `full_codes.R` | |
| `03_cells_per_cohort_human` | Donut: cells contributed per cohort | `full_codes.R` | |
| `03_proportion_test_*.csv` | Sample-level Wilcoxon tests of proportions | new (replaces cell-level eyeballing) | suppl. table |

## Module 04 — differential expression

| Output (stem) | Shows | Original script(s) | Manuscript panel |
|---|---|---|---|
| `04_deg_mouse_all.csv` | All mouse DEGs, all classes, both contrasts | `Fig.6_D.R`, `最新的修改figure.R` | suppl. table |
| `04_deg_counts_mouse` | DEG counts per class × contrast | `inhibit_5_28.1.R` | |
| `04_volcano_inhibitory_KOF_vs_WT` | Inhibitory DEGs with ASD genes labelled | `inhibit_4.R` | |
| `04_deg_human_all.csv` | Disease vs control DEGs per class | new (generalises `single cell-20250811.R`) | suppl. table |
| `04_shared_deg_mouse_human_inhibitory.csv` | Genes significant in both species | new | suppl. table |

## Module 05 — inhibitory neurons

| Output (stem) | Shows | Original script(s) | Manuscript panel |
|---|---|---|---|
| `05_umap_inhibitory_clusters_mouse` | Interneuron subclusters | `inhibit.R` | |
| `05_dotplot_inhibitory_markers_mouse` | VIP/SST/PV/TH/LAMP5/PAX6 panels | `inhibit.R`, `cell_chat_5.09.R` | |
| `05_umap_inhibitory_subtypes_mouse` | PVALB / SST / VIP assignment | `inhibit.R`, `inhibit_3.R` | |
| `05_subtype_composition_mouse` | Subtype proportions per group | `inhibit_3.R` | |
| `05_prmt8_violin_mouse` | Prmt8 per subtype × group | `full_codes.R` (human analogue) | |
| `05_prmt8_test_mouse.csv` | Sample-level Wilcoxon per subtype | new | suppl. table |
| `05_deg_inhibitory_subtypes_mouse.csv` | Subtype-specific DEGs | `inhibit_2.R`, `inhibit_4.R` | suppl. table |
| `05_volcano_<subtype>_Gq_vs_vector` | Volcano per subtype | `inhibit_3.R`, `inhibit_4.R` | |
| `05_venn_inhibitory_deg_mouse` | Overlap of the two contrasts | `inhibit_5_28.1.R` | |
| `05_hp_autism_score_mouse` | HP_AUTISM module score per group | `inhibit_6_29.1.R` | |
| `05_hp_autism_heatmap_mouse` | Mean score, subtype × group | `inhibit_6_29.1.R` | |
| `05_prmt8_violin_human_inhibitory` | PRMT8 in human inhibitory neurons | `single cell-20250811.R` | |

## Module 06 — excitatory neurons

| Output (stem) | Shows | Original script(s) | Manuscript panel |
|---|---|---|---|
| `06_umap_excitatory_subtypes_mouse` | L2/3 – L6 layers | `excit_2.R`, `all_neurons_analyse.R` | |
| `06_dotplot_layer_markers_mouse` | Layer marker panel | `excit_2.R` | |
| `06_dotplot_excitatory_ordered_mouse` | Clusters ordered by marker score | `all_neurons_analyse.R` | |
| `06_prmt8_by_layer_mouse` | Prmt8 across layers and groups | `all_neurons_analyse.R` | |
| `06_layer_composition_mouse` | Layer proportions per group | `excit_2.R` | |
| `06_excitatory_cluster_markers.xlsx` | One sheet per cluster | `all_neurons_analyse.R` | suppl. table |

## Module 07 — enrichment

| Output (stem) | Shows | Original script(s) | Manuscript panel |
|---|---|---|---|
| `07_go_<contrast>_<dir>.csv` | GO BP terms per contrast and direction | `DEG_enrichment_network.R` | suppl. table |
| `07_go_bar_genes_<contrast>_<dir>` | GO bar with member genes inside | `GO_Figure_6.R`, `GO_Figure_6_2 (1).R` | |
| `07_go_network_<contrast>_<dir>` | Term network with cell-type pies | `DEG_enrichment_network.R`, `Fig.6_D.R` | |
| `07_go_signed_bar_KOF_vs_WT` | Up vs down terms side by side | `最新的修改figure.R` | |
| `07_kegg_*.csv` | KEGG pathways | `GO_Figure_6_2 (1).R` | suppl. table |

## Module 08 — CellChat

| Output (stem) | Shows | Original script(s) | Manuscript panel |
|---|---|---|---|
| `08_umap_neuron_subtypes` | PN / PVALB / SST / VIP | `cell_chat_5.09.R` | |
| `08_cellchat_circle_interactions` | Interaction counts, four groups, common scale | `cell chat.R` | |
| `08_cellchat_circle_strength` | Interaction strength, four groups | `cell chat.R` | |
| `08_cellchat_overall_interactions` | Total number and weight per group | `PRMT8KO_snRNA_seq_cellchat.R` | |
| `08_cellchat_diff_KOF_vs_WT` | Differential interaction network | `cell_chat_5.09.R` | |
| `08_cellchat_information_flow` | Pathway information flow, stacked | `cell chat-6.29.R` | |
| `08_cellchat_signaling_roles` | Outgoing vs incoming scatter per group | `cell_chat_5.09.R`, `code_for_cellchat.R` | |
| `08_cellchat_synapse_bubble` | Synapse pathways to PN and PVALB | `PRMT8KO_snRNA_seq_cellchat.R` | |
| `08_cellchat_interactions_all.csv` | Every inferred interaction | `cell chat.R` | suppl. table |

## Module 09 — WGCNA

| Output (stem) | Shows | Original script(s) | Manuscript panel |
|---|---|---|---|
| `09_wgcna_soft_powers` | Scale-free topology fit | `WGCNA.R` | suppl. figure |
| `09_wgcna_module_dendrogram` | Module assignment | `WGCNA.R` | |
| `09_wgcna_module_network_PRMT8` | PRMT8 module hub network | `WGCNA.R` | |
| `09_wgcna_ggraph_PRMT8` | Cleaner redrawn neighbourhood | `WGCNA.R` | |
| `09_wgcna_module_trait` | Module–disease correlation | `WGCNA.R` | |
| `09_wgcna_hub_genes_*.csv` | Hub genes with kME | `WGCNA.R` | suppl. table |

## Module 10 — heat maps

| Output (stem) | Shows | Original script(s) | Manuscript panel |
|---|---|---|---|
| `10_heatmap_ion_channels_human` | Ion-channel genes by cohort | `Heatmap-last.R`, `Heatmap-last-2.R` | |
| `10_heatmap_gaba_receptors_human` | GABA/glutamate receptors by cohort | `热图绘制0916.R` | |
| `10_heatmap_synaptic_prmt8_human` | Synaptic/ASD genes incl. PRMT8 | `热图绘制0916.R` | |
| `10_prmt8_heatmap_human` | PRMT8 per class × cohort | `full_codes.R` | |
| `10_prmt8_brain_region_correlation` | PRMT8 correlation across brain regions | `heatmap_correlation_region.R` | |
| `10_dotplot_markers_ordered_human` | Dot plot with fixed class order | `marker gene.R` | |

## Module 11 — ASD organoids

| Output (stem) | Shows | Original script(s) | Manuscript panel |
|---|---|---|---|
| `11_umap_timepoint_organoid` | Organoids by developmental stage | `ASD_S1.R`, `single cell-ASD.R` | |
| `11_umap_group_organoid` | ASD vs control | `single cell-ASD.R` | |
| `11_volcano_ASD_vs_Control_D30` | DE at D30 | `ASD_S1.R` | |
| `11_prmt8_trajectory_organoid` | PRMT8 mean expression per stage | `ASD_S1.R` | |
| `11_gsea_dotplot_organoid` | GSEA of PRMT8-high vs PRMT8-low | `single cell-ASD.R` | |

---

## Conventions

- Numbers in file names are the module that produced them, so a file can
  always be traced back to its script.
- Any figure that was drawn once and then re-edited in Illustrator should keep
  the exported PDF under version control and be re-exported from the script
  rather than re-edited by hand.
- `docs/migration_notes.md` lists which of the 43 original scripts fed which
  module.
