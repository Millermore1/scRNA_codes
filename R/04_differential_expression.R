# =============================================================================
# R/04_differential_expression.R
#
# Differential expression with Seurat's Wilcoxon test:
#   * mouse atlas: KOF vs WT and KOF_Gq_CNO vs KOF_vector_CNO, within each of
#     the seven major cell classes;
#   * human atlas: each disease cohort vs controls, within each cell class;
#   * a focused volcano plot for inhibitory neurons (PRMT8 cKO vs WT) with the
#     autism-risk genes labelled.
#
# Outputs
#   results/tables/04_deg_<dataset>_<scope>.csv
#   results/figures/04_deg_counts_<dataset>.pdf
#   results/figures/04_volcano_*.pdf
# =============================================================================

if (!file.exists("R/00_setup.R")) stop("Run this script from the repository root.")
source("R/00_setup.R")

need("Seurat", "dplyr", "ggplot2", "ggrepel", "patchwork")

# -----------------------------------------------------------------------------
# 1. Mouse atlas: chemogenetic and genetic contrasts
# -----------------------------------------------------------------------------

step_header("Mouse DEG: KOF vs WT and Gq vs vector, per cell type")

mouse <- read_seurat(p_processed("mouse_ko_annotated.rds"), "annotated mouse atlas")

deg_mouse <- run_deg_by_celltype(
  mouse,
  contrasts    = cfg$design$mouse_contrasts,
  group_col    = "group",
  celltype_col = "celltype"
)

if (!is.null(deg_mouse)) {
  deg_mouse$dataset <- "mouse"
  save_table(filter_deg(deg_mouse), "04_deg_mouse_significant", digits = 4)
  save_table(deg_mouse, "04_deg_mouse_all", digits = 4)

  counts_mouse <- count_deg(deg_mouse)
  save_table(counts_mouse, "04_deg_counts_mouse")

  if (nrow(counts_mouse) > 0) {
    save_figure(
      plot_deg_tile(counts_mouse),
      "04_deg_counts_mouse", width = 4.6, height = 4
    )
  }
}

# -----------------------------------------------------------------------------
# 2. Mouse atlas: inhibitory neurons in detail
# -----------------------------------------------------------------------------

step_header("Mouse DEG: inhibitory neurons, KOF vs WT")

mouse_inh <- subset(mouse, cells = colnames(mouse)[
  as.character(mouse$celltype) == "Inhibitory neurons"
])
note("inhibitory neurons: ", ncol(mouse_inh), " cells")

deg_inh <- run_deg_pair(mouse_inh, ident1 = "KOF", ident2 = "WT",
                        group_col = "group")
if (!is.null(deg_inh)) {
  deg_inh$dataset <- "mouse_inhibitory"
  save_table(deg_inh, "04_deg_inhibitory_KOF_vs_WT", digits = 4)
  save_table(filter_deg(deg_inh), "04_deg_inhibitory_KOF_vs_WT_significant",
             digits = 4)

  save_figure(
    plot_volcano(deg_inh,
                 label_genes = c("Prmt8", cfg$genes$asd_label_genes),
                 title = "Inhibitory neurons: KOF vs WT"),
    "04_volcano_inhibitory_KOF_vs_WT", width = 5.4, height = 4.6
  )

  # Gene-level summary used in the text: how many genes go up and down.
  sig_inh <- filter_deg(deg_inh)
  note("significant DEGs -- up: ", sum(sig_inh$direction == "Up"),
       ", down: ", sum(sig_inh$direction == "Down"))
}

# -----------------------------------------------------------------------------
# 3. Human atlas: disease vs control
# -----------------------------------------------------------------------------

step_header("Human DEG: each cohort vs controls, per cell type")

human <- read_seurat(p_processed("human_atlas_annotated.rds"),
                     "annotated human atlas")

deg_human <- run_deg_by_celltype(
  human,
  contrasts    = cfg$design$human_contrasts,
  group_col    = "group",
  celltype_col = "celltype"
)

if (!is.null(deg_human)) {
  deg_human$dataset <- "human"
  save_table(filter_deg(deg_human), "04_deg_human_significant", digits = 4)
  save_table(deg_human, "04_deg_human_all", digits = 4)

  counts_human <- count_deg(deg_human)
  save_table(counts_human, "04_deg_counts_human")
  if (nrow(counts_human) > 0) {
    save_figure(
      plot_deg_tile(counts_human),
      "04_deg_counts_human", width = 6, height = 4
    )
  }
}

# -----------------------------------------------------------------------------
# 4. Cross-species check: are PRMT8 targets shared?
# -----------------------------------------------------------------------------

step_header("Compare mouse and human inhibitory DEGs")

if (!is.null(deg_inh) && !is.null(deg_human)) {
  mouse_sig <- unique(filter_deg(deg_inh)$gene)

  human_sig_df <- filter_deg(deg_human)
  human_sig_df <- human_sig_df[human_sig_df$celltype == "Inhibitory.neurons", ]
  human_sig <- unique(human_sig_df$gene)

  # Mouse symbols are title-case, human symbols upper-case; harmonise.
  shared <- intersect(toupper(mouse_sig), human_sig)

  save_table(
    data.frame(gene = sort(shared), stringsAsFactors = FALSE),
    "04_shared_deg_mouse_human_inhibitory"
  )
  note(length(shared), " genes are significant in both species")
}

note("done: DEG tables and figures written to results/")
