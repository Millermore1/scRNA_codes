# =============================================================================
# tests/integration_test.R
#
# Exercised the shared functions on a small synthetic Seurat object, so that
# interface mistakes (wrong argument names, missing columns, broken ggplot
# layers) are caught without needing the real data.
#
# Run with:  Rscript tests/integration_test.R
# Output goes to a temporary directory; nothing is written to results/.
# =============================================================================

local({
  args <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", args, value = TRUE)
  if (length(f) > 0) setwd(dirname(dirname(normalizePath(sub("^--file=", "", f[1])))))
})

source("R/00_setup.R")
need("Seurat", "ggplot2", "ComplexHeatmap", "circlize", "Matrix")

# Redirect every output to a scratch directory.
scratch <- file.path(tempdir(), "prmt8_test")
unlink(scratch, recursive = TRUE)
cfg$dirs$figures   <- file.path(scratch, "figures")
cfg$dirs$tables    <- file.path(scratch, "tables")
cfg$dirs$objects   <- file.path(scratch, "objects")
cfg$dirs$processed <- file.path(scratch, "processed")
for (d in cfg$dirs[c("figures", "tables", "objects", "processed")]) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

failures <- 0L
checks <- 0L
check <- function(label, condition) {
  checks <<- checks + 1L
  ok <- isTRUE(condition)
  if (!ok) failures <<- failures + 1L
  cat(sprintf("%s %s\n", if (ok) "[ ok ]" else "[FAIL]", label))
  invisible(ok)
}

# -----------------------------------------------------------------------------
# Synthetic object: 3 cell types x 2 groups, with real structure
# -----------------------------------------------------------------------------

set.seed(1)
n_cells <- 180
n_genes <- 400

mat <- matrix(stats::rpois(n_genes * n_cells, lambda = 2),
              nrow = n_genes, ncol = n_cells,
              dimnames = list(paste0("Gene", seq_len(n_genes)),
                              paste0("cell", seq_len(n_cells))))

ct <- rep(c("A", "B", "C"), length.out = n_cells)
grp <- rep(c("WT", "KOF"), each = n_cells / 2)

# Make marker genes and a group effect so DE and clustering have signal.
for (i in seq_len(n_genes)) {
  if (i %% 50 == 1) mat[i, ct == "A"] <- mat[i, ct == "A"] + 20
  if (i %% 50 == 11) mat[i, ct == "B"] <- mat[i, ct == "B"] + 20
  if (i %% 50 == 21) mat[i, ct == "C"] <- mat[i, ct == "C"] + 20
}
mat[1, grp == "KOF"] <- mat[1, grp == "KOF"] + 30   # Gene1 up in KOF

obj <- Seurat::CreateSeuratObject(mat, min.cells = 0, min.features = 0)
obj$group      <- grp
obj$celltype   <- ct
obj$orig.ident <- rep(c("s1", "s2", "s3", "s4"), length.out = n_cells)

obj <- Seurat::NormalizeData(obj, verbose = FALSE)
obj <- Seurat::FindVariableFeatures(obj, nfeatures = 200, verbose = FALSE)
obj <- Seurat::ScaleData(obj, verbose = FALSE)
obj <- Seurat::RunPCA(obj, npcs = 10, verbose = FALSE)
obj <- Seurat::FindNeighbors(obj, dims = 1:10, verbose = FALSE)
obj <- Seurat::FindClusters(obj, resolution = 0.5, verbose = FALSE)
obj <- Seurat::RunUMAP(obj, dims = 1:10, verbose = FALSE)

cat("\n== object handling ====================================================\n")

check("ensure_joined_layers returns a Seurat object",
      inherits(ensure_joined_layers(obj), "Seurat"))
check("cells_per_ident counts every cell",
      sum(cells_per_ident(obj, "celltype")) == n_cells)

cat("\n== relabel_clusters ===================================================\n")

mapping <- c("A", "B", "C")   # one label per level, matched by sorted level
lab <- relabel_clusters(obj, mapping, new_col = "celltype2",
                        cluster_col = "celltype")
check("relabel_clusters writes the new column",
      all(c("A", "B", "C") %in% as.character(lab$celltype2)))
check("relabel_clusters keeps the mapping order as factor levels",
      identical(levels(lab$celltype2), c("A", "B", "C")))
check("relabel_clusters sets the identities",
      identical(as.character(Seurat::Idents(lab)),
                as.character(lab$celltype2)))

bad <- tryCatch({
  relabel_clusters(obj, c("A", "B"), new_col = "x", cluster_col = "celltype")
  FALSE
}, error = function(e) TRUE)
check("relabel_clusters rejects a wrong-length mapping", bad)

named <- c(A = "Alpha", B = "Beta", C = "Gamma")
lab2 <- relabel_clusters(obj, named, new_col = "celltype3",
                         cluster_col = "celltype")
check("relabel_clusters accepts an explicit named mapping",
      all(as.character(lab2$celltype3) %in% named))

cat("\n== differential expression ===========================================\n")

deg <- run_deg_pair(obj, ident1 = "KOF", ident2 = "WT", group_col = "group")
check("run_deg_pair returns a data frame", is.data.frame(deg))
check("run_deg_pair keeps the expected columns",
      all(c("gene", "celltype", "comparison", "ident1", "ident2",
            "avg_log2FC", "p_val_adj") %in% colnames(deg)))
check("run_deg_pair finds the spiked-in gene",
      "Gene1" %in% deg$gene && deg$avg_log2FC[deg$gene == "Gene1"] > 1)

deg_ct <- run_deg_by_celltype(
  obj,
  contrasts = data.frame(comparison = "KOF_vs_WT", ident1 = "KOF",
                         ident2 = "WT", stringsAsFactors = FALSE),
  group_col = "group", celltype_col = "celltype"
)
check("run_deg_by_celltype covers all three cell types",
      length(unique(deg_ct$celltype)) == 3)

sig <- filter_deg(deg_ct)
check("filter_deg returns only significant rows",
      all(sig$p_val_adj < cfg$params$de$padj_cutoff))
check("count_deg returns one row per cell type and direction",
      nrow(count_deg(deg_ct)) >= 3)

cat("\n== composition =======================================================\n")

comp <- celltype_composition(obj[[]], "group", "celltype")
check("celltype_composition totals one per group",
      all(tapply(comp$proportion, comp$group, sum) == 1))

cat("\n== module score ======================================================\n")

scored <- add_module_score_auto(obj, paste0("Gene", 1:30), name = "TEST")
check("add_module_score_auto reports the number of genes",
      scored$n_genes >= 10)
check("add_module_score_auto adds a score column",
      scored$score_col %in% colnames(scored$object[[]]))

cat("\n== figures ===========================================================\n")

p1 <- plot_umap(obj, "celltype", pal_celltype_human,
                title = "test") +
  ggplot2::scale_colour_manual(
    values = pal_for(unique(obj$celltype), pal_celltype_human), guide = "none"
  )
check("plot_umap returns a ggplot", inherits(p1, "ggplot"))

out <- save_figure(p1, "it_umap", width = 3, height = 3)
check("save_figure writes a PDF", file.exists(p_fig("it_umap.pdf")))
check("save_figure writes a PNG", file.exists(p_fig("it_umap.png")))

p2 <- plot_umap_facets(obj, "celltype", "group", pal_celltype_human)
check("plot_umap_facets returns a ggplot", inherits(p2, "ggplot"))

p3 <- dotplot_ordered(obj, paste0("Gene", 1:12), "celltype",
                      order_by = paste0("Gene", 1:12))
check("dotplot_ordered returns a ggplot", inherits(p3, "ggplot"))

p4 <- plot_volcano(deg, label_genes = "Gene1", title = "test")
check("plot_volcano returns a ggplot", inherits(p4, "ggplot"))

p5 <- plot_composition(comp, pal_celltype_human,
                       levels = c("WT", "KOF"))
check("plot_composition returns a ggplot", inherits(p5, "ggplot"))

p6 <- plot_deg_tile(count_deg(deg_ct))
check("plot_deg_tile returns a ggplot", inherits(p6, "ggplot"))

top_markers <- data.frame(
  gene = paste0("Gene", c(1, 51, 101, 11, 61, 111, 21, 71, 121)),
  cluster = rep(c("A", "B", "C"), each = 3),
  avg_log2FC = rep(c(3, 2, 1), 3)
)
ht <- build_marker_heatmap(obj, top_markers, palette = pal_celltype_human)
check("build_marker_heatmap returns a Heatmap",
      methods::is(ht, "Heatmap"))
save_graphics("it_heatmap", 3, 3, ComplexHeatmap::draw(ht))
check("save_graphics writes a PDF", file.exists(p_fig("it_heatmap.pdf")))

cat("\n== tables and objects =================================================\n")

save_table(comp, "it_composition")
check("save_table writes a CSV", file.exists(p_tab("it_composition.csv")))

saved <- write_object(obj, "it_object")
check("write_object writes an RDS", file.exists(saved))
check("read_seurat reads it back",
      ncol(read_seurat(saved, "test object")) == n_cells)

missing <- tryCatch({
  read_seurat(file.path(scratch, "nope.rds"), "missing")
  FALSE
}, error = function(e) TRUE)
check("read_seurat fails with a helpful error", missing)

cat("\n== QC metrics ========================================================\n")

qc <- add_qc_metrics(obj, species = "mouse")
check("add_qc_metrics adds mt_percent and HB_percent",
      all(c("mt_percent", "HB_percent") %in% colnames(qc[[]])))

qc$group <- ordered_factor(qc$group, c("WT", "KOF"))
check("ordered_factor does not drop cells",
      !anyNA(qc$group) && length(qc$group) == n_cells)
check("ordered_factor sets the requested level order",
      identical(levels(qc$group), c("WT", "KOF")))

cat("\n", strrep("=", 70), "\n", sep = "")
cat(sprintf("%d checks, %d failure(s)\n", checks, failures))
cat("scratch directory:", scratch, "\n")
quit(save = "no", status = if (failures == 0) 0 else 1)
