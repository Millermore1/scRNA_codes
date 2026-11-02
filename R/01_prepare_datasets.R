# =============================================================================
# R/01_prepare_datasets.R
#
# Build the two Seurat objects that every downstream module reads:
#
#   human_atlas_integrated.rds   multi-cohort human brain snRNA-seq
#                                (Alzheimer's disease, Parkinson's disease,
#                                 multiple sclerosis, epilepsy, controls)
#   mouse_ko_integrated.rds      Prmt8 conditional-KO snRNA-seq with
#                                chemogenetic (DREADD) manipulation
#
# Steps: read -> QC -> merge -> normalise / HVG / PCA -> Harmony -> cluster.
# The GSE203201 organoid cohort is prepared separately in 11_asd_organoid.R.
#
# Outputs
#   data/processed/human_atlas_integrated.rds
#   data/processed/mouse_ko_integrated.rds
#   results/tables/01_qc_summary_<dataset>.csv
#   results/figures/01_qc_violin_<dataset>.pdf
# =============================================================================

if (!file.exists("R/00_setup.R")) stop("Run this script from the repository root.")
source("R/00_setup.R")

need("Seurat", "dplyr", "ggplot2", "patchwork", "harmony", "Matrix")

# -----------------------------------------------------------------------------
# 1. Readers
# -----------------------------------------------------------------------------

#' Read one Cell Ranger output directory into a Seurat object.
read_cellranger_dir <- function(dir, project = basename(dir),
                                min_cells = 3, min_features = 200) {
  counts <- Seurat::Read10X(data.dir = dir)
  Seurat::CreateSeuratObject(
    counts = counts, project = project,
    min.cells = min_cells, min.features = min_features
  )
}

#' Read every Cell Ranger output directory inside `parent_dir`.
read_cellranger_dirs <- function(parent_dir, pattern = NULL,
                                 min_cells = 3, min_features = 200) {
  if (!dir.exists(parent_dir)) {
    stop("Directory not found: ", parent_dir, call. = FALSE)
  }
  samples <- list.dirs(parent_dir, recursive = FALSE, full.names = FALSE)
  if (!is.null(pattern)) samples <- grep(pattern, samples, value = TRUE)
  if (length(samples) == 0) {
    stop("No sample directories inside ", parent_dir, call. = FALSE)
  }
  note("found ", length(samples), " sample directories")

  objs <- lapply(samples, function(s) {
    read_cellranger_dir(file.path(parent_dir, s), project = s,
                        min_cells = min_cells, min_features = min_features)
  })
  names(objs) <- samples
  objs
}

#' Read a Cell Ranger .h5 matrix.
read_h5_matrix <- function(h5_file, project = NULL,
                           min_cells = 3, min_features = 200) {
  if (!file.exists(h5_file)) stop("File not found: ", h5_file, call. = FALSE)

  if (is.null(project)) {
    project <- tools::file_path_sans_ext(basename(h5_file))
  }

  counts <- Seurat::Read10X_h5(h5_file)
  Seurat::CreateSeuratObject(
    counts = counts, project = project,
    min.cells = min_cells, min.features = min_features
  )
}

#' Read a dense/sparse expression matrix supplied as a delimited text file.
#'
#' Used for cohorts distributed outside Cell Ranger (e.g. the MS snRNA-seq
#' expression matrix).
read_custom_matrix <- function(file, sep = "\t", row_names = 1,
                               min_cells = 3, min_features = 200) {
  if (!file.exists(file)) stop("File not found: ", file, call. = FALSE)

  mat <- utils::read.delim(file, header = TRUE, row.names = row_names,
                           check.names = FALSE)
  mat <- as.matrix(mat)
  mat <- methods::as(Matrix::Matrix(mat, sparse = TRUE), "dgCMatrix")

  Seurat::CreateSeuratObject(
    counts = mat,
    project = tools::file_path_sans_ext(basename(file)),
    min.cells = min_cells, min.features = min_features
  )
}

#' Read a matrix distributed as MatrixMarket + barcodes/features text files.
#'
#' Handles the non-standard layout used by several GEO submissions where the
#' three files carry a study prefix instead of the Cell Ranger names.
read_mtx_triplet <- function(dir, matrix_file = "Human_counts.mtx",
                             genes_file = "Human_genes.txt",
                             meta_file  = "Human_cell_meta.txt") {
  if (!dir.exists(dir)) stop("Directory not found: ", dir, call. = FALSE)

  counts <- Matrix::readMM(file.path(dir, matrix_file))
  genes  <- readLines(file.path(dir, genes_file))
  meta   <- utils::read.table(file.path(dir, meta_file), header = TRUE,
                              sep = "\t", stringsAsFactors = FALSE)

  rownames(counts) <- genes
  # Prefer the study's own cell names over the barcode column.
  cell_names <- if ("cell_name" %in% colnames(meta)) meta$cell_name else meta$cell_id
  colnames(counts) <- cell_names

  counts <- methods::as(counts, "CsparseMatrix")  # keeps memory low
  Seurat::CreateSeuratObject(counts = counts, meta.data = meta)
}

# -----------------------------------------------------------------------------
# 2. Quality control
#
# add_qc_metrics(), apply_qc_filter(), qc_summary_table() and plot_qc_violin()
# live in R/utils/qc.R, because module 11 (organoids) uses them too.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# 3. Integration
# -----------------------------------------------------------------------------

#' Normalise, reduce and cluster a merged atlas.
integrate_atlas <- function(obj, params = cfg$params$integration,
                            batch_var = params$batch_var,
                            run_tsne = FALSE) {
  obj <- ensure_joined_layers(obj)

  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj <- Seurat::FindVariableFeatures(
    obj, selection.method = "vst", nfeatures = params$nfeatures, verbose = FALSE
  )
  obj <- Seurat::ScaleData(obj, verbose = FALSE)
  obj <- Seurat::RunPCA(obj, npcs = params$npcs, verbose = FALSE)

  obj <- harmony::RunHarmony(obj, group.by.vars = batch_var, verbose = FALSE)

  obj <- Seurat::FindNeighbors(obj, reduction = "harmony", dims = params$dims,
                               verbose = FALSE)
  obj <- Seurat::FindClusters(obj, resolution = params$resolution,
                              verbose = FALSE)
  obj <- Seurat::RunUMAP(obj, reduction = "harmony", dims = params$dims,
                         verbose = FALSE)
  if (run_tsne) {
    obj <- Seurat::RunTSNE(obj, reduction = "harmony", dims = params$dims)
  }
  obj
}

#' Score the cell cycle (human gene symbols).
add_cell_cycle <- function(obj) {
  s_genes   <- Seurat::CaseMatch(search = Seurat::cc.genes.updated.2019$s.genes,
                                 match = rownames(obj))
  g2m_genes <- Seurat::CaseMatch(search = Seurat::cc.genes.updated.2019$g2m.genes,
                                 match = rownames(obj))
  if (length(s_genes) < 5 || length(g2m_genes) < 5) {
    note("cell-cycle genes not found -- skipping CellCycleScoring")
    return(obj)
  }
  Seurat::CellCycleScoring(obj, s.features = s_genes, g2m.features = g2m_genes)
}

# -----------------------------------------------------------------------------
# 4. Human multi-cohort atlas
# -----------------------------------------------------------------------------

#' Assign cohort labels (AD / PD / MS / EPI / NC) from sample identifiers.
#'
#' data/metadata/cohort_mapping.csv lists every sample name and its cohort. Any
#' sample that is not listed can still be mapped by the GSE prefix contained in
#' its identifier, which is what the original script did with grepl().
assign_cohort_groups <- function(obj, mapping_file = cfg$inputs$cohort_mapping,
                                 id_col = "orig.ident") {
  ids <- as.character(obj[[id_col, drop = TRUE]])

  if (file.exists(mapping_file)) {
    map <- utils::read.csv(mapping_file, stringsAsFactors = FALSE)
    stopifnot(all(c("sample", "group") %in% colnames(map)))
    lookup <- stats::setNames(map$group, map$sample)
    groups <- unname(lookup[ids])
    note("mapped ", sum(!is.na(groups)), "/", length(ids),
         " cells using cohort_mapping.csv")
  } else {
    groups <- rep(NA_character_, length(ids))
    note("cohort_mapping.csv not found -- falling back to identifier patterns")
  }

  # Fallback: derive the cohort from the sample name pattern.
  fill <- is.na(groups)
  if (any(fill)) {
    groups[fill] <- ifelse(grepl("AD", ids[fill]), "AD",
                    ifelse(grepl("PD", ids[fill]), "PD",
                    ifelse(grepl("MS", ids[fill]), "MS",
                    ifelse(grepl("EPI|Epi", ids[fill]), "EPI",
                    ifelse(grepl("NC|CO[0-9]", ids[fill]), "NC", NA)))))
  }

  obj$group <- groups
  obj
}

#' Full build of the human atlas from Cell Ranger directories.
build_human_atlas <- function(input_dir = cfg$inputs$human_10x_dir,
                              output = cfg$inputs$human_atlas_rds) {
  step_header("Human multi-cohort atlas")

  if (!dir.exists(input_dir)) {
    stop(
      "Human 10x directory not found: ", input_dir, "\n",
      "Set PRMT8_HUMAN_10X_DIR or edit cfg$inputs$human_10x_dir.",
      call. = FALSE
    )
  }

  objs <- read_cellranger_dirs(input_dir)
  objs <- lapply(objs, add_qc_metrics, species = "human")
  objs <- lapply(objs, apply_qc_filter, species = "human")

  merged <- if (length(objs) == 1) objs[[1]] else merge(objs[[1]], objs[-1])
  merged <- ensure_joined_layers(merged)
  merged <- add_cell_cycle(merged)
  merged <- assign_cohort_groups(merged)

  save_table(qc_summary_table(merged, "group"), "01_qc_summary_human")

  qc_plot <- plot_qc_violin(merged, group_col = "group",
                            palette = pal_group_human,
                            group_levels = cfg$design$human_group_order)
  save_figure(qc_plot, "01_qc_violin_human", width = 8, height = 6)

  step_header("Harmony integration")
  atlas <- integrate_atlas(merged, run_tsne = TRUE)
  write_object(atlas, basename(tools::file_path_sans_ext(output)),
               dir = dirname(output))
  atlas
}

# -----------------------------------------------------------------------------
# 5. Mouse Prmt8 cKO atlas
# -----------------------------------------------------------------------------

#' Constrain the mouse atlas groups to the experimental design and order them.
prepare_mouse_atlas <- function(path = cfg$inputs$mouse_atlas_rds) {
  step_header("Mouse Prmt8 cKO atlas")

  obj <- read_seurat(path, "mouse atlas")
  obj <- ensure_joined_layers(obj)

  if (!"group" %in% colnames(obj[[]])) {
    stop("The mouse object has no 'group' column in its metadata.", call. = FALSE)
  }

  unknown <- setdiff(unique(as.character(obj$group)),
                     cfg$design$mouse_group_order)
  if (length(unknown) > 0) {
    warning("Groups not in cfg$design$mouse_group_order: ",
            paste(unknown, collapse = ", "), call. = FALSE)
  }
  obj$group <- ordered_factor(obj$group, cfg$design$mouse_group_order)

  if (!"mt_percent" %in% colnames(obj[[]])) {
    obj <- add_qc_metrics(obj, species = "mouse")
  }

  # The mouse atlas used a stricter mitochondrial threshold (5%).
  keep <- obj$mt_percent < cfg$params$qc$mouse_max_mt_percent
  if (any(!keep, na.rm = TRUE)) {
    note("removing ", sum(!keep, na.rm = TRUE), " cells with mt > ",
         cfg$params$qc$mouse_max_mt_percent, "%")
    obj <- subset(obj, cells = colnames(obj)[which(keep)])
  }

  save_table(qc_summary_table(obj, "group"), "01_qc_summary_mouse")
  save_figure(
    plot_qc_violin(obj, "group", pal_group_mouse,
                   cfg$design$mouse_group_order),
    "01_qc_violin_mouse", width = 8, height = 6
  )

  write_object(obj, "mouse_ko_integrated")
  obj
}

# -----------------------------------------------------------------------------
# 6. Run
#
# Sourcing this file executes the module. run_all.R sources the modules in
# order; see README.md.
# -----------------------------------------------------------------------------

human <- if (file.exists(cfg$inputs$human_atlas_rds)) {
  note("human atlas already integrated -- loading ", cfg$inputs$human_atlas_rds)
  read_seurat(cfg$inputs$human_atlas_rds, "human atlas")
} else {
  build_human_atlas()
}

mouse <- prepare_mouse_atlas()

note("human atlas : ", ncol(human), " cells")
note("mouse atlas : ", ncol(mouse), " cells, groups: ",
     paste(levels(mouse$group), collapse = ", "))
