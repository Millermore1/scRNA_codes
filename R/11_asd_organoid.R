# =============================================================================
# R/11_asd_organoid.R
#
# Re-analysis of GSE203201, a human cortical organoid snRNA-seq series from
# ASD patients and controls sampled at developmental time points
# (D0, D6, D12, D30, D60, D100). This dataset is used to ask whether PRMT8
# expression differs between ASD and control organoids during development.
#
# Steps: file normalisation -> read -> QC -> Harmony -> clustering ->
#        ASD vs control DE at D30 -> PRMT8 expression trajectory ->
#        PRMT8-high/low DE and GSEA.
#
# Outputs
#   data/processed/asd_organoid.rds
#   results/tables/11_*.csv
#   results/figures/11_*.pdf
# =============================================================================

if (!file.exists("R/00_setup.R")) stop("Run this script from the repository root.")
source("R/00_setup.R")

need("Seurat", "dplyr", "ggplot2", "patchwork", "harmony", "limma",
     "clusterProfiler", "org.Hs.eg.db", "pheatmap")

# -----------------------------------------------------------------------------
# 1. Normalise the GEO file names
# -----------------------------------------------------------------------------

#' Rename study-prefixed 10x files to the Cell Ranger convention.
#'
#' GEO submissions frequently look like
#'   ASD_1.barcodes.tsv.gz  ASD_1.features.tsv.gz  ASD_1.matrix.mtx.gz
#' which Read10X() cannot read. This creates one folder per sample containing
#' the three files under their canonical names.
normalise_10x_filenames <- function(dir, dry_run = FALSE) {
  files <- list.files(dir, pattern = "\\.gz$", full.names = FALSE)
  if (length(files) == 0) return(invisible(0))

  samples <- unique(sub("[\\._](barcodes|genes|features|matrix).*", "", files))
  note("detected ", length(samples), " samples: ",
       paste(utils::head(samples, 5), collapse = ", "),
       if (length(samples) > 5) ", ..." else "")

  moved <- 0L
  for (s in samples) {
    target_dir <- file.path(dir, s)
    if (!dry_run) dir.create(target_dir, showWarnings = FALSE)

    s_files <- list.files(dir, pattern = paste0("^", s, "[\\._]"),
                          full.names = FALSE)
    for (f in s_files) {
      new_name <- if (grepl("barcodes", f)) {
        "barcodes.tsv.gz"
      } else if (grepl("matrix", f)) {
        "matrix.mtx.gz"
      } else if (grepl("features|genes", f)) {
        "features.tsv.gz"
      } else {
        NA_character_
      }
      if (is.na(new_name)) next

      from <- file.path(dir, f)
      to   <- file.path(target_dir, new_name)
      if (!dry_run && !file.exists(to)) {
        file.rename(from, to)
        moved <- moved + 1L
      }
    }
  }
  note("moved ", moved, " files into per-sample folders")
  invisible(moved)
}

# -----------------------------------------------------------------------------
# 2. Read the cohort
# -----------------------------------------------------------------------------

#' Read every sample folder and derive group and time point metadata.
read_organoid_cohort <- function(dir) {
  samples <- list.dirs(dir, recursive = FALSE, full.names = FALSE)
  samples <- samples[samples != ""]
  if (length(samples) == 0) {
    stop("No sample folders in ", dir, " -- run normalise_10x_filenames() first.",
         call. = FALSE)
  }

  objs <- lapply(samples, function(s) {
    obj <- read_cellranger_dir(file.path(dir, s), project = s,
                               min_cells = 3, min_features = 200)

    obj$Group <- if (grepl("ASD", s)) {
      "ASD"
    } else if (grepl("IMR|NC", s)) {
      "Control"
    } else {
      "Other"
    }

    tp <- regmatches(s, regexpr("D[0-9]+", s))
    obj$Timepoint <- if (length(tp) > 0) tp else NA_character_

    obj
  })
  names(objs) <- samples

  if (length(objs) == 1) objs[[1]] else merge(objs[[1]], objs[-1])
}

# -----------------------------------------------------------------------------
# 3. Build the processed object
# -----------------------------------------------------------------------------

build_organoid_object <- function(input_dir = cfg$inputs$asd_10x_dir,
                                  output = p_processed("asd_organoid.rds")) {
  step_header("GSE203201 ASD organoids")

  if (!dir.exists(input_dir)) {
    stop("GSE203201 directory not found: ", input_dir, call. = FALSE)
  }
  if (length(list.dirs(input_dir, recursive = FALSE)) == 0) {
    normalise_10x_filenames(input_dir)
  }

  obj <- read_organoid_cohort(input_dir)
  note("merged object: ", ncol(obj), " cells from ", length(unique(obj$orig.ident)),
       " samples")

  obj <- add_qc_metrics(obj, species = "human")

  par <- cfg$params$organoid
  keep <- obj$nFeature_RNA > par$min_features &
    obj$nFeature_RNA < par$max_features &
    obj$mt_percent < par$max_mt
  note("QC kept ", sum(keep), "/", ncol(obj), " cells")
  obj <- subset(obj, cells = colnames(obj)[keep])

  obj$Timepoint <- factor(obj$Timepoint, levels = par$timepoints)
  save_table(qc_summary_table(obj, "Group"), "11_qc_summary_organoid")

  obj <- ensure_joined_layers(obj)
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj <- Seurat::FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
  obj <- Seurat::ScaleData(obj, verbose = FALSE)
  obj <- Seurat::RunPCA(obj, npcs = 30, verbose = FALSE)
  obj <- harmony::RunHarmony(obj, group.by.vars = "orig.ident", verbose = FALSE)
  obj <- Seurat::FindNeighbors(obj, reduction = "harmony", dims = 1:30,
                               verbose = FALSE)
  obj <- Seurat::FindClusters(obj, resolution = par$resolution, verbose = FALSE)
  obj <- Seurat::RunUMAP(obj, reduction = "harmony", dims = 1:30,
                         verbose = FALSE)

  save_figure(
    plot_umap(obj, "Timepoint",
              palette = setNames(recycle_colours(length(par$timepoints)),
                                 par$timepoints),
              title = "GSE203201 organoids by time point"),
    "11_umap_timepoint_organoid", width = 4.8, height = 4
  )

  save_figure(
    plot_umap(obj, "Group", pal_group_organoid,
              title = "ASD vs control"),
    "11_umap_group_organoid", width = 4.4, height = 4
  )

  write_object(obj, "asd_organoid")
  obj
}

# -----------------------------------------------------------------------------
# 4. Run
# -----------------------------------------------------------------------------

obj <- if (file.exists(p_processed("asd_organoid.rds"))) {
  read_seurat(p_processed("asd_organoid.rds"), "ASD organoid object")
} else {
  build_organoid_object()
}

# -----------------------------------------------------------------------------
# 5. Differential expression at the focus time point
# -----------------------------------------------------------------------------

step_header("ASD vs control at ", cfg$params$organoid$focus_timepoint)

focus <- cfg$params$organoid$focus_timepoint
obj_focus <- subset(obj, cells = colnames(obj)[
  as.character(obj$Timepoint) == focus & as.character(obj$Group) != "Other"
])

if (ncol(obj_focus) < 100) {
  warning("Only ", ncol(obj_focus), " cells at ", focus, " -- skipping DE.",
          call. = FALSE)
} else {
  save_table(qc_summary_table(obj_focus, "Group"),
             paste0("11_qc_summary_", tolower(focus)))

  deg <- run_deg_pair(obj_focus, ident1 = "ASD", ident2 = "Control",
                      group_col = "Group")
  if (!is.null(deg)) {
    save_table(deg, paste0("11_deg_ASD_vs_Control_", focus), digits = 4)
    save_table(filter_deg(deg),
               paste0("11_deg_ASD_vs_Control_", focus, "_significant"),
               digits = 4)

    save_figure(
      plot_volcano(deg, label_genes = c("PRMT8"),
                   title = paste0("ASD vs control (", focus, ")")),
      paste0("11_volcano_ASD_vs_Control_", focus), width = 4.8, height = 4.2
    )

    go <- run_enrich_go(unique(filter_deg(deg)$gene[filter_deg(deg)$direction == "Up"]),
                        "org.Hs.eg.db", "BP")
    if (!is.null(go)) {
      save_table(go, paste0("11_go_up_ASD_vs_Control_", focus), digits = 4)
    }
  }
}

# -----------------------------------------------------------------------------
# 6. PRMT8 expression across time points
# -----------------------------------------------------------------------------

step_header("PRMT8 trajectory")

if ("PRMT8" %in% rownames(obj)) {
  prmt8 <- Seurat::FetchData(obj, vars = c("PRMT8", "Timepoint", "Group",
                                           "orig.ident"))
  names(prmt8)[1] <- "expr"
  prmt8 <- prmt8[!is.na(prmt8$Timepoint) & prmt8$Group %in% c("ASD", "Control"), ]

  trend <- stats::aggregate(expr ~ Timepoint + Group, data = prmt8,
                            FUN = function(x) c(mean = mean(x), sd = stats::sd(x),
                                                n = length(x)))
  trend <- data.frame(
    Timepoint = trend$Timepoint,
    Group     = trend$Group,
    mean_expr = trend$expr[, "mean"],
    sd_expr   = trend$expr[, "sd"],
    n         = trend$expr[, "n"]
  )
  trend$se_expr <- trend$sd_expr / sqrt(trend$n)
  save_table(trend, "11_prmt8_trajectory_organoid", digits = 4)

  save_figure(
    ggplot2::ggplot(trend, ggplot2::aes(x = .data$Timepoint, y = .data$mean_expr,
                                        colour = .data$Group,
                                        group = .data$Group)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_point(size = 2.6) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$mean_expr - .data$se_expr,
                     ymax = .data$mean_expr + .data$se_expr),
        width = 0.15, linewidth = 0.5
      ) +
      ggplot2::scale_colour_manual(values = pal_group_organoid,
                                   name = "Group") +
      ggplot2::labs(x = "Time point", y = "PRMT8 mean expression ± SEM") +
      theme_pub(),
    "11_prmt8_trajectory_organoid", width = 4.6, height = 3.6
  )

  save_figure(
    plot_feature(obj, "PRMT8", limits = c("q5", "q95"), point_size = 0.05,
                 title = "PRMT8 in organoids"),
    "11_prmt8_featureplot_organoid", width = 4.6, height = 4
  )
}

# -----------------------------------------------------------------------------
# 7. PRMT8-high vs PRMT8-low cells (limma + GSEA)
# -----------------------------------------------------------------------------

step_header("PRMT8-high vs PRMT8-low")

if ("PRMT8" %in% rownames(obj)) {
  obj <- ensure_joined_layers(obj)
  expr <- Seurat::FetchData(obj, vars = "PRMT8")[, 1]
  obj$PRMT8_group <- ifelse(expr > stats::median(expr), "High", "Low")

  mat <- as.matrix(Seurat::GetAssayData(obj, assay = "RNA", slot = "data"))
  design <- stats::model.matrix(~ obj$PRMT8_group)
  fit <- limma::lmFit(mat, design)
  fit <- limma::eBayes(fit)
  tt <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "P")

  save_table(data.frame(gene = rownames(tt), tt, row.names = NULL),
             "11_limma_prmt8_high_vs_low_organoid", digits = 4)

  ranks <- stats::setNames(tt$logFC, rownames(tt))
  ranks <- sort(ranks[!is.na(ranks)], decreasing = TRUE)

  gsea <- tryCatch(
    clusterProfiler::gseGO(
      geneList     = ranks,
      OrgDb        = org.Hs.eg.db,
      ont          = "BP",
      keyType      = "SYMBOL",
      minGSSize    = 15,
      maxGSSize    = 500,
      pvalueCutoff = 0.05,
      verbose      = FALSE
    ),
    error = function(e) {
      note("gseGO failed: ", conditionMessage(e))
      NULL
    }
  )

  if (!is.null(gsea) && nrow(as.data.frame(gsea)) > 0) {
    gsea_df <- as.data.frame(gsea)
    save_table(gsea_df, "11_gsea_prmt8_high_vs_low_organoid", digits = 4)

    save_figure(clusterProfiler::dotplot(gsea, showCategory = 20) +
                  ggplot2::ggtitle("GSEA: PRMT8-high vs PRMT8-low"),
                "11_gsea_dotplot_organoid", width = 6, height = 5)
  }
}

note("done: organoid module complete")
