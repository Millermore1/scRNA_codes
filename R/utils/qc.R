# =============================================================================
# R/utils/qc.R
#
# Quality-control helpers shared by the atlas module (01) and the organoid
# module (11): metric calculation, filtering, summary tables and the QC panel.
# =============================================================================

#' Add mitochondrial and haemoglobin fractions, using the species' gene naming.
#'
#' Returns 0 for a metric when no matching gene is present, so that a
#' species/naming mismatch degrades into a warning instead of an error.
add_qc_metrics <- function(obj, species = c("human", "mouse")) {
  species <- match.arg(species)

  mt_pattern <- if (species == "human") "^MT-" else "^mt-"
  hb_genes <- if (species == "human") {
    c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM",
      "HBQ1", "HBZ")
  } else {
    c("Hba-a1", "Hba-a2", "Hbb-b1", "Hbb-b2", "Hbb-bh1", "Hbb-y")
  }

  mt_genes <- grep(mt_pattern, rownames(obj), value = TRUE)
  if (length(mt_genes) == 0) {
    warning("No mitochondrial genes matched '", mt_pattern,
            "' -- mt_percent set to 0. Check the gene naming of this object.",
            call. = FALSE)
    obj[["mt_percent"]] <- 0
  } else {
    obj[["mt_percent"]] <- Seurat::PercentageFeatureSet(obj, features = mt_genes)
  }

  hb_present <- intersect(hb_genes, rownames(obj))
  if (length(hb_present) == 0) {
    obj[["HB_percent"]] <- 0
  } else {
    obj[["HB_percent"]] <- Seurat::PercentageFeatureSet(obj,
                                                        features = hb_present)
  }
  obj
}

#' Apply the cell-level QC thresholds from config/config.R.
#'
#' The upper UMI bound is a per-sample quantile, so filtering must be applied to
#' each sample before merging (as in the original workflow).
apply_qc_filter <- function(obj, params = cfg$params$qc,
                            species = c("human", "mouse")) {
  species <- match.arg(species)
  max_mt <- if (species == "mouse") {
    params$mouse_max_mt_percent
  } else {
    params$max_mt_percent
  }

  n_before <- ncol(obj)
  counts   <- obj$nCount_RNA

  keep <- obj$nFeature_RNA > params$min_features &
    obj$nFeature_RNA < params$max_features &
    obj$mt_percent   < max_mt &
    counts > params$min_counts &
    counts < stats::quantile(counts, params$upper_count_quantile)

  obj <- subset(obj, cells = colnames(obj)[keep])
  note("QC kept ", ncol(obj), "/", n_before, " cells")
  obj
}

#' Table of QC metrics per sample or per group, for the supplementary table.
qc_summary_table <- function(obj, group_col = NULL) {
  meta <- obj[[]]
  key <- if (!is.null(group_col) && group_col %in% colnames(meta)) {
    group_col
  } else {
    "orig.ident"
  }

  groups <- split(meta, meta[[key]])
  out <- do.call(rbind, lapply(names(groups), function(g) {
    m <- groups[[g]]
    data.frame(
      sample        = g,
      n_cells       = nrow(m),
      median_genes  = stats::median(m$nFeature_RNA),
      median_umi    = stats::median(m$nCount_RNA),
      median_mt_pct = round(stats::median(m$mt_percent), 2),
      stringsAsFactors = FALSE
    )
  }))
  out[order(out$sample), , drop = FALSE]
}

#' Violin panels for the four QC metrics, split by group.
plot_qc_violin <- function(obj, group_col = "group",
                           palette = pal_group_mouse,
                           group_levels = NULL) {
  metrics <- c("nFeature_RNA", "nCount_RNA", "mt_percent", "HB_percent")
  labels  <- c("Detected genes", "UMI counts", "MT (%)", "HB (%)")
  keep    <- metrics %in% colnames(obj[[]])
  metrics <- metrics[keep]
  labels  <- labels[keep]

  df <- obj[[]]
  if (!is.null(group_levels)) {
    df[[group_col]] <- ordered_factor(df[[group_col]], group_levels)
  }

  plots <- lapply(seq_along(metrics), function(i) {
    m <- metrics[i]
    ggplot2::ggplot(df, ggplot2::aes(x = .data[[group_col]],
                                     y = .data[[m]],
                                     fill = .data[[group_col]])) +
      ggplot2::geom_violin(scale = "width", trim = TRUE, linewidth = 0.2) +
      ggplot2::geom_boxplot(width = 0.15, fill = "white", outlier.size = 0,
                            linewidth = 0.25) +
      ggplot2::scale_fill_manual(values = pal_for(
        unique(as.character(df[[group_col]])), palette)) +
      ggplot2::scale_y_continuous(labels = scales::comma) +
      ggplot2::labs(x = NULL, y = labels[i]) +
      theme_pub(base_size = 9) +
      ggplot2::theme(
        legend.position = "none",
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )
  })

  patchwork::wrap_plots(plots, nrow = 2) +
    patchwork::plot_annotation(tag_levels = "a") &
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 12))
}
