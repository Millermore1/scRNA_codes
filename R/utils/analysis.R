# =============================================================================
# R/utils/analysis.R
#
# Reusable analysis steps: cluster relabelling, the standard Seurat workflow,
# differential expression across cell types, module scores, composition tables
# and GO enrichment helpers. Modules call these functions instead of repeating
# the same 20 lines of Seurat code.
# =============================================================================

# -----------------------------------------------------------------------------
# Seurat object handling
# -----------------------------------------------------------------------------

#' Join per-sample layers (Seurat v5) so that FindMarkers / GetAssayData work.
#'
#' The atlas objects are merged from many samples and therefore carry one layer
#' per sample. Differential expression needs a single counts layer, so this is
#' called before every DE step. No-op on Seurat v4 objects.
ensure_joined_layers <- function(obj, assay = "RNA") {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) return(obj)
  layers <- tryCatch(SeuratObject::Layers(obj[[assay]]), error = function(e) NULL)
  if (length(layers) > 1) {
    obj <- SeuratObject::JoinLayers(obj, assay = assay)
  }
  obj
}

#' Number of cells per identity level.
cells_per_ident <- function(obj, group.by = NULL) {
  if (is.null(group.by)) {
    table(Seurat::Idents(obj))
  } else {
    table(obj[[group.by, drop = TRUE]])
  }
}

# -----------------------------------------------------------------------------
# Cluster relabelling
# -----------------------------------------------------------------------------

#' Rename clusters using a manual mapping, with strict validation.
#'
#' The original scripts did `names(new.ids) <- levels(object)` and then
#' `RenameIdents()`. If the clustering produced a different number of clusters
#' the assignment silently recycled or errored deep inside Seurat; here the
#' mismatch is caught before anything is modified.
#'
#' @param obj Seurat object.
#' @param mapping Character vector of new labels, or a named vector
#'   `c("<cluster>" = "<label>")` for an explicit mapping.
#' @param new_col Name of the metadata column that receives the labels.
#' @param cluster_col Metadata column holding the current clusters. Defaults to
#'   the active identities.
#' @param drop_unmapped If TRUE, clusters missing from `mapping` are set to NA
#'   instead of raising an error.
relabel_clusters <- function(obj, mapping, new_col, cluster_col = NULL,
                             drop_unmapped = FALSE) {
  clusters <- if (is.null(cluster_col)) {
    as.character(Seurat::Idents(obj))
  } else {
    as.character(obj[[cluster_col, drop = TRUE]])
  }
  present <- sort(unique(clusters))

  if (!is.null(names(mapping)) && any(nzchar(names(mapping)))) {
    lookup <- mapping
  } else {
    if (length(mapping) != length(present)) {
      stop(
        "Cluster mapping has ", length(mapping), " labels but the object has ",
        length(present), " clusters (", paste(present, collapse = ", "), ").\n",
        "Update the mapping in config/config.R (dir: cfg$design) to match this ",
        "clustering, or change the clustering resolution.",
        call. = FALSE
      )
    }
    lookup <- stats::setNames(as.character(mapping), present)
  }

  labels <- unname(lookup[clusters])
  if (anyNA(labels) && !drop_unmapped) {
    missing <- setdiff(present, names(lookup))
    stop("No label defined for cluster(s): ", paste(missing, collapse = ", "),
         ". Add them to the mapping or set drop_unmapped = TRUE.",
         call. = FALSE)
  }

  # Keep the order of the mapping as the factor order, so that legends and
  # axes follow the order in which subtypes are listed in config/config.R.
  label_levels <- unique(unname(lookup))
  label_levels <- label_levels[!is.na(label_levels)]

  obj[[new_col]] <- factor(labels, levels = label_levels)
  Seurat::Idents(obj) <- new_col
  note("relabelled ", length(present), " clusters -> ",
       length(label_levels), " cell types")
  obj
}

# -----------------------------------------------------------------------------
# Standard workflow
# -----------------------------------------------------------------------------

#' Normalise, reduce and cluster a (sub)set of cells.
#'
#' Mirrors the workflow used throughout the manuscript: LogNormalize -> 2000-
#' 3000 variable features -> scale -> PCA -> optional Harmony batch correction
#' -> nearest neighbours -> Louvain clustering -> UMAP.
#'
#' @param obj Seurat object (a subset).
#' @param nfeatures,npcs,dims,resolution Parameters, see config/config.R.
#' @param batch_var Metadata column to correct for, or NULL to skip Harmony.
run_seurat_standard <- function(obj, nfeatures, npcs, dims, resolution,
                                batch_var = NULL, verbose = FALSE) {
  obj <- Seurat::NormalizeData(obj, verbose = verbose)
  obj <- Seurat::FindVariableFeatures(
    obj, selection.method = "vst", nfeatures = nfeatures, verbose = verbose
  )
  obj <- Seurat::ScaleData(obj, verbose = verbose)
  obj <- Seurat::RunPCA(obj, npcs = npcs, verbose = verbose)

  reduction <- "pca"
  if (!is.null(batch_var) && batch_var %in% colnames(obj[[]])) {
    if (!requireNamespace("harmony", quietly = TRUE)) {
      stop("Package 'harmony' is required for batch correction.", call. = FALSE)
    }
    obj <- harmony::RunHarmony(obj, group.by.vars = batch_var,
                               verbose = verbose)
    reduction <- "harmony"
  }

  obj <- Seurat::FindNeighbors(obj, reduction = reduction, dims = dims,
                               verbose = verbose)
  obj <- Seurat::FindClusters(obj, resolution = resolution, verbose = verbose)
  obj <- Seurat::RunUMAP(obj, reduction = reduction, dims = dims,
                         verbose = verbose)
  obj
}

#' Add a module score from a gene list, matching the gene-name style of the data.
#'
#' The public cohorts mix human (SLC17A7) and mouse (Slc17a7) symbol
#' conventions. This picks whichever spelling matches more genes in the object.
#'
#' `Seurat::AddModuleScore()` draws control genes from expression bins and
#' fails when a bin holds fewer genes than `ctrl.size` (default 100) -- which
#' happens on small objects. The binning and control size are therefore scaled
#' to the number of genes in the object.
#'
#' @return List with the updated object, the number of genes used and the name
#'   of the new metadata column.
add_module_score_auto <- function(obj, genes, name, ctrl_size = 100, nbin = 24) {
  all_genes <- rownames(obj)

  mouse_style <- paste0(
    toupper(substring(genes, 1, 1)),
    tolower(substring(genes, 2))
  )
  genes_mouse <- intersect(all_genes, mouse_style)
  genes_human <- intersect(all_genes, genes)

  use <- if (length(genes_mouse) >= length(genes_human)) genes_mouse else genes_human
  if (length(use) < 10) {
    warning("Only ", length(use), " genes of '", name,
            "' were found in the object -- the score may be unreliable.",
            call. = FALSE)
  }

  # Scale the binning to the size of the assay so AddModuleScore cannot fail.
  n_total <- length(all_genes)
  nbin_use <- max(2, min(nbin, floor(n_total / 100)))
  ctrl_use <- max(10, min(ctrl_size, floor(n_total / nbin_use / 3)))

  obj <- Seurat::AddModuleScore(
    obj, features = list(use), name = name,
    nbin = nbin_use, ctrl = ctrl_use
  )
  list(object = obj, n_genes = length(use), score_col = paste0(name, "1"))
}

# -----------------------------------------------------------------------------
# Differential expression
# -----------------------------------------------------------------------------

#' One DE test between two groups, optionally within a single cell type.
run_deg_pair <- function(obj, ident1, ident2, group_col = "group",
                         celltype = NULL, celltype_col = "celltype",
                         params = cfg$params$de) {
  obj <- ensure_joined_layers(obj)

  args <- list(
    object          = obj,
    ident.1         = ident1,
    ident.2         = ident2,
    group.by        = group_col,
    test.use        = params$test,
    min.pct         = params$min_pct,
    logfc.threshold = params$logfc_threshold,
    verbose         = FALSE
  )
  if (!is.null(celltype)) {
    Seurat::Idents(obj) <- celltype_col
    args$object      <- obj
    args$subset.ident <- celltype
  }

  res <- tryCatch(do.call(Seurat::FindMarkers, args),
                  error = function(e) {
                    note("skipped ", ident1, " vs ", ident2,
                         if (!is.null(celltype)) paste0(" in ", celltype),
                         ": ", conditionMessage(e))
                    NULL
                  })
  if (is.null(res) || nrow(res) == 0) return(NULL)

  data.frame(
    gene       = rownames(res),
    celltype   = if (is.null(celltype)) "all" else celltype,
    comparison = paste0(ident1, "_vs_", ident2),
    ident1     = ident1,
    ident2     = ident2,
    res,
    stringsAsFactors = FALSE,
    row.names  = NULL
  )
}

#' Run every contrast in every cell type.
#'
#' @param contrasts Data frame with columns `comparison`, `ident1`, `ident2`
#'   (see cfg$design$mouse_contrasts).
#' @return A single long data frame, one row per gene and cell type, or NULL.
run_deg_by_celltype <- function(obj,
                                contrasts = cfg$design$mouse_contrasts,
                                group_col = "group",
                                celltype_col = "celltype",
                                params = cfg$params$de) {
  obj <- ensure_joined_layers(obj)

  grp    <- obj[[group_col, drop = TRUE]]
  ct_vec <- obj[[celltype_col, drop = TRUE]]
  celltypes <- sort(unique(stats::na.omit(as.character(ct_vec))))
  out <- list()

  for (i in seq_len(nrow(contrasts))) {
    c1 <- contrasts$ident1[i]
    c2 <- contrasts$ident2[i]

    for (ct in celltypes) {
      sel <- !is.na(ct_vec) & as.character(ct_vec) == ct
      groups_present <- unique(as.character(grp[sel]))

      if (!all(c(c1, c2) %in% groups_present)) next

      n1 <- sum(sel & grp == c1, na.rm = TRUE)
      n2 <- sum(sel & grp == c2, na.rm = TRUE)
      if (n1 < params$min_cells || n2 < params$min_cells) {
        note("skipped ", ct, " (", c1, " n=", n1, ", ", c2, " n=", n2, ")")
        next
      }

      res <- run_deg_pair(obj, ident1 = c1, ident2 = c2,
                          group_col = group_col, celltype = ct,
                          celltype_col = celltype_col, params = params)
      if (!is.null(res)) out[[length(out) + 1]] <- res
    }
  }

  if (length(out) == 0) {
    warning("No contrast produced any DEGs.", call. = FALSE)
    return(NULL)
  }
  do.call(rbind, out)
}

#' Add a Up/Down/NS call and optionally filter on adjusted p value and FC.
filter_deg <- function(deg, params = cfg$params$de, significant_only = TRUE) {
  if (is.null(deg) || nrow(deg) == 0) return(deg)

  deg$direction <- ifelse(deg$avg_log2FC > 0, "Up", "Down")
  deg$significant <- deg$p_val_adj < params$padj_cutoff &
    abs(deg$avg_log2FC) > params$fc_cutoff

  if (!significant_only) return(deg)
  deg[deg$significant, , drop = FALSE]
}

#' Count significant DEGs per comparison and cell type (for the summary tile plot).
count_deg <- function(deg) {
  if (is.null(deg) || nrow(deg) == 0) {
    return(data.frame(comparison = character(0), celltype = character(0),
                      direction = character(0), n_deg = integer(0)))
  }
  sig <- filter_deg(deg)
  stats::aggregate(
    gene ~ comparison + celltype + direction,
    data = sig,
    FUN = length
  ) |> stats::setNames(c("comparison", "celltype", "direction", "n_deg"))
}

# -----------------------------------------------------------------------------
# Composition
# -----------------------------------------------------------------------------

#' Long table of cell-type proportions per group.
celltype_composition <- function(meta, group_col = "group",
                                 celltype_col = "celltype",
                                 group_levels = NULL) {
  df <- data.frame(group = meta[[group_col]],
                   celltype = meta[[celltype_col]]) |>
    stats::na.omit()
  df <- df[!is.na(df$group) & !is.na(df$celltype), , drop = FALSE]

  tab <- as.data.frame(table(group = df$group, celltype = df$celltype),
                       stringsAsFactors = FALSE)
  names(tab) <- c("group", "celltype", "count")
  tab <- tab[tab$count > 0, , drop = FALSE]

  totals <- stats::aggregate(count ~ group, data = tab, FUN = sum)
  names(totals)[2] <- "total"
  tab <- merge(tab, totals, by = "group")
  tab$proportion <- tab$count / tab$total

  if (!is.null(group_levels)) {
    tab$group <- ordered_factor(tab$group, group_levels)
  }
  tab[order(tab$group, tab$celltype), , drop = FALSE]
}

#' Proportions of one cluster across groups (used for the "cluster of interest"
#' stacked bar in the inhibitory analysis).
cluster_proportion <- function(meta, cluster_col, cluster_value,
                               group_col = "group") {
  df <- data.frame(group = meta[[group_col]],
                   cluster = as.character(meta[[cluster_col]]))
  df <- df[!is.na(df$group), , drop = FALSE]

  per_group <- table(df$group, df$cluster == cluster_value)
  out <- data.frame(
    group      = rownames(per_group),
    n_cluster  = as.integer(per_group[, "TRUE"]),
    n_total    = as.integer(rowSums(per_group)),
    stringsAsFactors = FALSE
  )
  out$proportion <- out$n_cluster / out$n_total
  out
}

# -----------------------------------------------------------------------------
# Enrichment
# -----------------------------------------------------------------------------

#' GO over-representation analysis with redundancy reduction.
#'
#' @param genes Gene symbols (any case).
#' @param orgdb "org.Mm.eg.db" or "org.Hs.eg.db".
#' @param ont One of "BP", "CC", "MF", "ALL".
#' @return data.frame of enriched terms, or NULL when nothing is significant.
run_enrich_go <- function(genes, orgdb = "org.Mm.eg.db", ont = "BP",
                          params = cfg$params$enrichment) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Package 'clusterProfiler' is required.", call. = FALSE)
  }
  if (!requireNamespace(orgdb, quietly = TRUE)) {
    stop("Annotation package '", orgdb, "' is required.", call. = FALSE)
  }

  genes <- unique(stats::na.omit(as.character(genes)))
  if (length(genes) < 5) {
    note("only ", length(genes), " genes -- skipping enrichment")
    return(NULL)
  }

  org_obj <- get(orgdb, envir = asNamespace(orgdb))

  eg <- tryCatch(
    clusterProfiler::enrichGO(
      gene          = genes,
      OrgDb         = org_obj,
      keyType       = "SYMBOL",
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = params$pvalue_cutoff,
      qvalueCutoff  = params$qvalue_cutoff,
      readable      = TRUE
    ),
    error = function(e) {
      note("enrichGO failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(eg) || nrow(as.data.frame(eg)) == 0) return(NULL)

  eg <- tryCatch(
    clusterProfiler::simplify(eg, cutoff = params$simplify_cutoff,
                              by = "p.adjust", select_fun = min),
    error = function(e) eg
  )
  as.data.frame(eg)
}

#' Turn a term table into nodes and Jaccard-weighted edges for a network figure.
#'
#' @param enrich_df Data frame from run_enrich_go().
#' @param top_n Number of most significant terms to keep.
build_enrich_graph_data <- function(enrich_df, top_n = 10,
                                    edge_jaccard = 0.05,
                                    min_shared = 1) {
  top <- enrich_df[order(enrich_df$p.adjust), , drop = FALSE]
  top <- utils::head(top, top_n)

  term_genes <- strsplit(top$geneID, "/")
  names(term_genes) <- top$ID

  edges <- do.call(rbind, lapply(utils::combn(names(term_genes), 2,
                                              simplify = FALSE),
    function(pr) {
      a <- term_genes[[pr[1]]]
      b <- term_genes[[pr[2]]]
      shared <- length(intersect(a, b))
      if (shared < min_shared) return(NULL)
      jac <- shared / length(union(a, b))
      if (is.na(jac) || jac < edge_jaccard) return(NULL)
      data.frame(from = pr[1], to = pr[2], weight = jac,
                 stringsAsFactors = FALSE)
    }))

  list(terms = top, term_genes = term_genes, edges = edges)
}

#' KEGG over-representation, mirroring run_enrich_go().
run_enrich_kegg <- function(genes, organism = "mmu",
                            params = cfg$params$enrichment) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Package 'clusterProfiler' is required.", call. = FALSE)
  }
  genes <- unique(stats::na.omit(as.character(genes)))
  if (length(genes) < 5) return(NULL)

  eg <- tryCatch(
    clusterProfiler::enrichKEGG(
      gene          = genes,
      organism      = organism,
      keyType       = "kegg",
      pAdjustMethod = "BH",
      pvalueCutoff  = params$pvalue_cutoff,
      qvalueCutoff  = params$qvalue_cutoff
    ),
    error = function(e) {
      note("enrichKEGG failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(eg) || nrow(as.data.frame(eg)) == 0) return(NULL)

  df <- as.data.frame(eg)
  df$Description <- sub(" - Mus musculus \\(house mouse\\)$", "",
                        df$Description)
  df
}

#' Convert an enrichment table into the two-line bar figure layout.
#'
#' Each term becomes one bar: the bar length is the fold enrichment
#' (GeneRatio / BgRatio) and the term name plus its member genes are drawn
#' inside the bar. This reproduces the layout of the published GO panels.
enrich_to_bar_df <- function(enrich_df,
                             top_n = cfg$params$enrichment$top_n_terms,
                             max_genes = cfg$params$enrichment$max_genes_in_label) {
  if (is.null(enrich_df) || nrow(enrich_df) == 0) return(NULL)

  df <- enrich_df[order(enrich_df$p.adjust), , drop = FALSE]
  df <- utils::head(df, top_n)

  gene_list <- strsplit(df$geneID, "/")
  df$gene_label <- vapply(gene_list, function(g) {
    if (length(g) > max_genes) {
      paste0(paste(g[seq_len(max_genes)], collapse = ", "), ", ...")
    } else {
      paste(g, collapse = ", ")
    }
  }, character(1))

  ratio    <- as.numeric(sub("/.*", "", df$GeneRatio)) /
    as.numeric(sub(".*/", "", df$GeneRatio))
  bg_ratio <- as.numeric(sub("/.*", "", df$BgRatio)) /
    as.numeric(sub(".*/", "", df$BgRatio))

  df$score <- ratio / bg_ratio
  df$Description <- stringr::str_to_sentence(df$Description)
  df$y      <- rev(seq_len(nrow(df)))
  df$term_y <- df$y + 0.22
  df$gene_y <- df$y - 0.22
  df
}

#' Truncate a comma-separated gene string for use as a figure label.
truncate_gene_label <- function(genes, max_genes = 8) {
  vapply(genes, function(g) {
    x <- trimws(strsplit(g, ",\\s*")[[1]])
    if (length(x) > max_genes) {
      paste0(paste(x[seq_len(max_genes)], collapse = ", "), ", ...")
    } else {
      paste(x, collapse = ", ")
    }
  }, character(1), USE.NAMES = FALSE)
}

# -----------------------------------------------------------------------------
# Misc
# -----------------------------------------------------------------------------

#' Normalise a numeric vector to [to[1], to[2]]; constant vectors map to the mean.
norm01 <- function(v, to = c(0.1, 0.9)) {
  r <- range(v, na.rm = TRUE)
  if (diff(r) == 0) return(rep(mean(to), length(v)))
  (v - r[1]) / diff(r) * diff(to) + to[1]
}

#' Safe wilcoxon test that returns NA instead of erroring on degenerate input.
safe_wilcox <- function(x, y) {
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  if (length(x) < 3 || length(y) < 3) return(NA_real_)
  tryCatch(stats::wilcox.test(x, y)$p.value, error = function(e) NA_real_)
}
