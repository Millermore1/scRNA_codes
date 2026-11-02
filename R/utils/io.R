# =============================================================================
# R/utils/io.R
#
# Filesystem helpers: every output path, every figure and every table written
# by the pipeline goes through these functions so that (i) directories are
# created on demand and (ii) file names stay consistent.
# =============================================================================

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------

#' Build a path inside `data/processed/`.
p_processed <- function(...) file.path(cfg$dirs$processed, ...)

#' Build a path inside `results/figures/`.
p_fig <- function(...) file.path(cfg$dirs$figures, ...)

#' Build a path inside `results/tables/`.
p_tab <- function(...) file.path(cfg$dirs$tables, ...)

#' Build a path inside `results/objects/`.
p_obj <- function(...) file.path(cfg$dirs$objects, ...)

#' Build a path inside `data/metadata/`.
p_meta <- function(...) file.path(cfg$dirs$metadata, ...)

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

# -----------------------------------------------------------------------------
# Logging
# -----------------------------------------------------------------------------

#' Print a section header so long pipeline runs stay readable.
step_header <- function(text) {
  bar <- strrep("-", max(0, 74 - nchar(text)))
  message("\n", text, " ", bar)
  invisible(text)
}

#' Informative message that is silent when the script runs non-interactively.
note <- function(...) message("  * ", paste0(...))

# -----------------------------------------------------------------------------
# Reading
# -----------------------------------------------------------------------------

#' Read a Seurat object, failing with an actionable message when it is missing.
read_seurat <- function(path, what = "Seurat object") {
  if (!file.exists(path)) {
    stop(
      what, " not found:\n  ", path, "\n",
      "Run the module that creates it first, or point cfg$inputs at your copy ",
      "(see config/config.R).", call. = FALSE
    )
  }
  obj <- readRDS(path)
  if (!inherits(obj, "Seurat")) {
    stop(path, " does not contain a Seurat object (class: ",
         paste(class(obj), collapse = "/"), ").", call. = FALSE)
  }
  obj
}

#' Read any RDS, failing with an actionable message when it is missing.
read_rds_required <- function(path, what = "input object") {
  if (!file.exists(path)) {
    stop(what, " not found:\n  ", path, call. = FALSE)
  }
  readRDS(path)
}

#' Read a CSV that the pipeline itself produced in an earlier step.
read_table_required <- function(path, what = "input table") {
  if (!file.exists(path)) {
    stop(what, " not found:\n  ", path, call. = FALSE)
  }
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

#' Resolve an input path from cfg$inputs, checking that it exists.
input_file <- function(name) {
  path <- cfg$inputs[[name]]
  if (is.null(path)) {
    stop("cfg$inputs$", name, " is not defined.", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop("Input '", name, "' not found:\n  ", path,
         "\nSee docs/data_availability.md for how to obtain it.", call. = FALSE)
  }
  path
}

# -----------------------------------------------------------------------------
# Writing
# -----------------------------------------------------------------------------

#' Save a Seurat object (or any object) into results/objects/.
write_object <- function(object, name, dir = cfg$dirs$objects) {
  ensure_dir(dir)
  path <- file.path(dir, paste0(name, ".rds"))
  saveRDS(object, path)
  note("saved object: ", basename(path))
  invisible(path)
}

#' Save a ggplot / patchwork figure as PDF (vector, for submission) + PNG.
#'
#' @param plot A ggplot, patchwork or grid object.
#' @param name File stem, e.g. "Figure2_celltype_umap".
#' @param width,height Size in inches.
save_figure <- function(plot, name, width = 6, height = 5,
                        dpi = cfg$output$figure_dpi) {
  ensure_dir(cfg$dirs$figures)

  pdf_path <- p_fig(paste0(name, ".pdf"))
  grDevices::cairo_pdf(pdf_path, width = width, height = height, onefile = TRUE)
  print(plot)
  grDevices::dev.off()

  if (isTRUE(cfg$output$save_png)) {
    png_path <- p_fig(paste0(name, ".png"))
    ggplot2::ggsave(
      filename = png_path, plot = plot,
      width = width, height = height, dpi = dpi, bg = "white"
    )
  }

  note("saved figure: ", basename(pdf_path))
  invisible(pdf_path)
}

#' Save base-graphics output (CellChat circle plots, pheatmap, ...) to PDF+PNG.
#'
#' `code` is an expression that draws onto the open device, e.g.
#'   save_graphics("fig", 8, 6, { pheatmap(mat) })
save_graphics <- function(name, width, height, code,
                          dpi = cfg$output$figure_dpi) {
  ensure_dir(cfg$dirs$figures)
  expr <- substitute(code)
  env  <- parent.frame()

  pdf_path <- p_fig(paste0(name, ".pdf"))
  grDevices::cairo_pdf(pdf_path, width = width, height = height)
  eval(expr, envir = env)
  grDevices::dev.off()

  if (isTRUE(cfg$output$save_png)) {
    png_path <- p_fig(paste0(name, ".png"))
    grDevices::png(png_path, width = width, height = height,
                   units = "in", res = dpi, bg = "white")
    eval(expr, envir = env)
    grDevices::dev.off()
  }

  note("saved figure: ", basename(pdf_path))
  invisible(pdf_path)
}

#' Write a data frame to results/tables/ as CSV.
save_table <- function(df, name, digits = NULL) {
  ensure_dir(cfg$dirs$tables)
  if (!is.null(digits)) {
    num <- vapply(df, is.numeric, logical(1))
    df[num] <- lapply(df[num], round, digits = digits)
  }
  path <- p_tab(paste0(name, ".csv"))
  utils::write.csv(df, path, row.names = FALSE, fileEncoding = "UTF-8")
  note("saved table: ", basename(path), " (", nrow(df), " rows)")
  invisible(path)
}

#' Persist the exact software environment used for a run.
save_session_info <- function(name = "session_info") {
  ensure_dir(cfg$dirs$results)
  path <- file.path(cfg$dirs$results, paste0(name, ".txt"))
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(paste("R version:", R.version.string), con)
  writeLines(paste("Platform :", R.version$platform), con)
  writeLines(paste("Run at   :", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), con)
  writeLines("", con)
  writeLines("Attached / loaded packages:", con)
  writeLines(utils::capture.output(utils::sessionInfo()), con)
  invisible(path)
}
