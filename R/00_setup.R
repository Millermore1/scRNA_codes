# =============================================================================
# R/00_setup.R
#
# Entry point for every pipeline module. It
#   1. locates the repository root (so scripts can be run from anywhere),
#   2. loads the project configuration and colour palettes,
#   3. sources the utility functions,
#   4. attaches the R packages required by a module, with a readable error
#      message when one is missing,
#   5. creates the output directories.
#
# Every other script starts with:
#     if (!file.exists("R/00_setup.R")) stop("Run from the repository root")
#     source("R/00_setup.R")
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Repository root
# -----------------------------------------------------------------------------

find_repo_root <- function(start = getwd(), max_up = 8) {
  path <- normalizePath(start, winslash = "/", mustWork = FALSE)
  for (i in seq_len(max_up)) {
    if (file.exists(file.path(path, "config", "config.R"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) break
    path <- parent
  }
  stop(
    "Could not locate the repository root (no config/config.R found in any ",
    "parent directory of ", start, ").\n",
    "Set it explicitly with Sys.setenv(PRMT8_REPO_ROOT = '<path to repo>').",
    call. = FALSE
  )
}

.root <- Sys.getenv("PRMT8_REPO_ROOT", unset = "")
if (!nzchar(.root)) .root <- find_repo_root()
options(prmt8.root = normalizePath(.root, winslash = "/", mustWork = TRUE))

# -----------------------------------------------------------------------------
# 2. Configuration and utilities
# -----------------------------------------------------------------------------

source(file.path(getOption("prmt8.root"), "config", "config.R"))
source(file.path(getOption("prmt8.root"), "config", "palettes.R"))

local({
  utils_dir <- file.path(getOption("prmt8.root"), "R", "utils")
  for (f in c("io.R", "qc.R", "analysis.R", "plotting.R")) {
    source(file.path(utils_dir, f))
  }
})

# -----------------------------------------------------------------------------
# 3. Package loading
# -----------------------------------------------------------------------------

# Versions this project was developed and tested with. Used only to produce a
# helpful warning; the exact environment is recorded by save_session_info().
.tested_versions <- c(
  Seurat = "5.3.0", SeuratObject = "5.1.0", harmony = "1.2.0",
  tidyverse = "2.0.0", dplyr = "1.1.4", ggplot2 = "3.4.0",
  patchwork = "1.1.3", ComplexHeatmap = "2.18.0", circlize = "0.4.17",
  pheatmap = "1.0.12", CellChat = "1.6.1", clusterProfiler = "4.10.0",
  "org.Mm.eg.db" = "3.18.0", "org.Hs.eg.db" = "3.15.0",
  igraph = "2.0.1", ggraph = "2.1.0", tidygraph = "1.3.0",
  scatterpie = "0.2.1", ggforce = "0.4.1", ggnewscale = "0.4.9",
  ggrepel = "0.9.3", ggtext = "0.1.2", scales = "1.4.0",
  openxlsx = "4.2.8.1", writexl = "1.5.4", readxl = "1.4.3",
  DoubletFinder = "2.0.6", hdWGCNA = "0.4.12", WGCNA = "1.74",
  ggVennDiagram = "1.5.7", EnhancedVolcano = "1.20.0", limma = "3.58.1",
  ggpubr = "0.6.1", Matrix = "1.6.4", ggh4x = "0.3.1", viridis = "0.6.4"
)

#' Attach packages, failing with installation instructions if one is missing.
#'
#' @param ... Package names.
#' @param quietly If TRUE, suppress the package start-up messages.
need <- function(..., quietly = TRUE) {
  pkgs <- c(...)
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing) > 0) {
    stop(
      "Missing required package(s): ", paste(missing, collapse = ", "), "\n",
      "Install them with:\n  source('env/install_packages.R')\n",
      "or individually, e.g.:\n  install.packages('", missing[1], "')\n",
      "Bioconductor packages: BiocManager::install('", missing[1], "')",
      call. = FALSE
    )
  }

  for (p in pkgs) {
    suppressPackageStartupMessages(
      library(p, character.only = TRUE, quietly = quietly)
    )
    tested <- .tested_versions[[p]]
    if (!is.null(tested)) {
      installed <- as.character(utils::packageVersion(p))
      if (!identical(installed, tested)) {
        message("  ! ", p, " ", installed, " (this project was tested with ",
                tested, ")")
      }
    }
  }
  invisible(pkgs)
}

# -----------------------------------------------------------------------------
# 4. Output directories and session banner
# -----------------------------------------------------------------------------

invisible(lapply(cfg$dirs[c("processed", "figures", "tables", "objects")],
                 ensure_dir))

set.seed(cfg$params$seed)

fancy_banner <- function() {
  message(
    "\n", strrep("=", 74),
    "\n PRMT8 / interneuron scRNA-seq analysis pipeline",
    "\n repository : ", cfg$root,
    "\n R          : ", R.version.string,
    "\n started    : ", format(Sys.time(), "%Y-%m-%d %H:%M"),
    "\n", strrep("=", 74)
  )
}

if (interactive()) fancy_banner()

# Convenience aliases used throughout the modules.
path_raw <- function(...) file.path(cfg$dirs$raw, ...)
