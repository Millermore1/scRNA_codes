# =============================================================================
# env/install_packages.R
#
# Install every R package the pipeline needs.
#
#   Rscript env/install_packages.R
#
# The versions in env/package_versions.tsv are the ones the analysis was run
# with (R 4.3.0 / Bioconductor 3.18). Newer versions normally work, but Seurat
# and CellChat have both changed APIs across major releases, so pin them if you
# need to reproduce the original figures exactly.
# =============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = 600)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_packages <- c(
  "Seurat", "SeuratObject", "harmony",
  "tidyverse", "dplyr", "tidyr", "stringr", "tibble", "purrr",
  "ggplot2", "patchwork", "scales", "viridis", "ggh4x",
  "RColorBrewer", "pheatmap", "circlize",
  "igraph", "ggraph", "tidygraph",
  "scatterpie", "ggforce", "ggnewscale", "ggrepel", "ggtext",
  "ggVennDiagram", "openxlsx", "writexl", "readxl",
  "Matrix", "remotes", "devtools"
)

bioc_packages <- c(
  "ComplexHeatmap",
  "clusterProfiler", "org.Mm.eg.db", "org.Hs.eg.db", "AnnotationDbi",
  "limma", "EnhancedVolcano",
  "WGCNA", "hdWGCNA"
)

# Packages that only live on GitHub.
github_packages <- c(
  "chris-mcginnis-ucsf/DoubletFinder",
  "sqjin/CellChat"
)

install_missing <- function(pkgs, installer = install.packages,
                            label = "CRAN") {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) == 0) {
    cat("[", label, "] everything already installed\n", sep = "")
    return(invisible(NULL))
  }
  cat("[", label, "] installing: ", paste(missing, collapse = ", "), "\n",
      sep = "")
  installer(missing)
}

cat("\n== CRAN ==============================================================\n")
install_missing(cran_packages, install.packages, "CRAN")

cat("\n== Bioconductor ======================================================\n")
install_missing(bioc_packages, BiocManager::install, "Bioc")

cat("\n== GitHub ============================================================\n")
for (repo in github_packages) {
  pkg <- sub(".*/", "", repo)
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat("[GitHub] ", pkg, " already installed\n", sep = "")
    next
  }
  cat("[GitHub] installing ", repo, "\n", sep = "")
  tryCatch(
    remotes::install_github(repo, upgrade = "never", build_vignettes = FALSE),
    error = function(e) {
      message("Could not install ", repo, ": ", conditionMessage(e))
      message("CellChat and DoubletFinder often need the GitHub token unset: ",
              "Sys.unsetenv(c('GITHUB_TOKEN', 'GITHUB_PAT'))")
    }
  )
}

cat("\n== done ==============================================================\n")
installed <- c(cran_packages, bioc_packages, sub(".*/", "", github_packages))
status <- data.frame(
  package = installed,
  installed = vapply(installed, requireNamespace, logical(1), quietly = TRUE),
  row.names = NULL
)
print(status[!status$installed, , drop = FALSE], row.names = FALSE)
