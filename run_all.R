# =============================================================================
# run_all.R
#
# Run the whole pipeline, or a subset of modules, from a single command:
#
#   Rscript run_all.R              # every module, in order
#   Rscript run_all.R 01 04 07     # only modules whose file name starts with
#                                  # one of these prefixes
#   Rscript run_all.R --list       # show the modules and what they do
#
# Each module is a self-contained script that reads the objects written by the
# previous one, so modules can also be run individually (see README.md).
# =============================================================================

# Resolve the repository root even when called from another directory.
local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script <- normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
    setwd(dirname(script))
  } else if (file.exists("run_all.R")) {
    setwd(getwd())
  } else {
    stop("Run this script from the repository root, or call it with Rscript.",
         call. = FALSE)
  }
})

modules <- data.frame(
  file = c("01_prepare_datasets.R",
           "02_annotation.R",
           "03_cell_composition.R",
           "04_differential_expression.R",
           "05_inhibitory_neurons.R",
           "06_excitatory_neurons.R",
           "07_enrichment_networks.R",
           "08_cellchat.R",
           "09_wgcna.R",
           "10_figures_heatmaps.R",
           "11_asd_organoid.R"),
  description = c(
    "Read, QC and integrate the human and mouse atlases",
    "Annotate cell types, marker plots, UMAP panels",
    "Cell-type composition and sample-level proportion tests",
    "Differential expression across cell types and contrasts",
    "Inhibitory subtypes, PRMT8 expression, module score",
    "Excitatory layer subtypes",
    "GO/KEGG enrichment, GO-with-genes bar plots, term networks",
    "CellChat communication analysis and comparisons",
    "hdWGCNA co-expression modules and the PRMT8 module",
    "Expression heat maps, PRMT8 correlation panel",
    "GSE203201 ASD organoid re-analysis"
  ),
  stringsAsFactors = FALSE
)

args <- commandArgs(trailingOnly = TRUE)

if (any(args %in% c("--list", "-l"))) {
  cat("Modules:\n")
  for (i in seq_len(nrow(modules))) {
    cat(sprintf("  %-32s %s\n", modules$file[i], modules$description[i]))
  }
  quit(save = "no", status = 0)
}

select <- if (length(args) == 0 || all(grepl("^--", args))) {
  rep(TRUE, nrow(modules))
} else {
  prefixes <- args[!grepl("^--", args)]
  Reduce(`|`, lapply(prefixes, function(p) startsWith(modules$file, p)))
}

stop_on_error <- !("--keep-going" %in% args)

cat(strrep("=", 78), "\n")
cat("PRMT8 scRNA-seq pipeline\n")
cat("repository:", getwd(), "\n")
cat("R         :", R.version.string, "\n")
cat("modules   :", sum(select), "of", nrow(modules), "\n")
cat(strrep("=", 78), "\n\n")

results <- data.frame(module = character(0), status = character(0),
                      seconds = numeric(0), stringsAsFactors = FALSE)

for (i in which(select)) {
  file <- modules$file[i]
  cat("\n", strrep("#", 78), "\n", sep = "")
  cat("## ", file, " -- ", modules$description[i], "\n", sep = "")
  cat(strrep("#", 78), "\n", sep = "")

  start <- Sys.time()
  status <- tryCatch({
    source(file.path("R", file), echo = FALSE)
    "ok"
  }, error = function(e) {
    message("MODULE FAILED: ", conditionMessage(e))
    if (stop_on_error) stop(e) else "failed"
  })
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))

  results <- rbind(results, data.frame(module = file, status = status,
                                       seconds = round(elapsed, 1)))
  cat(sprintf("-- %s finished in %.1f s (%s)\n", file, elapsed, status))
}

cat("\n", strrep("=", 78), "\n", sep = "")
print(results, row.names = FALSE)
cat(strrep("=", 78), "\n")

if (any(results$status == "failed")) quit(save = "no", status = 1)
