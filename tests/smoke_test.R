# =============================================================================
# tests/smoke_test.R
#
# Fast checks that run without any single-cell data:
#   1. every R file in the repository parses,
#   2. the configuration is internally consistent,
#   3. the colour palettes cover all categories used by the figures,
#   4. the shared helper functions behave correctly on toy input.
#
# Run with:  Rscript tests/smoke_test.R
# Exit code 0 means everything passed.
# =============================================================================

# Resolve the repository root from the location of this file.
local({
  args <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", args, value = TRUE)
  if (length(f) > 0) setwd(dirname(dirname(normalizePath(sub("^--file=", "", f[1])))))
})

failures <- 0L
checks <- 0L

check <- function(label, condition) {
  checks <<- checks + 1L
  ok <- isTRUE(condition)
  if (!ok) failures <<- failures + 1L
  cat(sprintf("%s %s\n", if (ok) "[ ok ]" else "[FAIL]", label))
  invisible(ok)
}

cat("== 1. syntax =========================================================\n")

r_files <- c(
  list.files("R", pattern = "\\.R$", recursive = TRUE, full.names = TRUE),
  list.files("config", pattern = "\\.R$", full.names = TRUE),
  "run_all.R"
)

for (f in r_files) {
  parsed <- tryCatch({
    parse(f)
    TRUE
  }, error = function(e) {
    cat("       ", conditionMessage(e), "\n")
    FALSE
  })
  check(paste("parses:", f), parsed)
}

cat("\n== 2. setup ==========================================================\n")

source("R/00_setup.R")

check("cfg is a list", is.list(cfg))
check("repository root resolved", dir.exists(cfg$root))
check("config has the expected sections",
      all(c("dirs", "inputs", "params", "design", "genes", "output") %in% names(cfg)))
check("contrasts are well formed",
      all(c("comparison", "ident1", "ident2") %in%
            colnames(cfg$design$mouse_contrasts)))
check("human contrasts have a control group",
      all(cfg$design$human_contrasts$ident2 == "NC"))
check("module files exist",
      all(file.exists(file.path("R", c(
        "00_setup.R", "01_prepare_datasets.R", "02_annotation.R",
        "03_cell_composition.R", "04_differential_expression.R",
        "05_inhibitory_neurons.R", "06_excitatory_neurons.R",
        "07_enrichment_networks.R", "08_cellchat.R", "09_wgcna.R",
        "10_figures_heatmaps.R", "11_asd_organoid.R"
      )))))

cat("\n== 3. palettes =======================================================\n")

check("human cell types all have a colour",
      !anyNA(pal_celltype_human[cfg$design$human_celltypes]))
check("mouse cell types all have a colour",
      !anyNA(pal_celltype_mouse[cfg$design$mouse_celltypes]))
check("mouse groups all have a colour",
      !anyNA(pal_group_mouse[cfg$design$mouse_group_order]))
check("human groups all have a colour",
      !anyNA(pal_group_human[cfg$design$human_group_order]))
check("interneuron subtypes all have a colour",
      !anyNA(pal_inhib_subtype[unique(cfg$design$mouse_inhibitory_mapping)]))
check("excitatory subtypes all have a colour",
      !anyNA(pal_excit_subtype[unique(cfg$design$mouse_excitatory_mapping)]))
check("neuron subtypes all have a colour",
      !anyNA(pal_neuron_subtype[unique(cfg$design$mouse_neuron_mapping)]))

unknown <- pal_for(c("Excitatory.neurons", "not-a-celltype"), pal_celltype_human)
check("pal_for falls back for unknown categories",
      length(unknown) == 2 && !anyNA(unknown))

check("ordered_factor puts known levels first",
      identical(levels(ordered_factor(c("B", "A", "Z"), c("A", "B"))),
                c("A", "B", "Z")))

cat("\n== 4. helper functions ===============================================\n")

v <- norm01(c(1, 2, 3, 4), to = c(0, 1))
check("norm01 maps the range onto [0, 1]",
      isTRUE(all.equal(range(v), c(0, 1))))
check("norm01 handles constant input",
      all(norm01(rep(5, 3)) == mean(c(0.1, 0.9))))

labels <- truncate_gene_label(c("A, B, C", "D"), max_genes = 2)
check("truncate_gene_label truncates long lists",
      identical(labels[1], "A, B, ...") && identical(labels[2], "D"))

enrich <- data.frame(
  ID = c("GO:1", "GO:2", "GO:3"),
  Description = c("synapse organization", "neuron projection", "gliogenesis"),
  p.adjust = c(0.001, 0.01, 0.02),
  geneID = c("A/B/C/D", "A/B/E", "F/G"),
  GeneRatio = c("4/10", "3/10", "2/10"),
  BgRatio = c("40/1000", "30/1000", "20/1000"),
  Count = c(4, 3, 2),
  stringsAsFactors = FALSE
)

g <- build_enrich_graph_data(enrich, top_n = 3, edge_jaccard = 0.1)
check("build_enrich_graph_data returns nodes and edges",
      is.list(g) && nrow(g$terms) == 3)
check("Jaccard edges link terms that share genes",
      nrow(g$edges) == 1 &&
        ((g$edges$from == "GO:1" && g$edges$to == "GO:2") ||
         (g$edges$from == "GO:2" && g$edges$to == "GO:1")))

bar <- enrich_to_bar_df(enrich, top_n = 3)
check("enrich_to_bar_df computes fold enrichment",
      isTRUE(all.equal(bar$score[bar$ID == "GO:1"], 10)))
check("enrich_to_bar_df keeps two-line label layout",
      all(c("term_y", "gene_y", "gene_label") %in% colnames(bar)))

deg <- data.frame(
  gene = c("G1", "G2", "G3", "G4"),
  celltype = c("A", "A", "A", "A"),
  comparison = "X_vs_Y",
  ident1 = "X", ident2 = "Y",
  p_val = c(0.001, 0.001, 0.5, 0.001),
  avg_log2FC = c(1.5, -1.2, 0.1, 0.05),
  pct.1 = 0.5, pct.2 = 0.2,
  p_val_adj = c(0.001, 0.001, 0.5, 0.001),
  stringsAsFactors = FALSE
)
sig <- filter_deg(deg)
check("filter_deg keeps only significant genes", nrow(sig) == 2)
check("filter_deg labels direction",
      identical(sort(sig$direction), c("Down", "Up")))

counts <- count_deg(deg)
check("count_deg summarises one row per direction",
      nrow(counts) == 2 && sum(counts$n_deg) == 2)

meta <- data.frame(
  group = c("X", "X", "Y", "Y"),
  celltype = c("A", "B", "A", "B"),
  orig.ident = c("s1", "s1", "s2", "s2"),
  stringsAsFactors = FALSE
)
comp <- celltype_composition(meta, "group", "celltype")
check("celltype_composition proportions sum to one per group",
      isTRUE(all.equal(
        tapply(comp$proportion, comp$group, sum)[["X"]], 1
      )))

check("safe_wilcox returns NA for tiny groups",
      is.na(safe_wilcox(c(1, 2), c(3, 4, 5, 6))))

cat("\n== 5. metadata =======================================================\n")

check("cohort_mapping.csv exists", file.exists("data/metadata/cohort_mapping.csv"))
check("hp_autism_genes.txt exists", file.exists("data/metadata/hp_autism_genes.txt"))

if (file.exists("data/metadata/cohort_mapping.csv")) {
  cm <- utils::read.csv("data/metadata/cohort_mapping.csv",
                        stringsAsFactors = FALSE)
  check("cohort_mapping has the expected columns",
        all(c("sample", "group") %in% colnames(cm)))
  check("cohort_mapping groups are known",
        all(cm$group %in% cfg$design$human_group_order))
}

if (file.exists("data/metadata/hp_autism_genes.txt")) {
  hp <- readLines("data/metadata/hp_autism_genes.txt")
  hp <- hp[nzchar(trimws(hp))]
  check("HP_AUTISM gene list has more than 100 genes", length(hp) > 100)
  check("HP_AUTISM gene list is unique", !anyDuplicated(hp))
}

cat("\n", strrep("=", 70), "\n", sep = "")
cat(sprintf("%d checks, %d failure(s)\n", checks, failures))
quit(save = "no", status = if (failures == 0) 0 else 1)
