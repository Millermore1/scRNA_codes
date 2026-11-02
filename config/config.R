# =============================================================================
# config/config.R
#
# Central configuration for the PRMT8 scRNA-seq analysis pipeline.
# Everything that a reader might need to change lives here: input locations,
# analysis parameters and contrast definitions. The analysis scripts read
# from `cfg` and never hard-code a path or a threshold.
#
# Paths are resolved relative to the repository root, which is detected by
# R/00_setup.R. Absolute input locations can be overridden without editing
# this file by exporting environment variables (see `cfg$inputs` below).
# =============================================================================

.root <- getOption("prmt8.root")
if (is.null(.root)) {
  stop("Option 'prmt8.root' is not set. Source R/00_setup.R instead of ",
       "sourcing config/config.R directly.", call. = FALSE)
}

# -----------------------------------------------------------------------------
# 1. Directory layout
# -----------------------------------------------------------------------------
cfg <- list()

cfg$root <- .root

cfg$dirs <- list(
  config    = file.path(.root, "config"),
  code      = file.path(.root, "R"),
  data      = file.path(.root, "data"),
  raw       = file.path(.root, "data", "raw"),
  metadata  = file.path(.root, "data", "metadata"),
  processed = file.path(.root, "data", "processed"),
  results   = file.path(.root, "results"),
  figures   = file.path(.root, "results", "figures"),
  tables    = file.path(.root, "results", "tables"),
  objects   = file.path(.root, "results", "objects"),
  docs      = file.path(.root, "docs")
)

# -----------------------------------------------------------------------------
# 2. Input files
#
# Two ways to supply the data:
#   (a) put the files in data/raw/ and keep the defaults below, or
#   (b) export environment variables, e.g. in ~/.Renviron
#         PRMT8_INPUT_DIR = "D:/scrnadata"
#       and every dataset is looked up below that directory.
#
# The three primary datasets are:
#   human_atlas   multi-cohort human brain snRNA-seq  (AD, PD, MS, EPI, NC)
#   mouse_ko      Prmt8 cKO / DREADD snRNA-seq        (WT, KOF, KOF_*_CNO)
#   asd_organoid  GSE203201 human cortical organoids  (ASD, IMR, NC; D0-D100)
# -----------------------------------------------------------------------------

.input_dir <- Sys.getenv("PRMT8_INPUT_DIR", unset = cfg$dirs$raw)

cfg$inputs <- list(

  # -- human multi-cohort atlas -----------------------------------------------
  # Directory holding one Cell Ranger output folder per sample, or a single
  # .rds with the already-integrated Seurat object.
  human_10x_dir    = Sys.getenv("PRMT8_HUMAN_10X_DIR",
                                unset = file.path(.input_dir, "human_10x")),
  human_atlas_rds  = Sys.getenv("PRMT8_HUMAN_ATLAS_RDS",
                                unset = file.path(cfg$dirs$processed,
                                                  "human_atlas_integrated.rds")),
  # Non-standard matrices used for the MS / PD / EPI cohorts.
  human_ms_matrix  = Sys.getenv("PRMT8_MS_MATRIX",
                                unset = file.path(.input_dir,
                                                  "MSCtr_snRNA_ExpressionMatrix_R.txt")),
  gse186537_dir    = Sys.getenv("PRMT8_GSE186537_DIR",
                                unset = file.path(.input_dir, "GSE186537")),
  gse174367_h5     = Sys.getenv("PRMT8_GSE174367_H5",
                                unset = file.path(.input_dir,
                                                  "GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5")),

  # -- mouse Prmt8 cKO / chemogenetic atlas ----------------------------------
  mouse_atlas_rds  = Sys.getenv("PRMT8_MOUSE_ATLAS_RDS",
                                unset = file.path(cfg$dirs$processed,
                                                  "mouse_ko_integrated.rds")),

  # -- GSE203201 ASD organoids -----------------------------------------------
  asd_10x_dir      = Sys.getenv("PRMT8_ASD_10X_DIR",
                                unset = file.path(.input_dir, "GSE203201")),
  asd_h5_dir       = Sys.getenv("PRMT8_ASD_H5_DIR",
                                unset = file.path(.input_dir, "GSE203201_h5")),

  # -- annotation resources --------------------------------------------------
  # Downloaded from MSigDB / published papers; see docs/data_availability.md.
  hp_autism_gmt    = Sys.getenv("PRMT8_HP_AUTISM_GMT",
                                unset = file.path(.input_dir, "genesets",
                                                  "hp_autism.gmt")),
  brain_cor_xlsx   = Sys.getenv("PRMT8_BRAIN_COR_XLSX",
                                unset = file.path(.input_dir,
                                                  "correlations_PRMT8_brain_region.xlsx")),
  cohort_mapping   = file.path(cfg$dirs$metadata, "cohort_mapping.csv"),
  celltype_markers = file.path(cfg$dirs$metadata, "celltype_markers.csv")
)

# -----------------------------------------------------------------------------
# 3. Re-analysis parameters
#    These values reproduce the settings used for the manuscript. Each block is
#    referenced by exactly one pipeline module, named in the comment.
# -----------------------------------------------------------------------------

cfg$params <- list(

  # -- 01_prepare_datasets.R --------------------------------------------------
  qc = list(
    min_features          = 300,    # nFeature_RNA lower bound
    max_features          = 5000,   # nFeature_RNA upper bound
    min_counts            = 1000,   # nCount_RNA lower bound
    upper_count_quantile  = 0.97,   # drop cells above this per-sample UMI quantile
    max_mt_percent        = 10,     # mitochondrial reads (%)
    max_hb_percent        = 3,      # haemoglobin reads (%)  (human samples)
    mouse_max_mt_percent  = 5       # stricter mt filter applied to the mouse atlas
  ),

  integration = list(
    nfeatures    = 3000,
    npcs         = 30,
    batch_var    = "orig.ident",
    dims         = 1:10,
    resolution   = 0.05
  ),

  # -- 02_annotation.R -------------------------------------------------------
  annotation = list(
    # Marker panel shown in the cell-type dot plot (human).
    human_markers = c("CLDN5", "PECAM1", "ABCB1",
                      "P2RY12", "DOCK8", "CSF1R",
                      "SLC1A3", "GFAP", "AQP4",
                      "GAD1", "GAD2", "PVALB",
                      "MOG", "MBP", "OLIG1",
                      "GRIN1", "CAMK2A", "CBLN2"),
    mouse_markers = c("Slc17a7", "Slc17a6", "Camk2a",
                      "Gad1", "Gad2", "Slc6a1",
                      "Mog", "Mbp", "Plp1",
                      "Pdgfra", "Cspg4", "Olig2",
                      "Aqp4", "Gfap", "Aldh1l1",
                      "P2ry12", "Cx3cr1", "Tmem119",
                      "Pecam1", "Cldn5", "Kdr"),
    marker_heatmap_top_n = 8,       # top markers per cell type in the heat map
    marker_heatmap_min_fc = 1       # |avg_log2FC| threshold for marker selection
  ),

  # -- 05_inhibitory_neurons.R ----------------------------------------------
  inhibitory = list(
    nfeatures  = 2000,
    npcs       = 20,
    dims       = 1:20,
    resolution = 0.1
  ),

  neuron_subclustering = list(
    nfeatures  = 3000,
    npcs       = 30,
    dims       = 1:30,
    resolution = 0.03
  ),

  # -- 06_excitatory_neurons.R ----------------------------------------------
  excitatory = list(
    nfeatures  = 2000,
    npcs       = 20,
    dims       = 1:20,
    resolution = 0.1,
    remove_cluster = "3"            # cluster excluded before re-clustering
  ),

  # -- 04_differential_expression.R -----------------------------------------
  de = list(
    test            = "wilcox",
    min_pct         = 0.1,
    logfc_threshold = 0.25,
    padj_cutoff     = 0.05,
    fc_cutoff       = 0.25,
    min_cells       = 5             # cell types with fewer cells are skipped
  ),

  # -- 07_enrichment_networks.R ---------------------------------------------
  enrichment = list(
    pvalue_cutoff   = 0.05,
    qvalue_cutoff   = 0.2,
    simplify_cutoff = 0.6,
    top_n_terms     = 10,
    edge_jaccard    = 0.05,
    min_shared_genes = 1,
    max_genes_in_label = 8
  ),

  # -- 08_cellchat.R --------------------------------------------------------
  cellchat = list(
    db             = "mouse",
    min_cells      = 10,
    type           = "triMean",
    population_size = TRUE,
    # Pathways highlighted in the synapse-focused figure.
    synapse_pathways = c("GLUTAMATE", "GABA", "NRXN", "NLGN", "NCAM",
                         "CADM", "EPHA", "EPHB", "NRG", "PTN", "MK")
  ),

  # -- 09_wgcna.R -----------------------------------------------------------
  wgcna = list(
    gene_select  = "fraction",
    fraction     = 0.05,
    k            = 25,
    max_shared   = 10,
    network_type = "signed",
    soft_power   = 7,
    n_hubs       = 100,
    trait        = "PD_Status"
  ),

  # -- 11_asd_organoid.R ----------------------------------------------------
  organoid = list(
    min_features = 200,
    max_features = 6000,
    max_mt       = 10,
    dims         = 1:20,
    resolution   = 0.1,
    timepoints   = c("D0", "D6", "D12", "D30", "D60", "D100"),
    focus_timepoint = "D30"
  ),

  seed = 42
)

# -----------------------------------------------------------------------------
# 4. Experimental design and contrasts
#    `ident1` is the numerator of the log2 fold change.
# -----------------------------------------------------------------------------

cfg$design <- list(

  human_group_order  = c("NC", "AD", "PD", "MS", "EPI"),
  mouse_group_order  = c("WT", "KOF", "KOF_vector_CNO", "KOF_Gq_CNO"),
  organoid_group_order = c("Control", "ASD"),

  human_celltypes = c("Excitatory.neurons", "Inhibitory.neurons",
                      "Oligodendrocytes", "Astrocytes", "Microglia",
                      "Endothelial.cells"),

  mouse_celltypes = c("Excitatory neurons", "Inhibitory neurons",
                      "Oligodendrocyte", "OPC", "Astrocyte", "Microglia",
                      "Endothelial cells"),

  # Manual cluster -> cell type mappings.
  # Lengths must equal the number of clusters; the pipeline checks this and
  # stops with an informative message if a run produces a different number.
  human_cluster_mapping = c(
    "Excitatory.neurons", "Oligodendrocytes", "Inhibitory.neurons",
    "Astrocytes", "Inhibitory.neurons", "Microglia", "Inhibitory.neurons",
    "Endothelial.cells"
  ),

  # Neuron compartment (module 08): 15 clusters at resolution 0.05.
  mouse_neuron_mapping = c("PN", "PVALB", "PN", "PN", "SST", "PVALB", "PN",
                           "VIP", "SST", "PN", "PN", "PN", "PN", "SST", "PN"),

  mouse_inhibitory_mapping = c("PVALB", "SST", "PVALB", "VIP", "PVALB",
                               "PVALB", "VIP", "PVALB", "VIP", "SST", "SST"),

  mouse_excitatory_mapping = c("Exc_L2/3", "Exc_L6", "Exc_L4", "Exc_L4",
                               "Exc_L5", "Exc_L2/3", "Exc_L5", "Exc_L5",
                               "Exc_L5", "Exc_L5", "Exc_L6", "Exc_L2/3",
                               "Exc_L2/3"),

  # Differentially expressed gene contrasts (mouse atlas).
  mouse_contrasts = data.frame(
    comparison = c("KOF_vs_WT", "KOF_Gq_CNO_vs_KOF_vector_CNO"),
    ident1     = c("KOF", "KOF_Gq_CNO"),
    ident2     = c("WT", "KOF_vector_CNO"),
    stringsAsFactors = FALSE
  ),

  # Disease contrasts for the human atlas.
  human_contrasts = data.frame(
    comparison = c("AD_vs_NC", "PD_vs_NC", "MS_vs_NC", "EPI_vs_NC"),
    ident1     = c("AD", "PD", "MS", "EPI"),
    ident2     = c("NC", "NC", "NC", "NC"),
    stringsAsFactors = FALSE
  ),

  # ASD organoid contrast.
  organoid_contrasts = data.frame(
    comparison = "ASD_vs_Control",
    ident1     = "ASD",
    ident2     = "Control",
    stringsAsFactors = FALSE
  )
)

# -----------------------------------------------------------------------------
# 5. Cell-type-specific gene panels used across figures
# -----------------------------------------------------------------------------

cfg$genes <- list(
  prmt8        = "PRMT8",
  prmt8_mouse  = "Prmt8",

  # Interneuron subtype markers (mouse, from the inhibitory subclustering).
  inhib_subtypes = list(
    VIP   = c("Vip", "Crh", "Igf1", "Reln", "Vipr2"),
    SST   = c("Sst", "Gpr149", "Thsd7a", "Chodl"),
    PV    = c("Pvalb", "Gabra1", "Gabra4", "Kcnc1"),
    TH    = c("Th", "Slc6a3", "Ddc", "Slc18a2"),
    LAMP5 = c("Lamp5", "Reln", "Lhx1", "Lhx2", "Nos1"),
    PAX6  = c("Pax6", "Vip", "Reln", "Igf1")
  ),

  # Human interneuron panel.
  human_inhib = c("PRMT8", "PVALB", "SST", "VIP", "CR", "NPY", "GAD1", "GAD2"),

  # Mouse interneuron panel (title-case symbols).
  mouse_inhib = c("Prmt8", "Pvalb", "Sst", "Vip", "Calb2", "Npy", "Gad1",
                  "Gad2"),

  # Cortical layer markers for the excitatory subclustering (mouse).
  excit_layers = list(
    L2_3 = c("Cux1", "Cux2", "Reln"),
    L4   = c("Rorb", "Rprm"),
    L5   = c("Bcl11b", "Fezf2", "Tbr1"),
    L6   = c("Tle4", "Foxp2", "Pcp4")
  ),

  # Excitatory neuron ordering panel (dot plot sorted by summed expression).
  excit_ordering = c("Satb2", "Slc17a7", "Bcl11b", "Tle4", "Foxp2",
                     "Rxfp1", "Rorb", "Cux2", "Prmt8"),

  # Ion-channel / receptor panels for the expression heat maps (human).
  ion_channels = c("CACNA1A", "CACNA1B", "CACNA1E",
                   "HCN1", "HCN2", "HCN3", "HCN4",
                   "SCN1A", "SCN1B", "SCN8A",
                   "VAMP2", "SLC32A1"),

  gaba_receptors = c("GABRA1", "GABRA2", "GABRB2", "GABRB3", "GABRD", "GABRG2",
                     "GABBR1", "GABBR2",
                     "GRIN1", "GRIN2A", "GRIN2B", "GRIN3A", "GRIN3B",
                     "GRIA1", "GRIA2", "GRIA3", "GRIA4",
                     "GRIP1", "NTRK2"),

  # ASD / synaptic risk genes highlighted in the volcano plot.
  asd_label_genes = c("Auts2", "Tcf4", "Ctnnb1", "Chd7", "Zeb2", "Pcdh10",
                      "Cdh10", "Dscam", "Dscaml1", "Kirrel3", "Ncam2", "Ntng1",
                      "Cacna1c", "Cacna1d", "Scn3a", "Gabrg2", "Gad1", "Plp1",
                      "Slc4a10", "Pgap1", "Dlx1", "Dclk2", "Cntnap4", "Sema5a",
                      "Sema6a", "Pcdh7", "Cdh4", "Cpa6")
)

# -----------------------------------------------------------------------------
# 6. Output behaviour
# -----------------------------------------------------------------------------

cfg$output <- list(
  figure_dpi    = 300,
  figure_device = "cairo_pdf",   # vector output for journal submission
  save_png      = TRUE,          # always write a PNG preview next to the PDF
  table_digits  = 6,
  write_session = TRUE
)

invisible(cfg)
