# =============================================================================
# config/palettes.R
#
# Colour schemes, cell-type orders and the small helpers that keep figures
# consistent across modules. Every palette is a *named* character vector so a
# colour follows a category no matter how the factor levels are ordered.
# =============================================================================

# -----------------------------------------------------------------------------
# Cell types
# -----------------------------------------------------------------------------

# Human multi-cohort atlas: six major classes, mirroring the manuscript palette.
pal_celltype_human <- c(
  "Excitatory.neurons"  = "#919ac2",
  "Inhibitory.neurons"  = "#ffac98",
  "Oligodendrocytes"    = "#70a4c8",
  "Astrocytes"          = "#a5a9af",
  "Microglia"           = "#63917d",
  "Endothelial.cells"   = "#dbd1b4"
)

# Mouse Prmt8 cKO atlas: seven major classes (Nature-style qualitative palette).
pal_celltype_mouse <- c(
  "Excitatory neurons" = "#E64B35",
  "Inhibitory neurons" = "#4DBBD5",
  "Oligodendrocyte"    = "#00A087",
  "OPC"                = "#3C5488",
  "Microglia"          = "#F39B7F",
  "Astrocyte"          = "#8491B4",
  "Endothelial cells"  = "#91D1C2"
)

# Alternative mouse palette used in the GO network figure.
pal_celltype_mouse_alt <- c(
  "Excitatory neurons" = "#3B73B9",
  "Inhibitory neurons" = "#D95F5F",
  "Oligodendrocyte"    = "#5AB4AC",
  "OPC"                = "#8C6BB1",
  "Microglia"          = "#C49A45",
  "Astrocyte"          = "#6A994E",
  "Endothelial cells"  = "#E9C46A"
)

# Neuron subtype colours used in the CellChat panels.
pal_neuron_subtype <- c(
  "PN"    = "#DE3533",
  "PVALB" = "#00A0E9",
  "VIP"   = "#8957A1",
  "SST"   = "#78B360"
)

# Interneuron subtype colours.
pal_inhib_subtype <- c(
  "PVALB" = "#0A8CED",
  "SST"   = "#D62729",
  "VIP"   = "#D7517C"
)

# Excitatory layer colours.
pal_excit_subtype <- c(
  "Exc_L2/3" = "#0092BE",
  "Exc_L4"   = "#0058A3",
  "Exc_L5"   = "#037F41",
  "Exc_L6"   = "#263A83"
)

# -----------------------------------------------------------------------------
# Experimental groups
# -----------------------------------------------------------------------------

# Human cohorts.
pal_group_human <- c(
  "NC"  = "#ED2225",
  "AD"  = "#1f78b4",
  "PD"  = "#F2A1A7",
  "MS"  = "#7DC69B",
  "EPI" = "#D8EEFB"
)

# Mouse groups: sequential blue scale, WT darkest.
pal_group_mouse <- c(
  "WT"             = "#08306B",
  "KOF_vector_CNO" = "#2171B5",
  "KOF_Gq_CNO"     = "#6BAED6",
  "KOF"            = "#C6DBEF"
)

# ASD organoid groups.
pal_group_organoid <- c(
  "Control" = "#4DBBD5",
  "ASD"     = "#E64B35"
)

# -----------------------------------------------------------------------------
# Differential expression / enrichment
# -----------------------------------------------------------------------------

pal_direction <- c(Up = "#E64B35", Down = "#4DBBD5", NS = "#BDBDBD")

pal_enrichment <- c(Up = "#A65628", Down = "#1B9E8E")

# Diverging scale used for Z-scored heat maps.
col_zscore <- c(low = "#3B4CC0", mid = "#F7F7F7", high = "#B40426")

# Viridis-like scale used for single-cell heat maps.
col_scaled_expr <- c(low = "#440154", mid = "#3B528B", high = "#FDE725")

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

#' Look up colours for a set of categories, filling gaps from a fallback ramp.
#'
#' @param levels Character vector of categories to colour.
#' @param palette Named colour vector (see the pal_* objects above).
#' @param fallback_hue Optional base hue used to generate colours for any
#'   category missing from `palette`. Use `NA` to leave them grey.
#' @return Named character vector of colours, in the order of `levels`.
pal_for <- function(levels, palette, fallback_hue = NA) {
  levels <- as.character(levels)
  out <- unname(palette[levels])

  missing <- is.na(out) | !nzchar(out)
  if (any(missing)) {
    n_missing <- sum(missing)
    if (is.na(fallback_hue)) {
      out[missing] <- "grey70"
    } else {
      out[missing] <- grDevices::hcl(
        h = fallback_hue,
        c = 55,
        l = seq(70, 40, length.out = n_missing)
      )
    }
    warning(
      "No colour defined for: ", paste(levels[missing], collapse = ", "),
      " -- using a fallback.", call. = FALSE
    )
  }

  stats::setNames(out, levels)
}

#' Order a factor using the project-wide level order, appending unknown levels.
ordered_factor <- function(x, levels) {
  present <- unique(as.character(x))
  lv <- c(intersect(levels, present), setdiff(present, levels))
  factor(as.character(x), levels = lv)
}

#' Fill-in helper: replace NA colours in a ggplot scale by recycling.
recycle_colours <- function(n) {
  base <- unlist(list(
    "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
    "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC"
  ), use.names = FALSE)
  rep(base, length.out = n)
}
