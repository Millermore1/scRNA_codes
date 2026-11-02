# results/

Everything the pipeline produces lands here. The folder is git-ignored; only
the empty structure is tracked.

```
results/
├── figures/    <module>_<name>.pdf   vector, cairo_pdf, for submission
│               <module>_<name>.png   300 dpi preview
├── tables/     <module>_<name>.csv   UTF-8, one header row
└── objects/    <name>.rds            intermediate objects, read by later modules
```

The numeric prefix of every file is the module that produced it, so a figure can
always be traced to its script:

| Prefix | Module |
|---|---|
| 01 | `R/01_prepare_datasets.R` |
| 02 | `R/02_annotation.R` |
| 03 | `R/03_cell_composition.R` |
| 04 | `R/04_differential_expression.R` |
| 05 | `R/05_inhibitory_neurons.R` |
| 06 | `R/06_excitatory_neurons.R` |
| 07 | `R/07_enrichment_networks.R` |
| 08 | `R/08_cellchat.R` |
| 09 | `R/09_wgcna.R` |
| 10 | `R/10_figures_heatmaps.R` |
| 11 | `R/11_asd_organoid.R` |

`figures/` is deleted at the start of a full run only if you delete it — the
pipeline overwrites files by name, so a partial re-run will not erase panels
from modules you did not re-run.

`docs/figure_mapping.md` lists every expected output with the manuscript panel
it corresponds to.
