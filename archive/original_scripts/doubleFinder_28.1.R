library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(devtools)
library(remotes)
library(CellChat)
library(DoubletFinder)
# cellchat & doublefinder installment
Sys.getenv(c("GITHUB_TOKEN", "GITHUB_PAT"))
Sys.unsetenv("GITHUB_TOKEN")
Sys.unsetenv("GITHUB_PAT")
remotes::install_github(
  "chris-mcginnis-ucsf/DoubletFinder",
  auth_token = NULL
)
remotes::install_github("sqjin/CellChat")

remove.packages("DoubletFinder", lib = "E:/R/R-4.3.0/library")

## 确认已经卸载掉
find.package("DoubletFinder")
## 这一步应该报错：there is no package called 'DoubletFinder'


Sys.unsetenv("GITHUB_TOKEN")
Sys.unsetenv("GITHUB_PAT")

remotes::install_github(
  "chris-mcginnis-ucsf/DoubletFinder",
  force           = TRUE,
  upgrade         = "never",
  dependencies    = FALSE,
  build_vignettes = FALSE
)

library(DoubletFinder)

packageVersion("DoubletFinder")
ls("package:DoubletFinder")

DimPlot(scRNA_HKU)


# 随便选 5–6 个存在的基因
marker_genes <- rownames(scRNA_HKU)[1:6]
marker_genes

Idents(scRNA_HKU) <- "celltype"
Idents(scRNA_HKU)
DotPlot(
  scRNA_HKU,
  features = marker_genes,
  group.by = "celltype"
) + RotatedAxis()
table(scRNA_HKU$celltype)
marker_genes <- c(
  # Excitatory neurons
  "Slc17a7", "Camk2a", "Satb2",
  # Inhibitory neurons
  "Gad1", "Gad2", "Slc32a1",
  # Oligodendrocyte
  "Mbp", "Mog", "Plp1",
  # Microglia
  "P2ry12", "Cx3cr1", "C1qa",
  # OPC
  "Pdgfra", "Cspg4", "Olig1",
  # Astrocyte
  "Gfap", "Aqp4", "Aldh1l1",
  # Endothelial cells
  "Cdh5", "Kdr", "Flt1"
)

Idents(scRNA_HKU) <- "celltype"

DotPlot(
  scRNA_HKU,
  features = marker_genes,
  group.by = "celltype"
) + RotatedAxis()

DotPlot(
  scRNA_HKU,
  features = marker_genes,
  group.by = "celltype",
  dot.scale = 9,                     # 点变大
  cols = c("#D5EAED", "#0787C3")       # 低表达到高表达的颜色
) + RotatedAxis()






