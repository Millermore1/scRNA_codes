pv <- subset(inh, subset = inhibit_celltype == "PVALB")
colnames(pv@meta.data)         
table(pv$group)               
Idents(pv) <- "group"         
pv.DEG <- FindMarkers(
  object = pv,
  ident.1 = "KOF_Gq_CNO",
  ident.2 = "KOF_vector_CNO",
  min.pct = 0.1,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)
rm(mycolors,group_col,cell.cols,p)
table(inh$group)

head(pv.DEG)
deg.sig <- pv.DEG[
  pv.DEG$p_val_adj < 0.05 & abs(pv.DEG$avg_log2FC) > 0.25, 
]
deg.sig <- cbind(gene = rownames(deg.sig), deg.sig)
write.csv(
  deg.sig,
  file = "PV_KOF_Gq_vs_Vector_DEGs.csv",
  row.names = FALSE
)
deg.up.KOF <- subset(deg.sig, avg_log2FC > 0)
deg.up.WT  <- subset(deg.sig, avg_log2FC < 0)

write.csv(deg.up.KOF, "PV_DEGs_up_in_KOF_Gq.csv", row.names = FALSE)
write.csv(deg.up.WT,  "PV_DEGs_up_in_vector.csv",  row.names = FALSE)




library(clusterProfiler)
library(org.Mm.eg.db)   # 小鼠 GO 注释
library(ggplot2)


##====================================================##
## 1. 从 Seurat 差异分析结果中筛选显著 DEGs
##    假设你的差异结果叫 pv.DEG (FindMarkers 的输出)
##====================================================##

# 根据需要修改阈值：p_val_adj、avg_log2FC
deg.sig <- pv.DEG[
  pv.DEG$p_val_adj < 0.05 & abs(pv.DEG$avg_log2FC) > 0.25,
]

# 把行名（基因名）变成一列 gene
deg.sig <- cbind(gene = rownames(deg.sig), deg.sig)
rownames(deg.sig) <- NULL

# 分成 KOF 上调、WT 上调两部分（可选）
deg.up.KOF <- subset(deg.sig, avg_log2FC > 0)   # 在 KOF 组上调
deg.up.WT  <- subset(deg.sig, avg_log2FC < 0)   # 在 WT 组上调


##====================================================##
## 2. 基因 ID 转换：Symbol -> EntrezID (clusterProfiler 推荐)
##====================================================##

##---------- 2.1 KOF 上调基因 ----------##
genes_up_KOF <- unique(deg.up.KOF$gene)

gene_up_KOF_entrez <- bitr(
  genes_up_KOF,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Mm.eg.db
)

# 去掉可能重复的 ENTREZID
gene_up_KOF_entrez <- gene_up_KOF_entrez[!duplicated(gene_up_KOF_entrez$ENTREZID), ]

##---------- 2.2 WT 上调基因（可选） ----------##
genes_up_WT <- unique(deg.up.WT$gene)

gene_up_WT_entrez <- bitr(
  genes_up_WT,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Mm.eg.db
)

gene_up_WT_entrez <- gene_up_WT_entrez[!duplicated(gene_up_WT_entrez$ENTREZID), ]


##====================================================##
## 3. GO 富集分析（小鼠）
##    这里以 Biological Process(BP) 为主
##    想看 CC / MF：把 ont 改成 "CC" 或 "MF"
##====================================================##

##---------- 3.1 KOF 上调基因 GO 富集 ----------##
ego_up_KOF_BP <- enrichGO(
  gene          = gene_up_KOF_entrez$ENTREZID,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",      # "BP" / "CC" / "MF" / "ALL"
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE       # 把 ENTREZID 显示为 Symbol
)

# 查看前几条结果
head(as.data.frame(ego_up_KOF_BP))

##---------- 3.2 WT 上调基因 GO 富集（可选） ----------##
ego_up_WT_BP <- enrichGO(
  gene          = gene_up_WT_entrez$ENTREZID,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

head(as.data.frame(ego_up_WT_BP))


##====================================================##
## 4. 可视化（点图 / 柱状图）
##====================================================##

# KOF 上调 BP 富集 dotplot
pdf("GO_BP_dotplot_KOF_up_in_PV.pdf", width = 8, height = 6)
dotplot(ego_up_KOF_BP, showCategory = 20) +
  ggtitle("GO BP enrichment (KOF up in PV)")
dev.off()

# WT 上调 BP 富集 dotplot（可选）
pdf("GO_BP_dotplot_WT_up_in_PV.pdf", width = 8, height = 6)
dotplot(ego_up_WT_BP, showCategory = 20) +
  ggtitle("GO BP enrichment (WT up in PV)")
dev.off()

# 也可以画 barplot（可选）
# barplot(ego_up_KOF_BP, showCategory = 20) + ggtitle("GO BP (KOF up)")
# barplot(ego_up_WT_BP,  showCategory = 20) + ggtitle("GO BP (WT up)")


##====================================================##
## 5. 结果导出为 CSV
##====================================================##

# KOF 上调基因 GO BP 富集结果
ego_up_KOF_BP_df <- as.data.frame(ego_up_KOF_BP)
write.csv(ego_up_KOF_BP_df,
          file = "GO_BP_enrich_KOF_up_in_PV.csv",
          row.names = FALSE)

# WT 上调基因 GO BP 富集结果
ego_up_WT_BP_df <- as.data.frame(ego_up_WT_BP)
write.csv(ego_up_WT_BP_df,
          file = "GO_BP_enrich_WT_up_in_PV.csv",
          row.names = FALSE)

