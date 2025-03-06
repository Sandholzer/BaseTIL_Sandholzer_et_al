## -----------------------------------------------------------------------------------------------------------
#load dependencies

suppressPackageStartupMessages({
  library(scRepertoire)
  library(readr)
  library(dplyr)
  library(cowplot)
  library(ggpubr)
  library(tidyr)
  library(ggeasy)
  library(Seurat)
  library(ggplot2)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(circlize)
  library(ggplotify)
  library(EnhancedVolcano)
})

set.seed(1234567)
setwd("~/BaseTIL_code/CD4_Plotting/")

gex.CD4 <- readRDS(file = "../saveFiles/CD4_annotated.rds")
gex.CD4$RNR <- ifelse(gex.CD4$Patient %in% c("UPN001", "UPN003"), "R", "NR")


## -----------------------------------------------------------------------------------------------------------
#Differential gene expression of Treges in NR pre and post-ACT

#Subset 
gex.treg <- subset(gex.CD4, Celltype %in% c("Treg","activeTreg") & RNR == "NR" & Type %in% c("Tumor", "Rebiopsy1", "Rebiopsy2")) #
gex.treg$PrePost <- ifelse(gex.treg$Type %in% c("Rebiopsy1", "Rebiopsy2"), "Post", "Pre")

Idents(gex.treg) <- gex.treg$PrePost

#Differential gene expression
markers <- FindMarkers(gex.treg, ident.1 = "Post", ident.2 = "Pre", logfc.threshold = 0) 
markers$gene <- rownames(markers)

write_csv(markers, "Graphs/dge_Treg_NR_PrevsPost.csv")
markers <- read_csv("Graphs/dge_Treg_NR_PrevsPost.csv")


#Filter unwanted genes
load("~/BaseTIL_code/data/Zheng_Japrin/exclude.gene.misc.human.v4.RData")
f.feat <- !(markers$gene %in% all.gene.ignore.df[["seu.id"]]) &
  !(grepl("^(RP[LS]|Rp[ls])",markers$gene,perl=T)) &
  !(grepl("^MT-|^MTRN|^MTAT|^MTND|^MRP",markers$gene,perl=T)) &
  !(grepl("^RP[0-9]+-",markers$gene,perl=T)) &
  !(markers$gene %in% c("MALAT1", 'XIST', 'NEAT1', "LINC02694", "LINC01619", "LINC02506", "AL136456.1", "MIR4435-2HG"))
markers.flt <- markers[f.feat,]
print(head(markers.flt))
print(paste("Length of anchors:", nrow(markers)) )
print(paste("Length of filtered anchors:", nrow(markers.flt)) )

markers.flt$gene[markers.flt$p_val_adj < 0.01 & markers.flt$avg_log2FC > 0.5]

#Select genes for annotation
show.genes_pos <- c("CCR1", "CCR5", "CCR8", "PDCD1","CTLA4", "LGALS1", "TNFRSF4", "TNFRSF18","PRF1","GZMA", "IL10", "IL32", "DUSP2", "S100A11", "S100A4", "S100A6", "CCL5", "CCL4", "CD33", "IL2RG", "CXCR3","CXCR6", "ENTPD1", "ITGAM")
show.genes_neg <- c("MAL", "LEF1", "SELL", "TCF7", "S1PR1", "MYB", "CCR7", "IL7R", "CH25H", "FOXP1", "BACH2")
markers.flt$selection <- ifelse(markers.flt$gene %in% c(show.genes_pos, show.genes_neg),markers.flt$gene, NA)


#plot the graph
pdf("Fig/DGE_Treg_NR_PrevsPost.pdf",
    width = 5,
    height = 5)
EnhancedVolcano::EnhancedVolcano(
  markers.flt,
  x = "avg_log2FC",
  y = "p_val",
  rownames(markers.flt),
  FCcutoff = 0.25,
  pCutoff = 10e-20,
  xlim = c(-3.5, 3.5),
  ylim = c(0, 300),
  selectLab = NA,
  parseLabels = T,
  drawConnectors = T,
  subtitle = "",
  caption = "",
  title = "",
  pointSize = 1.0,
  labSize = 3.0,
  col = c("#D6D6D6", "#D6D6D6", "#D6D6D6", "#CD2626"),
  colAlpha = 1,
  boxedLabels = T,
  max.overlaps = 2,
  arrowheads = F,
  maxoverlapsConnectors =
    15
) &
  theme_classic() & NoLegend() &
  theme(
    legend.position = "none",
    axis.title = element_text(size = 15, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color = "black"),
    axis.text.x = element_text(color = "black")
  ) & geom_label_repel(
    aes(label = selection),
    parse = T,
    size = 3,
    box.padding = 0.25,
    segment.color = "black",
    min.segment.length = 0,
    show.legend = FALSE,
    max.overlaps = 15,
    point.padding = 0,
    label.padding = 0.15
  )
dev.off()

