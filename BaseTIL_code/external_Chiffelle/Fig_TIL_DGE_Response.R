## -----------------------------------------------------------------------------------------------------------
#load dependencies

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(Seurat)
  library(EnhancedVolcano)
  library(stringr)
  
})

set.seed(12345678)

setwd("~/BaseTIL_code/external_Chiffelle")

acttil <- readRDS("saveFiles/Chiffelle_TIL_act.rds")



## -----------------------------------------------------------------------------------------------------------
#Differential expression analysis

Idents(acttil) <- acttil$RNR

markers <- FindMarkers(acttil, ident.1 = "R", ident.2 = "NR", logfc.threshold = 0, min.pct = 0.005, min.diff.pct = 0)
markers$gene <- rownames(markers)

write_csv(markers, "Graphs/dge_TIL_RNR.csv")


#filter irrelevant genes
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

#select genes to annotate
show.genes_pos <- c("PDCD1", "CTLA4", "HLA-DRB5", "ITGAE","ETV1", "KIR2DL4", "VCAM1","TOX", "KLRC4","FCER2", "KIR2DL1", "TOX2")
show.genes_neg <- c("IL26","RORC", "KLRB1", "IL17RE", "CCL3", "CCL4", "CX3CR1", "ISG15", "CD69", "AQP3", "CXCR6", "IFIT1", "IFI27", "IRF7", "LTB", "GNLY", "LAIR2", "GZMK", "MX1")



markers.flt$selection <- ifelse(markers.flt$gene %in% c(show.genes_pos, show.genes_neg),markers.flt$gene, NA)


#plot the graph
pdf("Fig/DGE_TIL_RNR.pdf", width = 5, height = 5)

EnhancedVolcano::EnhancedVolcano(
  markers.flt,
  x = "avg_log2FC",
  y = "p_val",
  rownames(markers.flt),
  FCcutoff = 0.5,
  pCutoff = 1e-15,
  xlim = c(-4.5, 4.5),
  ylim = c(0, 310),
  selectLab = NA,
  parseLabels = T,
  drawConnectors = T, 
  subtitle = "", caption = "", 
  title = "", 
  pointSize = 1.0,
  labSize = 3.0,
  col=c( "#D6D6D6", "#D6D6D6", "#D6D6D6","#CD2626"),
  colAlpha = 1,
  boxedLabels = T, 
  max.overlaps = 2, 
  arrowheads = F, 
  maxoverlapsConnectors =
    15) & 
  theme_classic() & NoLegend() & 
  theme(
    legend.position = "none",
    axis.title = element_text(size = 15,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) & geom_label_repel(aes(label = selection),
                       parse = T,
                       size = 3,
                       box.padding = 0.25,
                       segment.color = "black",min.segment.length = 0,
                       show.legend = FALSE, max.overlaps = 15, point.padding = 0, label.padding = 0.15)


EnhancedVolcano(markers.flt, x = "avg_log2FC", y="p_val_adj", lab = rownames(markers.flt), 
                pCutoff = 1e-18, 
                FCcutoff = 0.5, drawConnectors = T, 
                #selectLab = c(show.genes_pos, show.genes_neg),
                pointSize = 1.0,
                labSize = 4.0,
                col=c( "#D6D6D6", "#D6D6D6", "#D6D6D6","#CD2626"),
                colAlpha = 1,
                title = "",
                subtitle = "", boxedLabels = F, max.overlaps = 2, arrowheads = F, maxoverlapsConnectors=20, parseLabels=T, 
                xlim = c(-4.5,4.5)
                ) + NoLegend() +  theme_classic() +
  scale_x_continuous(limits = c(-4.5,4.5), breaks = seq(-4.5, 4.5, 1.5)) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 15,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  )
dev.off()

