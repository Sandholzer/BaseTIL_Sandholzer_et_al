suppressPackageStartupMessages({
  library(Seurat) 
  library(ggplot2)
  library(dplyr)
  library(cowplot)
  library(ggpubr)
  library(readr)
  library(ggeasy)
  library(readr)
  library(dittoSeq)
  library(RColorBrewer)
  
  Colors <- RColorBrewer::brewer.pal(n=12, "Paired")
})
set.seed(1234567)
setwd("~/BaseTIL_code/CD8_Plotting")

gex.CD8 <- readRDS(file = "../saveFiles/CD8_integrated_3.rds")


gex.CD8 <- Seurat::SetIdent(gex.CD8, value = gex.CD8$seurat_clusters)

gex.CD8 <-
  RenameIdents(
    gex.CD8,
    `1` = "Teff",
    `2` = "Tm",
    `3` = "Tprol1",
    `4` = "Temra",
    `5` = "Tex",
    `6` = "Tem",
    `7` = "Tprol2",
    `8` = "Tprol3",
    `9` = "Tn",
    `10` = "FOXP3",
    `11` = "Tc17",
    `12` = "Teff")


gex.CD8$Celltype.fine <- Idents(gex.CD8)


gex.CD8 <- Seurat::SetIdent(gex.CD8, value = gex.CD8$seurat_clusters)

gex.CD8 <-
  RenameIdents(
    gex.CD8,
    `1` = "Teff",
    `2` = "Tm",
    `3` = "Tprol1",
    `4` = "Temra",
    `5` = "Tex",
    `6` = "Tem",
    `7` = "Tprol2",
    `8` = "Tprol3",
    `9` = "Tn",
    `10` = "FOXP3",
    `11` = "Tc17",
    `12` = "Teff")


gex.CD8$Celltype <- Idents(gex.CD8)
gex.CD8$Celltype <- factor(gex.CD8$Celltype ,levels = c("Tn", "Tm", "Tem", "Teff", "Temra", "Tex", "Tc17", "FOXP3", "Tprol1", "Tprol2", "Tprol3"))


gex.CD8 <- Seurat::SetIdent(gex.CD8, value = gex.CD8$seurat_clusters)

gex.CD8 <-
  RenameIdents(
    gex.CD8,
    `1` = "Teff",
    `2` = "Tm",
    `3` = "Tprol1",
    `4` = "Temra",
    `5` = "Tex",
    `6` = "Tem",
    `7` = "Tprol2",
    `8` = "Tprol3",
    `9` = "Tn",
    `10` = "FOXP3",
    `11` = "Tc17",
    `12` = "Teff")

gex.CD8$Celltype.simple <- Idents(gex.CD8)


pdf(file = "Fig/CD8_Annotated_UMAP.pdf", height = 8, width = 10)
DimPlot(gex.CD8, reduction = "umap", group.by = "Celltype", label = T, repel = F, raster = F, cols = Colors, label.size = 11,label.box = T)+NoLegend()+ ggeasy::easy_remove_axes() + ggtitle(NULL)
dev.off()

gex.CD8$prepost <-  ifelse(as.character(gex.CD8$Type_new) %in% c("post-ACT.1", "post-ACT.2"), "post-ACT", as.character(gex.CD8$Type_new))
gex.CD8$prepost <- factor(gex.CD8$prepost, level = c("pre-ACT", "Intermediate product", "TIL product", "PBMC 7dpt", "post-ACT"))

png(filename = "Graphs/CD8_Umap_splitby.png", height = 400, width = 2200, type = "cairo")
DimPlot(gex.CD8, reduction = "umap", split.by = "prepost", group.by = "Celltype", label = T, repel = F, raster = F, cols = Colors, label.size = 7,label.box = T)+NoLegend()+ ggeasy::easy_remove_axes() + ggtitle(NULL)
dev.off()


pdf("Fig/CD8_Umap_splitby.pdf", height = 4, width = 22)
DimPlot(gex.CD8, reduction = "umap", split.by = "prepost", group.by = "Celltype",
        label = F, repel = F, raster = F, cols = Colors, label.size = 5.5,label.box = T, alpha = 0.65) &
  ggeasy::easy_remove_axes() & 
  ggtitle(NULL) & theme(strip.text.x = element_text(size=20, face = "bold"),
                        legend.text = element_text(size=20))
dev.off()



saveRDS(gex.CD8, file = "../saveFiles/CD8_annotated.rds")


annotation_markersCD8 <- c("CCR7", "LEF1","TCF7", "SELL", "IL7R",
                           "ZNF683", "CD52", "GNLY","PRF1","GZMK", 
                           "CD69", "TNF", "IFNG", "CXCR3",
                           "CX3CR1", "TBX21", "ASCL2", "KLRG1","FGFBP2","FCGR3A",
                           "LAG3", "PDCD1", "HAVCR2", "TOX","ITGAE","TIGIT","CXCR6",
                           "RORC", "KLRB1", "CCR6","IL23R", "IL26",
                           "FOXP3","CTLA4", "IL2RA", "IKZF2",
                           "MKI67", "TOP2A")

pdf("Fig/CD8_Annotation_markers.pdf", width = 12, height = 5)
DotPlot(gex.CD8, features = annotation_markersCD8, assay = "RNA", group.by= "Celltype", cols= "RdBu",col.min = -1.6, col.max = 1.6)+
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 1.0,
    hjust = 1
  )
  ) +
  xlab(NULL) +
  ylab(NULL)+scale_y_discrete(limits= rev(c("Tn", "Tm", "Tem", "Teff", "Temra", "Tex", "Tc17", "FOXP3", "Tprol1", "Tprol2", "Tprol3")))
#+ coord_flip()
dev.off()


DefaultAssay(gex.CD8) <- "RNA"
gex.CD8 <- Seurat::SetIdent(gex.CD8, value = gex.CD8$Celltype)

plt <- list()
allMarkers <- FindAllMarkers(gex.CD8,  only.pos = T,  assay = "RNA") 

allMarkers <- subset(allMarkers, allMarkers$p_val_adj<0.05)
write.csv(allMarkers, file = "Fig/CD8_FindAllMarkers_ann.csv", row.names = F, quote = F)

allMarkers %>%
  group_by(cluster) %>%
  #dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 3) %>%
  ungroup() -> top3

plt <- dittoSeq::dittoDotPlot(gex.CD8, vars = unique(top3$gene),
                              group.by = "Celltype")


pdf(file = "Graphs/CD8_dottplot_top3_Markers_ann.pdf", width = 15, height = 8)
plt
dev.off()



allMarkers %>%
  group_by(cluster) %>%
  #dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10


tmp.sup <- subset(x = gex.CD8, downsample = 250)
tmp.sup <- ScaleData(tmp.sup, features = top10$gene)

png(file = paste0("Graphs/CD8_Heatmap_Markers_An.png"), width = 2000, height = 3000, res = 120, type = "cairo")
DoHeatmap(tmp.sup, features = top10$gene, group.by = "Celltype" ) + NoLegend()
dev.off()
