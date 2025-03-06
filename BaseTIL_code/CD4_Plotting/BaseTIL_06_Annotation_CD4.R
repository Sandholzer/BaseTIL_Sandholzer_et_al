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
setwd("~/BaseTIL_code/CD4_Plotting")

gex.CD4 <- readRDS(file = "../saveFiles/CD4_integrated_3.rds")


# merge specific cluster to reflect biological phenotype
gex.CD4@meta.data[gex.CD4@meta.data$RNA_snn_res.0.6 == 13 ,"seurat_clusters"] <- 9
gex.CD4@meta.data[gex.CD4@meta.data$RNA_snn_res.0.6 == 9 ,"seurat_clusters"] <- 7
gex.CD4 <- Seurat::SetIdent(gex.CD4, value = gex.CD4$seurat_clusters)


gex.CD4 <-
  RenameIdents(
    gex.CD4,
    `1` = "Tm",
    `2` = "Teff",
    `3` = "Tn",
    `4` = "Treg",
    `5` = "Th1/Tex",
    `6` = "Treg",
    `7` = "Th17",
    `8` = "Tprol",
    `9` = "Tfh",
    `10` = "Th2")

gex.CD4$Celltype.fine <- Idents(gex.CD4)

gex.CD4 <- Seurat::SetIdent(gex.CD4, value = gex.CD4$seurat_clusters)
gex.CD4 <-
  RenameIdents(
    gex.CD4,
    `1` = "Tm",
    `2` = "Teff",
    `3` = "Tn",
    `4` = "Treg",
    `5` = "Th1/Tex",
    `6` = "Treg",
    `7` = "Th17",
    `8` = "Tprol",
    `9` = "Tfh",
    `10` = "Th2")
gex.CD4$Celltype <- Idents(gex.CD4)

gex.CD4$Celltype <- factor(gex.CD4$Celltype, levels = c("Tn", "Tm", "Teff", "Th1/Tex", "Tfh", "Th2", "Treg", "Th17", "Tprol"))




pdf(file = "Fig/CD4_Annotated_UMAP.pdf", height = 8, width = 10)
DimPlot(gex.CD4, reduction = "umap", group.by = "Celltype", label = T, repel = F, raster = F, cols = Colors, label.size = 11,label.box = T)+NoLegend()+ ggeasy::easy_remove_axes() + ggtitle(NULL)
dev.off()


# 
# gex.CD4$prepost <-  ifelse(as.character(gex.CD4$Type_new) %in% c("post-ACT.1", "post-ACT.2"), "post-ACT", as.character(gex.CD4$Type_new))
# gex.CD4$prepost <- factor(gex.CD4$prepost, level = c("pre-ACT", "Intermediate product", "TIL product", "PBMC 7dpt", "post-ACT"))

png(filename = "Graphs/CD4_Umap_splitby.png", height = 400, width = 2200, type = "cairo")
DimPlot(gex.CD4, reduction = "umap", split.by = "prepost", group.by = "Celltype", label = T, repel = F, raster = F, cols = Colors, label.size = 7)+NoLegend()+ ggeasy::easy_remove_axes() + ggtitle(NULL)
dev.off()

pdf("Fig/CD4_Umap_splitby.pdf", height = 4, width = 22)
DimPlot(gex.CD4, reduction = "umap", split.by = "prepost", group.by = "Celltype",
        label = F, repel = F, raster = F, cols = Colors, label.size = 5.5,label.box = T, alpha = 0.65) &
  #NoLegend() & 
  ggeasy::easy_remove_axes() & 
  ggtitle(NULL) & theme(strip.text.x = element_text(size=20, face = "bold"),
                        legend.text = element_text(size=20))
dev.off()



annotation_markersCD4 <- c("CCR7", "LEF1","TCF7", "SELL", "IL7R",
                           "CD69", "FOS", "JUN",
                           "GNLY","PRF1","GZMK","TNF", "IFNG", "NKG7",
                           "LAG3", "HAVCR2","PDCD1",  "TOX","ITGAE","TIGIT",
                           "CXCL13","TOX2",
                           "IL17RB", "IL13","CCR3", 
                           "FOXP3","CTLA4", "IL2RA", "IKZF2", "TNFRSF9",
                           "RORC", "IL17A", "CCL20","IL23R", "IL26",
                           "RUNX1", "MKI67", "TOP2A")

pdf("Fig/CD4_Annotation_markers.pdf", width = 14, height = 5)
DotPlot(gex.CD4, features = annotation_markersCD4, assay = "RNA", group.by= "Celltype", cols= "RdBu",col.min = -1.7, col.max = 1.7)+
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 1.0,
    hjust = 1
  )
  ) +
  xlab(NULL) +
  ylab(NULL)+scale_y_discrete(limits= rev(c("Tn", "Tm", "Teff", "Th1/Tex", "Tfh", "Th2", "Treg", "Th17", "Tprol")))
#+ coord_flip()
dev.off()



gc()

saveRDS(gex.CD4, file = "../saveFiles/CD4_annotated.rds")




DefaultAssay(gex.CD4) <- "RNA"
gex.CD4 <- Seurat::SetIdent(gex.CD4, value = gex.CD4$Celltype)

plt <- list()
allMarkers <- FindAllMarkers(gex.CD4,  only.pos = T,  assay = "RNA") 

allMarkers <- subset(allMarkers, allMarkers$p_val_adj<0.05)
write.csv(allMarkers, file = "Fig/CD4_FindAllMarkers_ann.csv", row.names = F, quote = F)

allMarkers %>%
  group_by(cluster) %>%
  #dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 3) %>%
  ungroup() -> top3

plt <- dittoSeq::dittoDotPlot(gex.CD4, vars = unique(top3$gene),
                              group.by = "Celltype")


pdf(file = "Graphs/CD4_dottplot_top3_Markers_ann.pdf", width = 15, height = 8)
plt
dev.off()



allMarkers %>%
  group_by(cluster) %>%
  #dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10


tmp.sup <- subset(x = gex.CD4, downsample = 250)
tmp.sup <- ScaleData(tmp.sup, features = top10$gene)

png(file = paste0("Graphs/CD4_Heatmap_Markers_Annotation.png"), width = 2000, height = 3000, res = 120, type = "cairo")
DoHeatmap(tmp.sup, features = top10$gene, group.by = "Celltype" ) + NoLegend()
dev.off()


