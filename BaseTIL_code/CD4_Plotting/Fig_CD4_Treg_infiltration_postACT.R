## -----------------------------------------------------------------------------------------------------------
#load dependencies

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggeasy)
})

set.seed(1234567)
setwd("~/BaseTIL_code/CD4_Plotting/")

gex.CD4 <- readRDS(file = "../saveFiles/CD4_annotated.rds")



## -----------------------------------------------------------------------------------------------------------
#highlights the Tregs found in PBMC 7dpt in post-ACT lesions

plt <- list()
tracked_cluster <- "Treg"
selected_type <- "PBMC_7dpt"


clones <- gex.CD4@meta.data[gex.CD4@meta.data$Type == selected_type & gex.CD4@meta.data$Celltype == tracked_cluster, "CTaa"]
clones <-  clones[!is.na(clones)]

gex.CD4$Tracked_clones <- "None"
gex.CD4$Tracked_clones <- gex.CD4$CTaa %in%  clones


pdf("Fig/Treg_infiltration.pdf", width = 4, height = 4)
DimPlot(
  subset(gex.CD4, Type  %in% c("Rebiopsy1", "Rebiopsy2")),
  label = F,
  cells.highlight = WhichCells(gex.CD4, expression = Tracked_clones != FALSE),
  sizes.highlight = 0.75, raster=FALSE
) +
  scale_color_manual(
    labels = c("Unselected", "NeoTCRs"),
    values = c("grey80",  "#F8766D")
  ) + NoLegend() + easy_remove_axes() + ggtitle("Treg infiltration in post-ACT") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()



## -----------------------------------------------------------------------------------------------------------
#plots CXCR3 expression in PBMC 7dpt


pdf("Fig/CXCR3_in_Blood.pdf",  height = 4 ,width = 5)
Nebulosa::plot_density(subset(gex.CD4, Type %in% c("PBMC_7dpt")), features = "CXCR3") & 
  ggtitle("CXCR3 expression\nin PBMC 7dpt") &
  ggeasy::easy_remove_axes() & 
  ggeasy::easy_remove_legend_title() &
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 15))

FeaturePlot(subset(gex.CD4, Type %in% c("PBMC_7dpt")), features = "CXCR3") + 
  easy_remove_axes() + 
  ggtitle("CXCR3 expression\nin PBMC 7dpt")  
dev.off()

