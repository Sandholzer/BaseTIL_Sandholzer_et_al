## -----------------------------------------------------------------------------------------------------------
#load dependencies

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(stringr)
  library(ggplot2)
  library(ggalluvial)
  library(ggpubr)
  library(readr)
  library(ggeasy)
  library(purrr)
  library(ComplexHeatmap)
  library(colorRamp2)
})

set.seed(1234567)
setwd("~/BaseTIL_code/CD4_Plotting")

gex.CD4 <- readRDS(file = "../saveFiles/CD4_annotated.rds")

#remove control PBMC sample from analysis
gex.CD4 <- subset(gex.CD4, Type != "PBMC_Ctrl")


#Aggregate per cell type
aggro_CD4 <- AggregateExpression(gex.CD4, return.seurat = T, group.by = "Celltype")



#Load for Annotation with Pan-cancer Zheng et al. 
  Zheng_CD4 <- read_csv("../data/Pan_cancer_Zheng_CD4.csv")
  load("../data/Zheng_Japrin/exclude.gene.misc.human.v4.RData")
  
  
  Zheng_high_ES <- subset(Zheng_CD4, subset =  Zheng_CD4$comb.positive.freq > 20) #Filter based on Enrichment Score
  f.feat <- !(Zheng_high_ES$geneSymbol %in% all.gene.ignore.df[["seu.id"]]) &
    !(grepl("^(RP[LS]|Rp[ls])",Zheng_high_ES$geneSymbol,perl=T)) &
    !(grepl("^MT-|^MTRN|^MTAT|^MTND|^MRP",Zheng_high_ES$geneSymbol,perl=T)) &
    !(grepl("^RP[0-9]+-",Zheng_high_ES$geneSymbol,perl=T)) &
    !(Zheng_high_ES$geneSymbol %in% c("MALAT1", 'XIST', 'NEAT1', "LINC02694", "LINC01619", "LINC02506", "AL136456.1"))
  Zheng_high_ES <- Zheng_high_ES[f.feat,]
  Zheng_high_ES %>% group_by(cluster.name) %>% top_n(n = 15, wt = comb.ES) -> Zheng_high_ES
  gene_list <- split(Zheng_high_ES, f=Zheng_high_ES$cluster.name) #Create per Cluster list object
  gene_list <- lapply(gene_list, FUN = function(X){X <- X$geneSymbol}) #Take only geneSymbols for ModuleScore
  
  #Calculate UCell score per signature
  aggro_CD4 <- AddModuleScore_UCell(aggro_CD4, features = gene_list, name = "")
 

  mat<- aggro_CD4@meta.data[,names(gene_list)] %>% as.matrix()
  

mat<- t(scale(mat))

quantile(mat, c(0.05, 0.95))

library(circlize)
library(ComplexHeatmap)

pdf("Fig/CD4_Cluster_annotation_heatmap.pdf", height = 8, width = 6)
Heatmap(mat, name = "Z score",  
        show_column_dend = T,
        cluster_columns = F,
        cluster_column_slices = TRUE,
        cluster_rows = F,
        column_title_gp = gpar(fontsize = 12),
        column_gap = unit(0.75, "mm"),
        show_row_dend = FALSE,
        col=colorRamp2(c(-1.5, 0, 1.7), c("#3C4DA0", "white", "#CD1F2D")),
        row_names_gp = gpar(fontsize = 12),
        column_title_rot = 45,
        column_names_rot = 45,
        column_names_side = "top",
        show_column_names = T,
        use_raster = TRUE,
        raster_quality = 4,
        row_names_side = "left",
)
dev.off()
