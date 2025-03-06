## -----------------------------------------------------------------------------------------------------------
#load dependencies

suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(ggplot2)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(circlize)
  library(GSEABase)
})

set.seed(1234567)
setwd("~/BaseTIL_code/CD8_Plotting")

gex.CD8 <- readRDS(file = "../saveFiles/CD8_annotated.rds")

## -----------------------------------------------------------------------------------------------------------
#subset only tumor-reactive T cells and remove sampels with to low numbers
TR.CD8 <- subset(gex.CD8, overall_reactive == TRUE & Patient %in% c("UPN001","UPN002", "UPN003", "UPN006" , "UPN011", "UPN009", "UPN008")) #
TR.CD8 <- subset(TR.CD8, Sample_Name %in%  c("UPN008 Rebiopsy1", "UPN001 Rebiopsy2"), invert = T)


# aggregate based on Patient and timepoint
aggro <- AggregateExpression(TR.CD8, return.seurat = T, group.by = c("Patient", "Type"))

aggro$Type_post <- ifelse(aggro$Type %in% c("Rebiopsy1", "Rebiopsy2"), "Post", aggro$Type)
Idents(aggro) <- aggro$Type_post


# define markers to show
tumor_markers <- c("CXCL13", "SNX9", "CTLA4", "NR4A2", "NR4A3", "ZNF331", "PDCD1",  "TNFAIP3", "TOX")
preREP_markers <- c("CD38", "MX1", "IL2RA", "IL18R1", "CD300C",  "HLA-DQA1", "HLA-DOA","HLA-DQB1", "CD70", "CD74", "CD80", "CD86")
ExpTIL_markers <- c("CD33","ZNF683", "CXCR3", "TNF", "CD40LG", "IL9R", "CD70", "TNFSF8", "LTB", "IL13", "IL6")
PBMC_markers <- c("KLF2", "KLF3", "IL7R", "GNLY", "SELL", "CD52", "CX3CR1", "TCF7", "SOX13", "CCR3", "S1PR1")
post_markers <- c("VCAM1", "KIR2DL4",  "CD27", "KLRD1", "KLRC4", "TSC22D3","TNFSF4")



all_markers <- unique(c(tumor_markers, post_markers, preREP_markers, ExpTIL_markers, PBMC_markers))

#check if all available in the assay
length(all_markers)
all_markers <- all_markers[all_markers %in% rownames(aggro@assays[["RNA"]]$data)]
length(all_markers)

mat<- aggro@assays[["RNA"]]$data[all_markers, ] %>% as.matrix()
mat<- t(scale(t(mat)))

quantile(mat, c(0.05, 0.95))

cluster_anno<- paste(aggro@meta.data$Type, aggro@meta.data$Patient, sep = "_")

## -----------------------------------------------------------------------------------------------------------
#plot heatmap
pdf("Fig/TR_Heatmap_genes_selection.pdf", height = 11, width = 10)
Heatmap(mat, name = "Z score",  
        column_split = factor(aggro@meta.data$Type_post, levels = c("Tumor", "PreREP", "Expanded-TILs", "PBMC-7dpt", "Post"),
                              labels= c("pre-ACT","Intermediate product","TIL product", "PBMC 7dpt", "post-ACT")),
        cluster_columns = F,
        show_column_dend = FALSE,
        cluster_column_slices = TRUE,
        column_title_gp = gpar(fontsize = 14, fontface = "bold"),
        column_gap = unit(0.75, "mm"),
        cluster_rows = F,
        show_row_dend = FALSE,
        col=colorRamp2(c(-1, 0, 2), c("#3C4DA0", "white", "#CD1F2D")),
        row_names_gp = gpar(fontsize = 12),
        column_title_rot = 0,
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("#E64B35B2", "#4DBBD5B2", "#00A087B2", "#3C5488B2", "#F39B7FB2")))),
        show_column_names = F,
        use_raster = TRUE,
        raster_quality = 4,
        row_names_side = "left",
)
dev.off()





