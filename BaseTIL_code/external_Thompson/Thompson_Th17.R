## -----------------------------------------------------------------------------------------------------------
#load dependencies
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(ggeasy)
  library(Seurat)  
  library(UCell)
  library(harmony)
  library(readr)
  library(ggeasy)
  library(ggpubr)
  library(viridis)
  library(purrr)
  library(tidyr)
  library(dplyr)
  library(stringr)
  library(ComplexHeatmap)
  library(circlize)
  library(DESeq2)
  library(ggrepel)
  library(tibble)
  
})

set.seed(12345678)

setwd("~/BaseTIL_code/external_Thompson")

## -----------------------------------------------------------------------------------------------------------
#load the raw data


counts <- read.csv("~/BaseTIL_code/external_Thompson/data/GSE218004_rna_counts.csv", row.names = "Gene")

df_DP <- counts %>% select(contains("_DP"))
df_CD4 <- counts %>% select(contains("_CD4"))
df_CD8 <- counts %>% select(contains("_CD8"))

colnames(df_DP) <- gsub("_DP", x=colnames(df_DP), replacement = "")
colnames(df_CD4) <- gsub("_CD4", x=colnames(df_CD4), replacement = "")
colnames(df_CD8) <- gsub("_CD8", x=colnames(df_CD8), replacement = "")



meta <- read.csv("~/BaseTIL_code/external_Thompson/data/Thompson_meta.csv", row.names = "Patient")

meta$RNR <- ifelse(meta$Response %in% c("PD"), "NRs", "Rs") 

## -----------------------------------------------------------------------------------------------------------
# Plot Heatmap for interesting genes form signature for CD4

selection <- c("RORC","KLRB1", "ADAM12","AQP3", "CCR6",  "IL23R", "IL17RE",   "CCL20", "IL26","IL22", "IL21","IL4I1") 

length(selection)

mat<- as.matrix(df_CD4)
mat<- t(scale(t(mat)))
mat <- mat[selection,]


selection <- selection[selection %in% rownames(mat)]
length(selection)

quantile(mat, c(0.05, 0.95))

#define order of columns
col.order <- c("Pt15","Pt20","Pt22","Pt07","Pt11", "Pt08", "Pt18", "Pt23", "Pt21", "Pt24", "Pt09", "Pt16", "Pt14")
mat<- mat[ , col.order]

pdf("Graphs/Heatmap_genes_Thompson_CD4.pdf", height = 4, width = 6)

Heatmap(mat, name = "Z score",  
        column_split = factor(meta[col.order,]$RNR, levels = c("NRs", "Rs")),
        cluster_columns = F,
        show_column_dend = FALSE,
        cluster_column_slices = TRUE,
        column_title_gp = gpar(fontsize = 12),
        column_gap = unit(0.75, "mm"),
        cluster_rows = F,
        show_row_dend = FALSE,
        col=colorRamp2(c(-0.9, 0, 2.0), c("#3C4DA0", "white", "#CD1F2D")),
        row_names_gp = gpar(fontsize = 12),
        column_title_rot = 0,
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("#CD2626", "#009ACD")))),
        show_column_names = F,
        use_raster = TRUE,
        raster_quality = 4,
        row_names_side = "left",

)

dev.off()


## -----------------------------------------------------------------------------------------------------------
# Differential gene expression for CD4

colData_meta <- meta[rownames(meta) %in% colnames(df_CD4), ]
colData_meta$Patient <- rownames(colData_meta)


#create DESeq2 object
dds <- DESeqDataSetFromMatrix(df_CD4, 
                              colData = colData_meta, 
                              design = ~ RNR + LDH )
rld <- rlog(dds, blind=TRUE)

dds <- DESeq(dds)
plotDispEsts(dds)

#filtering
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]


resultsNames(dds)

# Generate results object
res <- results(dds, 
               name = "RNR_R_vs_NR",
               alpha = 0.1)


res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  arrange(padj)


#filter unwanted genes
res_table_thres <- res_tbl[!is.na(res_tbl$pvalue), ] 

load("~/BaseTIL_code/data/Zheng_Japrin/exclude.gene.misc.human.v4.RData")
f.feat <- !(res_table_thres$gene %in% all.gene.ignore.df[["seu.id"]]) &
  !(grepl("^RP[LS]",res_table_thres$gene,perl=T)) &
  !(grepl("^MT-|^MTRN|^MTAT|^MTND|^MRP",res_table_thres$gene,perl=T)) &
  !(grepl("^RP[0-9]+-",res_table_thres$gene,perl=T)) &
  !(res_table_thres$gene %in% c("MALAT1", 'XIST', 'NEAT1'))
markers_bulk.flt <- res_table_thres[f.feat,]
print(head(markers_bulk.flt))
print(paste("Length of anchors:", nrow(res_table_thres)) )
print(paste("Length of filtered anchors:", nrow(markers_bulk.flt)) )


show.genes_pos <- c("CCR7", "TNFRSF9", "HLA-DRB1")
show.genes_neg <- c("RORC", "IL26", "IL23R", "IL17RE", "ADAM12") 

markers_bulk.flt$selection <- ifelse(markers_bulk.flt$gene %in% c(show.genes_pos, show.genes_neg),markers_bulk.flt$gene, NA)


#create the plot
pdf(paste0("Graphs/DGE_response_Thompson_CD4.pdf"), width = 5, height = 5)
EnhancedVolcano::EnhancedVolcano(
  markers_bulk.flt,
  x = "log2FoldChange",
  y = "pvalue",
  rownames(markers_bulk.flt),
  FCcutoff = 1,
  pCutoff = 0.05,
  ylim = c(0, 7),
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
  )  & geom_label_repel(aes(label = selection),
                        parse = T,
                        size = 3,
                        box.padding = 0.25,
                        segment.color = "black",min.segment.length = 0,
                        show.legend = FALSE, max.overlaps = 15, point.padding = 0, label.padding = 0.15)


dev.off()






## -----------------------------------------------------------------------------------------------------------
# Plot Heatmap for interesting genes form signature for CD8


selection <- c("RORC", "IL26", "IL17RE", "IL4I1","IL23R",  "AQP3","CCL20" ,"IL22", "IL21","LTK"
               , "KLRB1", "ADAM12", "CCR6", "LTB", "TNF", "LTA")

length(selection)

mat<- as.matrix(df_CD8)
mat<- t(scale(t(mat)))
mat <- mat[selection,]


selection <- selection[selection %in% rownames(mat)]
length(selection)


col.order <- c("Pt15","Pt20","Pt22","Pt07","Pt11", "Pt08", "Pt18", "Pt23", "Pt21", "Pt24", "Pt09", "Pt16", "Pt14")
mat<- mat[ , col.order]

pdf("Graphs/Heatmap_genes_Thompson_CD8.pdf", height = 4, width = 6)

Heatmap(mat, name = "Z score",  
        column_split = factor(meta[rownames(meta) %in% colnames(mat), ]$RNR, levels = c("NRs", "Rs")),
        cluster_columns = F,
        show_column_dend = FALSE,
        cluster_column_slices = TRUE,
        column_title_gp = gpar(fontsize = 12),
        column_gap = unit(0.75, "mm"),
        cluster_rows = T,
        show_row_dend = FALSE,
        col=colorRamp2(c(-0.9, 0, 2.0), c("#3C4DA0", "white", "#CD1F2D")),
        row_names_gp = gpar(fontsize = 12),
        column_title_rot = 0,
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("#CD2626", "#009ACD")))),
        show_column_names = T,
        use_raster = TRUE,
        raster_quality = 4,
        row_names_side = "left",

)
dev.off()




## -----------------------------------------------------------------------------------------------------------
# Differential gene expression for CD8

colData_meta <- meta[rownames(meta) %in% colnames(df_CD8), ]
colData_meta$Patient <- rownames(colData_meta)


#create DESeq2 object
dds <- DESeqDataSetFromMatrix(df_CD8, 
                              colData = colData_meta, 
                              design = ~ RNR + LDH )

rld <- rlog(dds, blind=TRUE)

dds <- DESeq(dds)
plotDispEsts(dds)

#filtering
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]


resultsNames(dds)

# Generate results object
res <- results(dds, 
               name = "RNR_R_vs_NR",
               alpha = 0.1)


res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  arrange(padj)


#filter unwanted genes
res_table_thres <- res_tbl[!is.na(res_tbl$pvalue), ] 

load("~/BaseTIL_code/data/Zheng_Japrin/exclude.gene.misc.human.v4.RData")
f.feat <- !(res_table_thres$gene %in% all.gene.ignore.df[["seu.id"]]) &
  !(grepl("^RP[LS]",res_table_thres$gene,perl=T)) &
  !(grepl("^MT-|^MTRN|^MTAT|^MTND|^MRP",res_table_thres$gene,perl=T)) &
  !(grepl("^RP[0-9]+-",res_table_thres$gene,perl=T)) &
  !(res_table_thres$gene %in% c("MALAT1", 'XIST', 'NEAT1'))
markers_bulk.flt <- res_table_thres[f.feat,]
print(head(markers_bulk.flt))
print(paste("Length of anchors:", nrow(res_table_thres)) )
print(paste("Length of filtered anchors:", nrow(markers_bulk.flt)) )


show.genes_pos <- c("CCR7", "TNFRSF9", "HLA-DRB1")
show.genes_neg <- c("RORC", "IL26", "IL23R", "IL17RE", "ADAM12") #"KLRB1"

markers_bulk.flt$selection <- ifelse(markers_bulk.flt$gene %in% c(show.genes_pos, show.genes_neg),markers_bulk.flt$gene, NA)

pdf(paste0("Graphs/DGE_response_Thompson_CD8.pdf"), width = 5, height = 5)
EnhancedVolcano::EnhancedVolcano(
  markers_bulk.flt,
  x = "log2FoldChange",
  y = "pvalue",
  rownames(markers_bulk.flt),
  FCcutoff = 1,
  pCutoff = 0.05,
  ylim = c(0, 7),
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
  )  & geom_label_repel(aes(label = selection),
                        parse = T,
                        size = 3,
                        box.padding = 0.25,
                        segment.color = "black",min.segment.length = 0,
                        show.legend = FALSE, max.overlaps = 15, point.padding = 0, label.padding = 0.15)


dev.off()


