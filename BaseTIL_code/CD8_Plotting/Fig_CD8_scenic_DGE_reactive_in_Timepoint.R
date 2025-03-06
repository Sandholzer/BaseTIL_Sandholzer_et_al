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
  library(dplyr)
  library(ggrepel)
  library(purrr)
  library(AUCell)
})

set.seed(1234567)
setwd("~/BaseTIL_code/CD8_Plotting")

gex.CD8 <- readRDS(file = "../saveFiles/CD8_annotated.rds")

#prevents issues in case of NAs
gex.CD8$reactive_clones <- ifelse(gex.CD8$overall_reactive == TRUE & !is.na(gex.CD8$overall_reactive), TRUE, FALSE)



## -----------------------------------------------------------------------------------------------------------
# perform DGE on regulon activity for baseline tumor (pre-ACT)

arg_type <- "Tumor"


if (arg_type %in% c("post", "Post")) {
  arg_type2 <- c("Rebiopsy1", "Rebiopsy2")
} else{
  arg_type2 <- arg_type
}


#use same subset file as for the scenic run
gex.sub <- subset(gex.CD8, Type %in% arg_type2)

#load corresponding loom file
scenicLoomPath <- paste0("~/BaseTIL/Scenic/pySenic_data/output/pyscenic_",arg_type, "_output.loom")
loom <- open_loom(scenicLoomPath)

# Read information from loom file
regulons_incidMat <- get_regulons(loom, column.attr.name='Regulons')
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
regulonsAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)

AUCmat <- AUCell::getAUC(regulonsAUC)
rownames(AUCmat) <- gsub("\\(\\+\\)", "", rownames(AUCmat))

# add regulon activity into new assay
gex.sub@assays[["AUC"]] <- CreateAssayObject(data = AUCmat)


# only consider sampels with sufficiently high numbers of identified tumor-reactive T cells
tum.CD8 <- subset(gex.sub, subset= Type %in% arg_type2 & Patient %in% c("UPN001", "UPN003", "UPN006","UPN008", "UPN009", "UPN011")) #,"UPN002"
tum.CD8 <- subset(tum.CD8, Sample_Name %in% c("UPN008 Rebiopsy1", "UPN001 Rebiopsy2"), invert = T)


# scale new regulon assay
DefaultAssay(tum.CD8) <- "AUC"
tum.CD8 <- ScaleData(tum.CD8, assay = 'AUC', features = rownames(AUCmat))
Idents(tum.CD8) <- tum.CD8$reactive_clones

# perform Differential expression analysis on regulon assay
markers <- FindMarkers(tum.CD8, ident.1 = "TRUE",
                       ident.2 = "FALSE", logfc.threshold = 0, slot = "scale.data", fc.name = "avg_log2FC", fc.slot="scale.data") #, latent.vars ="Patient" , test.use = "MAST", min.pct = 0, min.diff.pct = 0

markers$gene <- rownames(markers)


# define annotation
selection <- c("STAT3", "HIVEP1", "TRPS1", "EOMES", "ETV1", "IRF9", 
               "GATA3", "PRDM1", "IRF1", "KLF3", "NFKB2", "ETV7", "BATF", "TBX21", "KLF13", 
               "IRF4", "IRF2", "NFAT5", "HOXD1", "IKZF3", "MAF", "CRX", 
               "MYC", "BACH2", "ETV3", "TCF7", "FOS", "LEF1", "TCF7L2", "FOXP1", "ZNF300", "KLF11", "JUNB")
markers$selection <- ifelse(rownames(markers) %in% selection, rownames(markers), NA)


#define color scale
keyvals <- ifelse(
  markers$avg_log2FC < -0.25, '#7E8CBE',
  ifelse(markers$avg_log2FC > 0.25, '#F8766D',
         'grey80'))
keyvals[is.na(keyvals)] <- 'grey80'
names(keyvals)[keyvals == '#7E8CBE'] <- 'Downregulated'
names(keyvals)[keyvals == 'grey80'] <- 'NS'
names(keyvals)[keyvals == '#F8766D'] <- 'Upregulated'


# plotting
pdf("Graphs/CD8_TR_Tumor_vulc.pdf", height = 5, width = 5)
plt1 <- EnhancedVolcano::EnhancedVolcano(
  markers,
  x = "avg_log2FC",
  y = "p_val",
  rownames(markers),
  FCcutoff = 0.25,
  pCutoff = 10e-50,
  xlim = c(-1.2, 1.2),
  ylim = c(0, 350),
  xlab = "single-cell log2FC",
  selectLab = NA,
  parseLabels = T,
  drawConnectors = T, 
  subtitle = NULL, caption = NULL, 
  title = NULL, 
  pointSize = 1.0,
  labSize = 3.0,
  #col=c( "#D6D6D6", "#D6D6D6", "#D6D6D6","#CD2626"),
  colCustom = keyvals,
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
plt1
dev.off()






## -----------------------------------------------------------------------------------------------------------
# perform DGE on regulon activity for preREP (intermediate product)

arg_type <- "PreREP"

if (arg_type %in% c("post", "Post")) {
  arg_type2 <- c("Rebiopsy1", "Rebiopsy2")
} else{
  arg_type2 <- arg_type
}

#use same subset file as for the scenic run
gex.sub <- subset(gex.CD8, Type %in% arg_type2)

#load corresponding loom file
scenicLoomPath <- paste0("~/BaseTIL/Scenic/pySenic_data/output/pyscenic_",arg_type, "_output.loom")
loom <- open_loom(scenicLoomPath)

# Read information from loom file
regulons_incidMat <- get_regulons(loom, column.attr.name='Regulons')
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
regulonsAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)
AUCmat <- AUCell::getAUC(regulonsAUC)
rownames(AUCmat) <- gsub("\\(\\+\\)", "", rownames(AUCmat))

gex.sub@assays[["AUC"]] <- CreateAssayObject(data = AUCmat)

# only consider sampels with sufficiently high numbers of identified tumor-reactive T cells
tum.CD8 <- subset(gex.sub, subset= Type %in% arg_type2 & Patient %in% c("UPN001", "UPN003", "UPN008", "UPN009", "UPN011")) #,"UPN002"
tum.CD8 <- subset(tum.CD8, Sample_Name %in% c("UPN008 Rebiopsy1", "UPN001 Rebiopsy2"), invert = T)

# scale new regulon assay
DefaultAssay(tum.CD8) <- "AUC"
tum.CD8 <- ScaleData(tum.CD8, assay = 'AUC', features = rownames(AUCmat))
Idents(tum.CD8) <- tum.CD8$reactive_clones

# perform Differential expression analysis on regulon assay
markers <- FindMarkers(tum.CD8, ident.1 = "TRUE",
                       ident.2 = "FALSE", logfc.threshold = 0, slot = "scale.data", fc.name = "avg_log2FC", fc.slot="scale.data") #, latent.vars ="Patient" , test.use = "MAST", min.pct = 0, min.diff.pct = 0
markers$gene <- rownames(markers)

# define annotation
selection <- c("ETV1" ,"IRF5","FOS", "JUNB", "REL", "FOSB", "JUND", "CEBPB",
               "KLF9", "NFKB1", "PRDM1", "MAF", "KLF6",  "SOX5", "RELB",
               "CLOCK", "STAT2", "STAT1", "PURA", "NFE2L1", "RUNX3", "ETV6", "TCF7L1", "FOXO3", "ETV3", "RUNX1", "STAT6", "NR2C2")

markers$selection <- ifelse(rownames(markers) %in% selection, rownames(markers), NA)



#define color scale
keyvals <- ifelse(
  markers$avg_log2FC < -0.25, '#7E8CBE',
  ifelse(markers$avg_log2FC > 0.25, '#F8766D',
         'grey80'))
keyvals[is.na(keyvals)] <- 'grey80'
names(keyvals)[keyvals == '#7E8CBE'] <- 'Downregulated'
names(keyvals)[keyvals == 'grey80'] <- 'NS'
names(keyvals)[keyvals == '#F8766D'] <- 'Upregulated'




# plotting
pdf(paste0("Graphs/CD8_TR_",arg_type, "_vulc.pdf"), height = 5, width = 5)
plt2 <- EnhancedVolcano::EnhancedVolcano(
  markers,
  x = "avg_log2FC",
  y = "p_val",
  rownames(markers),
  FCcutoff = 0.25,
  pCutoff = 10e-50,
  xlim = c(-1.2, 1.2),
  ylim = c(0, 350),
  xlab = "single-cell log2FC",
  selectLab = NA,
  parseLabels = T,
  drawConnectors = T, 
  subtitle = NULL, caption = NULL, 
  title = NULL, 
  pointSize = 1.0,
  labSize = 3.0,
  #col=c( "#D6D6D6", "#D6D6D6", "#D6D6D6","#CD2626"),
  colCustom = keyvals,
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
plt2
dev.off()






## -----------------------------------------------------------------------------------------------------------
# perform DGE on regulon activity for Expanded TILs (TIL product)

arg_type <- "Expanded_TILs"

if (arg_type %in% c("post", "Post")) {
  arg_type2 <- c("Rebiopsy1", "Rebiopsy2")
} else{
  arg_type2 <- arg_type
}

#use same subset file as for the scenic run
gex.sub <- subset(gex.CD8, Type %in% arg_type2)

#load corresponding loom file
scenicLoomPath <- paste0("~/BaseTIL/Scenic/pySenic_data/output/pyscenic_",arg_type, "_output.loom")
loom <- open_loom(scenicLoomPath)

# Read information from loom file
regulons_incidMat <- get_regulons(loom, column.attr.name='Regulons')
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
regulonsAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)
AUCmat <- AUCell::getAUC(regulonsAUC)
rownames(AUCmat) <- gsub("\\(\\+\\)", "", rownames(AUCmat))

gex.sub@assays[["AUC"]] <- CreateAssayObject(data = AUCmat)

# only consider sampels with sufficiently high numbers of identified tumor-reactive T cells
tum.CD8 <- subset(gex.sub, subset= Type %in% arg_type2 & Patient %in% c("UPN001", "UPN003", "UPN008", "UPN009", "UPN011")) #,"UPN002"
tum.CD8 <- subset(tum.CD8, Sample_Name %in% c("UPN008 Rebiopsy1", "UPN001 Rebiopsy2"), invert = T)

# scale new regulon assay
DefaultAssay(tum.CD8) <- "AUC"
tum.CD8 <- ScaleData(tum.CD8, assay = 'AUC', features = rownames(AUCmat))
Idents(tum.CD8) <- tum.CD8$reactive_clones

# perform Differential expression analysis on regulon assay
markers <- FindMarkers(tum.CD8, ident.1 = "TRUE",
                       ident.2 = "FALSE", logfc.threshold = 0, slot = "scale.data", fc.name = "avg_log2FC", fc.slot="scale.data") #, latent.vars ="Patient" , test.use = "MAST", min.pct = 0, min.diff.pct = 0
markers$gene <- rownames(markers)

# define annotation
selection <- c("ETV1" ,"IRF5","BCL3", "ETV7", "JUND", "CEBPB",
               "IRF5", "IRF8", "TBX21", "EOMES", "RORC", "NFKB2", "PPARG", "IRF7",
               "SKI", "MYB", "HOXB5", "RUNX3", "NFKB1", "IKZF1", "ETV3", "RUNX1")

markers$selection <- ifelse(rownames(markers) %in% selection, rownames(markers), NA)


#define color scale
keyvals <- ifelse(
  markers$avg_log2FC < -0.25, '#7E8CBE',
  ifelse(markers$avg_log2FC > 0.25, '#F8766D',
         'grey80'))
keyvals[is.na(keyvals)] <- 'grey80'
names(keyvals)[keyvals == '#7E8CBE'] <- 'Downregulated'
names(keyvals)[keyvals == 'grey80'] <- 'NS'
names(keyvals)[keyvals == '#F8766D'] <- 'Upregulated'




# plotting
pdf(paste0("Graphs/CD8_TR_",arg_type, "_vulc.pdf"), height = 5, width = 5)
plt3 <- EnhancedVolcano::EnhancedVolcano(
  markers,
  x = "avg_log2FC",
  y = "p_val",
  rownames(markers),
  FCcutoff = 0.25,
  pCutoff = 10e-50,
  xlim = c(-1.2, 1.2),
  ylim = c(0, 350),
  xlab = "single-cell log2FC",
  selectLab = NA,
  parseLabels = T,
  drawConnectors = T, 
  subtitle = NULL, caption = NULL, 
  title = NULL, 
  pointSize = 1.0,
  labSize = 3.0,
  #col=c( "#D6D6D6", "#D6D6D6", "#D6D6D6","#CD2626"),
  colCustom = keyvals,
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
                       show.legend = FALSE, max.overlaps = 20, point.padding = 0, label.padding = 0.15)
plt3
dev.off()








## -----------------------------------------------------------------------------------------------------------
# perform DGE on regulon activity for PBMC 7dpt

arg_type <- "PBMC_7dpt"

if (arg_type %in% c("post", "Post")) {
  arg_type2 <- c("Rebiopsy1", "Rebiopsy2")
} else{
  arg_type2 <- arg_type
}

#use same subset file as for the scenic run
gex.sub <- subset(gex.CD8, Type %in% arg_type2)

#load corresponding loom file
scenicLoomPath <- paste0("~/BaseTIL/Scenic/pySenic_data/output/pyscenic_",arg_type, "_output.loom")
loom <- open_loom(scenicLoomPath)

# Read information from loom file
regulons_incidMat <- get_regulons(loom, column.attr.name='Regulons')
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
regulonsAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)
AUCmat <- AUCell::getAUC(regulonsAUC)
rownames(AUCmat) <- gsub("\\(\\+\\)", "", rownames(AUCmat))

gex.sub@assays[["AUC"]] <- CreateAssayObject(data = AUCmat)

# only consider sampels with sufficiently high numbers of identified tumor-reactive T cells
tum.CD8 <- subset(gex.sub, subset= Type %in% arg_type2 & Patient %in% c("UPN001", "UPN003", "UPN008", "UPN009", "UPN011")) #,"UPN002"
tum.CD8 <- subset(tum.CD8, Sample_Name %in% c("UPN008 Rebiopsy1", "UPN001 Rebiopsy2"), invert = T)

# scale new regulon assay
DefaultAssay(tum.CD8) <- "AUC"
tum.CD8 <- ScaleData(tum.CD8, assay = 'AUC', features = rownames(AUCmat))
Idents(tum.CD8) <- tum.CD8$reactive_clones

# perform Differential expression analysis on regulon assay
markers <- FindMarkers(tum.CD8, ident.1 = "TRUE",
                       ident.2 = "FALSE", logfc.threshold = 0, slot = "scale.data", fc.name = "avg_log2FC", fc.slot="scale.data") #, latent.vars ="Patient" , test.use = "MAST", min.pct = 0, min.diff.pct = 0
markers$gene <- rownames(markers)

# define annotation
selection <- c("ETV7" ,"HTATIP2","BACH2", "IRF7",  "BACH2", "HTATIP2",
                "KLF2", "ETV6", "IKZF1", "MECP2")

markers$selection <- ifelse(rownames(markers) %in% selection, rownames(markers), NA)


#define color scale
keyvals <- ifelse(
  markers$avg_log2FC < -0.25, '#7E8CBE',
  ifelse(markers$avg_log2FC > 0.25, '#F8766D',
         'grey80'))
keyvals[is.na(keyvals)] <- 'grey80'
names(keyvals)[keyvals == '#7E8CBE'] <- 'Downregulated'
names(keyvals)[keyvals == 'grey80'] <- 'NS'
names(keyvals)[keyvals == '#F8766D'] <- 'Upregulated'


# plotting
pdf(paste0("Graphs/CD8_TR_",arg_type, "_vulc.pdf"), height = 5, width = 5)
plt4 <- EnhancedVolcano::EnhancedVolcano(
  markers,
  x = "avg_log2FC",
  y = "p_val",
  rownames(markers),
  FCcutoff = 0.25,
  pCutoff = 10e-50,
  xlim = c(-1.2, 1.2),
  ylim = c(0, 350),
  xlab = "single-cell log2FC",
  selectLab = NA,
  parseLabels = T,
  drawConnectors = T, 
  subtitle = NULL, caption = NULL, 
  title = NULL, 
  pointSize = 1.0,
  labSize = 3.0,
  #col=c( "#D6D6D6", "#D6D6D6", "#D6D6D6","#CD2626"),
  colCustom = keyvals,
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
plt4
dev.off()





## -----------------------------------------------------------------------------------------------------------
# perform DGE on regulon activity for post treatment lesions (post-ACT)

arg_type <- "Post"

if (arg_type %in% c("post", "Post")) {
  arg_type2 <- c("Rebiopsy1", "Rebiopsy2")
} else{
  arg_type2 <- arg_type
}

#use same subset file as for the scenic run
gex.sub <- subset(gex.CD8, Type %in% arg_type2)

#load corresponding loom file
scenicLoomPath <- paste0("~/BaseTIL/Scenic/pySenic_data/output/pyscenic_Rebiopsy1_output.loom")
loom <- open_loom(scenicLoomPath)

# Read information from loom file
regulons_incidMat <- get_regulons(loom, column.attr.name='Regulons')
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
regulonsAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)
AUCmat <- AUCell::getAUC(regulonsAUC)
rownames(AUCmat) <- gsub("\\(\\+\\)", "", rownames(AUCmat))

AUCmatR1 <- AUCmat

#load corresponding loom file
scenicLoomPath <- paste0("~/BaseTIL/Scenic/pySenic_data/output/pyscenic_Rebiopsy2_output.loom")
loom <- open_loom(scenicLoomPath)

# Read information from loom file
regulons_incidMat <- get_regulons(loom, column.attr.name='Regulons')
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
regulonsAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)
AUCmat <- AUCell::getAUC(regulonsAUC)
rownames(AUCmat) <- gsub("\\(\\+\\)", "", rownames(AUCmat))

AUCmatR2 <- AUCmat

data_frame_merge <- merge(as.data.frame(AUCmatR1), as.data.frame(AUCmatR2), 
                          by = 'row.names', all = TRUE) 

rownames(data_frame_merge) <- data_frame_merge$Row.names
data_frame_merge$Row.names <- NULL

data_frame_merge[is.na(data_frame_merge)] <- 0


gex.sub@assays[["AUC"]] <- CreateAssayObject(data = as.matrix(data_frame_merge), key = "AUC_")

# only consider sampels with sufficiently high numbers of identified tumor-reactive T cells
tum.CD8 <- subset(gex.sub, subset= Type %in% arg_type2 & Patient %in% c("UPN001","UPN003", "UPN008", "UPN009", "UPN011")) #,"UPN002"
tum.CD8 <- subset(tum.CD8, Sample_Name %in% c("UPN008 Rebiopsy1", "UPN001 Rebiopsy2"), invert = T)

# scale new regulon assay
DefaultAssay(tum.CD8) <- "AUC"
tum.CD8 <- ScaleData(tum.CD8, assay = 'AUC', features = rownames(AUCmat))
Idents(tum.CD8) <- tum.CD8$reactive_clones

# perform Differential expression analysis on regulon assay
markers <- FindMarkers(tum.CD8, ident.1 = "TRUE",
                       ident.2 = "FALSE", logfc.threshold = 0, slot = "scale.data", fc.name = "avg_log2FC", fc.slot="scale.data") #, latent.vars ="Patient" , test.use = "MAST", min.pct = 0, min.diff.pct = 0
markers$gene <- rownames(markers)

# define annotation
selection <- c("RFX5" ,"ETV1", "GFI1","VDR", "NFATC1", "NR5A2", "TRPS1", "IRF9", "EOMES", "MYB", "PRDM1",
               "RARA", "HMGA1", "KLF2", "NFKB1", "RUNX3", "LEF1", "NFATC2", "JUNB", "KLF3", "RORC", "REL", "KLF13")

markers$selection <- ifelse(rownames(markers) %in% selection, rownames(markers), NA)


#define color scale
keyvals <- ifelse(
  markers$avg_log2FC < -0.25, '#7E8CBE',
  ifelse(markers$avg_log2FC > 0.25, '#F8766D',
         'grey80'))
keyvals[is.na(keyvals)] <- 'grey80'
names(keyvals)[keyvals == '#7E8CBE'] <- 'Downregulated'
names(keyvals)[keyvals == 'grey80'] <- 'NS'
names(keyvals)[keyvals == '#F8766D'] <- 'Upregulated'


# plotting
pdf(paste0("Graphs/CD8_TR_",arg_type, "_vulc.pdf"), height = 5, width = 5)
plt5 <- EnhancedVolcano::EnhancedVolcano(
  markers,
  x = "avg_log2FC",
  y = "p_val",
  rownames(markers),
  FCcutoff = 0.25,
  pCutoff = 10e-50,
  xlab = "single-cell log2FC",
  xlim = c(-1.2, 1.2),
  ylim = c(0, 350),
  selectLab = NA,
  parseLabels = T,
  drawConnectors = T, 
  subtitle = NULL, caption = NULL, 
  title = NULL, 
  pointSize = 1.0,
  labSize = 3.0,
  #col=c( "#D6D6D6", "#D6D6D6", "#D6D6D6","#CD2626"),
  colCustom = keyvals,
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

plt5
dev.off()



#plot all together
pdf(paste0("Graphs/CD8_TR_combined_vulc.pdf"), height = 5, width = 25)
ggarrange(plt1,plt2,plt3,plt4,plt5, ncol = 5, nrow = 1)
dev.off()
