## -----------------------------------------------------------------------------------------------------------
#load dependencies

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
})

set.seed(1234567)
setwd("~/BaseTIL_code")


gex.CD4  <- readRDS(file = "saveFiles/CD4_integrated_2.rds")
gex.CD8  <- readRDS(file = "saveFiles/CD8_integrated_2.rds")


## -----------------------------------------------------------------------------------------------------------
# identify impurity cells based on gene expression 


gex.CD4$isCD8 <-
  ifelse(
    (gex.CD4@assays[["RNA"]]$data["CD4",] < 0.1 &
       gex.CD4@assays[["RNA"]]$data["CD8A",] > 0.1) |
      (gex.CD4@assays[["RNA"]]$data["CD4",] < 0.1 &
         gex.CD4@assays[["RNA"]]$data["CD8B",] > 0.1)
    ,TRUE,FALSE
  )

gex.CD8$isCD4 <-
  ifelse(
    (gex.CD8@assays[["RNA"]]$data["CD4",] > 0.1 &
       gex.CD8@assays[["RNA"]]$data["CD8A",] < 0.1) |
      (gex.CD8@assays[["RNA"]]$data["CD4",] > 0.1 &
         gex.CD8@assays[["RNA"]]$data["CD8B",] < 0.1)
    ,TRUE, FALSE
  )


gex.CD8fromCD4 <- subset(gex.CD4, subset = isCD8 == TRUE)
gex.CD4fromCD8 <- subset(gex.CD8, subset = isCD4 == TRUE)

gex.trueCD4 <- subset(gex.CD4, subset = isCD8 == FALSE)
gex.trueCD8 <- subset(gex.CD8, subset = isCD4 == FALSE)


gex.CD4 <- merge(x= gex.trueCD4, y= gex.CD4fromCD8)
gex.CD4[["RNA"]] <- JoinLayers(gex.CD4[["RNA"]])
rm(gex.trueCD4, gex.CD4fromCD8)
gc()

gex.CD8 <- merge(x= gex.trueCD8, y= gex.CD8fromCD4)
gex.CD8[["RNA"]] <- JoinLayers(gex.CD8[["RNA"]])
rm(gex.trueCD8, gex.CD8fromCD4)
gc()



## -----------------------------------------------------------------------------------------------------------
# Repeat integration for CD4

gex.CD4[["RNA"]] <- split(gex.CD4[["RNA"]], f = gex.CD4$orig.ident)


gex.CD4 <- gex.CD4 %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) 


gex.features <- VariableFeatures(gex.CD4)
load("data/Zheng_Japrin/exclude.gene.misc.human.v4.RData")


f.feat <- !(gex.features %in% all.gene.ignore.df[["seu.id"]]) &
  !(grepl("^RP[LS]",gex.features,perl=T)) &
  !(grepl("^MT-",gex.features,perl=T)) &
  !(gex.features %in% c("MALAT1", "NEAT1"))
gex.features.flt <- gex.features[f.feat]
print(head(gex.features.flt))
print(paste("Length of anchors:", length(gex.features)) )
print(paste("Length of filtered anchors:", length(gex.features.flt)) )

rm(gex.features)

VariableFeatures(gex.CD4) <- gex.features.flt 



gex.CD4  <- ScaleData(gex.CD4, 
                      vars.to.regress = c("S.Score","G2M.Score" ,"ifn_genes1","percent_mito" ), 
                      assay= "RNA")

print("Normalization done.")

gex.CD4 <- RunPCA(gex.CD4, features = gex.features.flt, npcs = 30)

print("start Harmony") 


gex.CD4 <- IntegrateLayers(
  object = gex.CD4, 
  method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  assay = "RNA",
  verbose = T
)

gex.CD4[["RNA"]] <- JoinLayers(gex.CD4[["RNA"]])


png("Graphs/CD4_3_Elbow.png", type = "cairo")
ElbowPlot(gex.CD4, ndims = 30) + geom_hline(yintercept=2, linetype="dashed", color = "red")
dev.off()

pdf("Graphs/CD4_3_pca.heatmap.pdf",width=14,height=18)
DimHeatmap(object = gex.CD4, dims = 1:15, cells = 500, balanced = TRUE,fast = TRUE)
DimHeatmap(object = gex.CD4, dims = 16:30, cells = 500, balanced = TRUE,fast = TRUE)
dev.off()

gex.CD4 <- RunUMAP(gex.CD4, reduction = "harmony", assay = "RNA", dims = 1:20)
gex.CD4 <- FindNeighbors(gex.CD4, dims = 1:20, reduction = "harmony")

resolutions <- seq(0.3, 0.6, by=0.1)

gex.CD4 <- Seurat::FindClusters(gex.CD4, resolution = resolutions, algorithm = 4, method = "igraph")


pdf("Graphs/CD4_3_Resolutions.pdf", height = 9, width = 16)
DimPlot(gex.CD4, reduction = "umap", raster = F, group.by = paste("RNA_snn_res.", resolutions, sep = ""))
dev.off()


#Define cluster resolution here!
gex.CD4$seurat_clusters <- gex.CD4$RNA_snn_res.0.5


gc()


pdf(file = "Graphs/CD4_3_UMAP.pdf", width = 8, height = 7)
DimPlot(gex.CD4, reduction = "umap", group.by = "orig.ident", raster = F)
DimPlot(gex.CD4, reduction = "umap", group.by = "seurat_clusters", label = T, raster = F)
DimPlot(gex.CD4, group.by = "batchV", raster = F)
DimPlot(gex.CD4, group.by = "Patient", raster = F)
DimPlot(gex.CD4, group.by = "Type", raster = F)
DimPlot(gex.CD4, group.by = "Phase", raster = F)
FeaturePlot(gex.CD4, "MKI67", min.cutoff = 0, raster = F)
FeaturePlot(gex.CD4, "ifn_genes1", min.cutoff = 0, raster = F)
FeaturePlot(gex.CD4, "dissoc_genes1", min.cutoff = 0, raster = F)
FeaturePlot(gex.CD4, c("CD3D", "CD4","CD8A",
                       "CD8B","CD14", "CXCL13"), raster = F)
dev.off()


pdf(file = "Graphs/CD4_3_Violin.pdf", width = 14, height = 7)
VlnPlot(gex.CD4, features = "CD4", group.by = "seurat_clusters", assay= "RNA", pt.size = 0)
VlnPlot(gex.CD4, features = "CD8A", group.by = "seurat_clusters", assay= "RNA", pt.size = 0)
VlnPlot(gex.CD4, features = "CD8B", group.by = "seurat_clusters", assay= "RNA", pt.size = 0)
dev.off()


saveRDS(gex.CD4, file = "saveFiles/CD4_integrated_3.rds", compress = "gzip")



ls()
rm(list = ls()[-match(c("gex.CD8"), ls())])
ls()
gc()

## -----------------------------------------------------------------------------------------------------------
# Repeat integration for CD8

gex.CD8[["RNA"]] <- split(gex.CD8[["RNA"]], f = gex.CD8$orig.ident)


gex.CD8 <- gex.CD8 %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) 


gex.features <- VariableFeatures(gex.CD8)
load("data/Zheng_Japrin/exclude.gene.misc.human.v4.RData")


f.feat <- !(gex.features %in% all.gene.ignore.df[["seu.id"]]) &
  !(grepl("^RP[LS]",gex.features,perl=T)) &
  !(grepl("^MT-",gex.features,perl=T)) &
  !(gex.features %in% c("MALAT1", "NEAT1"))
gex.features.flt <- gex.features[f.feat]
print(head(gex.features.flt))
print(paste("Length of anchors:", length(gex.features)) )
print(paste("Length of filtered anchors:", length(gex.features.flt)) )

rm(gex.features)

VariableFeatures(gex.CD8) <- gex.features.flt 



gex.CD8  <- ScaleData(gex.CD8, 
                      vars.to.regress = c("S.Score","G2M.Score" ,"ifn_genes1","percent_mito" ), 
                      assay= "RNA")

print("RNA done.")

gex.CD8 <- RunPCA(gex.CD8, features = gex.features.flt, npcs = 30)

print("start Harmony")


gex.CD8 <- IntegrateLayers(
  object = gex.CD8, 
  method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  assay = "RNA",
  verbose = T
)

gex.CD8[["RNA"]] <- JoinLayers(gex.CD8[["RNA"]])


png("Graphs/CD8_3_Elbow.png", type = "cairo")
ElbowPlot(gex.CD8, ndims = 30) + geom_hline(yintercept=2, linetype="dashed", color = "red")
dev.off()

pdf("Graphs/CD8_3_pca.heatmap.pdf",width=14,height=18)
DimHeatmap(object = gex.CD8, dims = 1:15, cells = 500, balanced = TRUE,fast = TRUE)
DimHeatmap(object = gex.CD8, dims = 16:30, cells = 500, balanced = TRUE,fast = TRUE)
dev.off()

gex.CD8 <- RunUMAP(gex.CD8, reduction = "harmony", assay = "RNA", dims = 1:20)
gex.CD8 <- FindNeighbors(gex.CD8, dims = 1:20, reduction = "harmony")

resolutions <- seq(0.3, 0.8, by=0.1)

gex.CD8 <- Seurat::FindClusters(gex.CD8, resolution = resolutions, algorithm = 4, method = "igraph")


pdf("Graphs/CD8_3_Resolutions.pdf", height = 9, width = 16)
DimPlot(gex.CD8, reduction = "umap", raster = F, group.by = paste("RNA_snn_res.", resolutions, sep = ""))
dev.off()


#Define cluster resolution here!
gex.CD8$seurat_clusters <- gex.CD8$RNA_snn_res.0.6


gc()


pdf(file = "Graphs/CD8_3_UMAP.pdf", width = 8, height = 7)
DimPlot(gex.CD8, reduction = "umap", group.by = "orig.ident", raster = F)
DimPlot(gex.CD8, reduction = "umap", group.by = "seurat_clusters", label = T, raster = F)
DimPlot(gex.CD8, group.by = "batchV", raster = F)
DimPlot(gex.CD8, group.by = "Patient", raster = F)
DimPlot(gex.CD8, group.by = "Type", raster = F)
DimPlot(gex.CD8, group.by = "Phase", raster = F)
FeaturePlot(gex.CD8, "MKI67", min.cutoff = 0, raster = F)
FeaturePlot(gex.CD8, "ifn_genes1", min.cutoff = 0, raster = F)
FeaturePlot(gex.CD8, "dissoc_genes1", min.cutoff = 0, raster = F)
FeaturePlot(gex.CD8, c("CD3D", "CD4","CD8A",
                       "CD8B","CD14", "CXCL13"), raster = F)
dev.off()


pdf(file = "Graphs/CD8_3_Violin.pdf", width = 14, height = 7)
VlnPlot(gex.CD8, features = "CD4", group.by = "seurat_clusters", assay= "RNA", pt.size = 0)
VlnPlot(gex.CD8, features = "CD8A", group.by = "seurat_clusters", assay= "RNA", pt.size = 0)
VlnPlot(gex.CD8, features = "CD8B", group.by = "seurat_clusters", assay= "RNA", pt.size = 0)
dev.off()



saveRDS(gex.CD8, file = "saveFiles/CD8_integrated_3.rds", compress = "gzip")






rm(list = ls())
gc()
