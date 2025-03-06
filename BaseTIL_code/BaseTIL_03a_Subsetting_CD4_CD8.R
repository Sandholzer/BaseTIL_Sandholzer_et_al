## -----------------------------------------------------------------------------------------------------------
#load dependencies

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
})

set.seed(1234567)
setwd("~/BaseTIL_code")
gex.Tcells <- readRDS(file = "saveFiles/BaseTIL_Integrated_Tcells.rds")


## -----------------------------------------------------------------------------------------------------------
#remove doublets
gex.Tcells <- subset(gex.Tcells, chain_pairing %in% c("extra VDJ", "two full chains"), invert=T)


#remove non T cell cluster
gex.Tcells <- subset(gex.Tcells, seurat_clusters %in% c(10, 11, 12), invert=T)


## -----------------------------------------------------------------------------------------------------------
#define CD4 cluster
CD4_cluster <- c(2,6,7)
gex.Tcells$CD4_cells <- FALSE
gex.Tcells$CD4_cells <- gex.Tcells$seurat_clusters %in% CD4_cluster

#Take remaining clusters
# CD8_cluster <- c(1,10,11,12,13,14,15,16)
gex.Tcells$CD8_cells <- FALSE
gex.Tcells$CD8_cells <- !gex.Tcells$seurat_clusters %in% CD4_cluster

gex.CD8 <- subset(gex.Tcells, subset = CD8_cells == TRUE)
gex.CD4 <- subset(gex.Tcells, subset = CD4_cells == TRUE)


## -----------------------------------------------------------------------------------------------------------
#Integration for CD4
gex.CD4[["RNA"]] <- split(gex.CD4[["RNA"]], f = gex.CD4$orig.ident)


gex.CD4 <- gex.CD4 %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) 


#remove unwanted genes from clustering
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

VariableFeatures(gex.CD4) <- gex.features.flt 



gex.CD4  <- ScaleData(gex.CD4, 
                      vars.to.regress = c("S.Score","G2M.Score" ,"ifn_genes1","percent_mito" ), 
                      assay= "RNA")


gex.CD4 <- RunPCA(gex.CD4, features = gex.features.flt, npcs = 30)


gex.CD4 <- IntegrateLayers(
  object = gex.CD4, 
  method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  assay = "RNA",
  verbose = T
)

gex.CD4[["RNA"]] <- JoinLayers(gex.CD4[["RNA"]])


png("Graphs/CD4_Elbow.png", type = "cairo")
ElbowPlot(gex.CD4, ndims = 30) + geom_hline(yintercept=2, linetype="dashed", color = "red")
dev.off()

pdf("Graphs/CD4_pca.heatmap.pdf",width=14,height=18)
DimHeatmap(object = gex.CD4, dims = 1:15, cells = 500, balanced = TRUE,fast = TRUE)
DimHeatmap(object = gex.CD4, dims = 16:30, cells = 500, balanced = TRUE,fast = TRUE)
dev.off()

gex.CD4 <- RunUMAP(gex.CD4, reduction = "harmony", assay = "RNA", dims = 1:20)
gex.CD4 <- FindNeighbors(gex.CD4, dims = 1:20, reduction = "harmony")

resolutions <- seq(0.2, 0.8, by=0.2)

gex.CD4 <- Seurat::FindClusters(gex.CD4, resolution = resolutions, algorithm = 4, method = "igraph")


pdf("Graphs/CD4_Resolutions.pdf", height = 9, width = 16)
DimPlot(gex.CD4, reduction = "umap", raster = F, group.by = paste("RNA_snn_res.", resolutions, sep = ""))
dev.off()


#Define cluster resolution here!
gex.CD4$seurat_clusters <- gex.CD4$RNA_snn_res.0.4


pdf(file = "Graphs/CD4_UMAP.pdf", width = 8, height = 7)
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


pdf(file = "Graphs/CD4_Violin.pdf", width = 14, height = 7)
VlnPlot(gex.CD4, features = "CD4", group.by = "seurat_clusters", assay= "RNA", pt.size = 0)
VlnPlot(gex.CD4, features = "CD8A", group.by = "seurat_clusters", assay= "RNA", pt.size = 0)
VlnPlot(gex.CD4, features = "CD8B", group.by = "seurat_clusters", assay= "RNA", pt.size = 0)
dev.off()


saveRDS(gex.CD4, file = "saveFiles/CD4_integrated.rds", compress = "gzip")


# Clean up environment
ls()
rm(list = ls()[-match(c("gex.CD8"), ls())])
ls()
gc()


## -----------------------------------------------------------------------------------------------------------
##Integration for CD8
gex.CD8[["RNA"]] <- split(gex.CD8[["RNA"]], f = gex.CD8$orig.ident)


gex.CD8 <- gex.CD8 %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) 


gex.features <- VariableFeatures(gex.CD8)
load("data/Zheng_Japrin/exclude.gene.misc.human.v4.RData")


f.feat <- !(gex.features %in% all.gene.ignore.df[["seu.id"]]) &
  !(grepl("^(RP[LS]|Rp[ls])",gex.features,perl=T)) &
  !(grepl("^MT-",gex.features,perl=T)) &
  !(gex.features %in% c("MALAT1", "Malat1", "NEAT1"))
gex.features.flt <- gex.features[f.feat]
print(head(gex.features.flt))
print(paste("Length of anchors:", length(gex.features)) )
print(paste("Length of filtered anchors:", length(gex.features.flt)) )

rm(gex.features)

VariableFeatures(gex.CD8) <- gex.features.flt 



gex.CD8  <- ScaleData(gex.CD8, 
                      vars.to.regress = c("S.Score","G2M.Score" ,"ifn_genes1","percent_mito" ), 
                      assay= "RNA")


gex.CD8 <- RunPCA(gex.CD8, features = gex.features.flt, npcs = 30)


gex.CD8 <- IntegrateLayers(
  object = gex.CD8, 
  method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  assay = "RNA",
  verbose = T
)

gex.CD8[["RNA"]] <- JoinLayers(gex.CD8[["RNA"]])


png("Graphs/CD8_Elbow.png", type = "cairo")
ElbowPlot(gex.CD8, ndims = 30) + geom_hline(yintercept=2, linetype="dashed", color = "red")
dev.off()

pdf("Graphs/CD8_pca.heatmap.pdf",width=14,height=18)
DimHeatmap(object = gex.CD8, dims = 1:15, cells = 500, balanced = TRUE,fast = TRUE)
DimHeatmap(object = gex.CD8, dims = 16:30, cells = 500, balanced = TRUE,fast = TRUE)
dev.off()

gex.CD8 <- RunUMAP(gex.CD8, reduction = "harmony", assay = "RNA", dims = 1:20)
gex.CD8 <- FindNeighbors(gex.CD8, dims = 1:20, reduction = "harmony")

resolutions <- seq(0.2, 0.8, by=0.2)

gex.CD8 <- Seurat::FindClusters(gex.CD8, resolution = resolutions, algorithm = 4, method = "igraph")

pdf("Graphs/CD8_Resolutions.pdf", height = 9, width = 16)
DimPlot(gex.CD8, reduction = "umap", raster = F, group.by = paste("RNA_snn_res.", resolutions, sep = ""))
dev.off()


#Define cluster resolution here!
gex.CD8$seurat_clusters <- gex.CD8$RNA_snn_res.0.4


gc()


pdf(file = "Graphs/CD8_UMAP.pdf", width = 8, height = 7)
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


pdf(file = "Graphs/CD8_Violin.pdf", width = 14, height = 7)
VlnPlot(gex.CD8, features = "CD4", group.by = "seurat_clusters", assay= "RNA", pt.size = 0)
VlnPlot(gex.CD8, features = "CD8A", group.by = "seurat_clusters", assay= "RNA", pt.size = 0)
VlnPlot(gex.CD8, features = "CD8B", group.by = "seurat_clusters", assay= "RNA", pt.size = 0)
dev.off()

saveRDS(gex.CD8, file = "saveFiles/CD8_integrated.rds", compress = "gzip")
