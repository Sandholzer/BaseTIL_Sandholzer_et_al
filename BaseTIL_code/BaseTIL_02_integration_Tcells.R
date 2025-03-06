## -----------------------------------------------------------------------------------------------------------
#load dependencies

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
})

set.seed(1234567)
setwd("~/BaseTIL_code")
gex.Tcells <- readRDS(file = "saveFiles/BaseTIL_QC.rds")
## -----------------------------------------------------------------------------------------------------------
# Normalizaiton, Feature selection and Scaling

gex.Tcells <- gex.Tcells %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

gex.Tcells  <- ScaleData(gex.Tcells)

gex.Tcells[["RNA"]] <- JoinLayers(gex.Tcells[["RNA"]])
DefaultAssay(gex.Tcells) <- "RNA"


## -----------------------------------------------------------------------------------------------------------
#Calculating cell cycle, dissociation and interferon response contribution
Seurat::cc.genes.updated.2019
s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes

gex.Tcells <- Seurat::CellCycleScoring(gex.Tcells,
                                     s.features = s.genes,
                                     g2m.features = g2m.genes, search = T)


dissoc_genes <- readLines("data/dissocation_genes.txt")
gex.Tcells <- AddModuleScore(gex.Tcells,
                           features = list(dissoc_genes),
                           name = "dissoc_genes", search = T)

ifn_genes <- readLines("data/interferon_responsive_genes.txt")
gex.Tcells <- AddModuleScore(gex.Tcells,
                           features = list(ifn_genes),
                           name = "ifn_genes", search = T)



gex.Tcells[["RNA"]] <- split(gex.Tcells[["RNA"]], f = gex.Tcells$orig.ident)
DefaultAssay(gex.Tcells) <- "RNA"


## -----------------------------------------------------------------------------------------------------------
## Remove unwanted genes like TCR genes, ribo genes and histone genes from variable features
gex.features <- VariableFeatures(gex.Tcells, assay = "RNA")
load("data/Zheng_Japrin/exclude.gene.misc.human.v4.RData")

print(paste("Length of anchors:", length(gex.features)) )

f.feat <- !(gex.features %in% all.gene.ignore.df[["seu.id"]]) &
  !(grepl("^RP[LS]",gex.features,perl=T)) 
gex.features <- gex.features[f.feat]
print(gex.features)
print(paste("Length of filtered anchors:", length(gex.features)) )

VariableFeatures(gex.Tcells, assay = "RNA") <- gex.features 


## -----------------------------------------------------------------------------------------------------------
#regress out disturbing factors

gex.Tcells  <- ScaleData(gex.Tcells, 
                       vars.to.regress = c("S.Score","G2M.Score" ,"ifn_genes1","percent_mito" ), 
                       assay= "RNA")



## -----------------------------------------------------------------------------------------------------------
#Harmony batch correction

gex.Tcells <- RunPCA(gex.Tcells, npcs = 30, assay = "RNA") 

gex.Tcells <- IntegrateLayers(
  object = gex.Tcells, 
  method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  assay = "RNA",
  verbose = T
)
gex.Tcells[["RNA"]] <- JoinLayers(gex.Tcells[["RNA"]])

## -----------------------------------------------------------------------------------------------------------
#Evaluate dimensions for clustering

png("Graphs/Harm_Elbow.png", type = "cairo")
ElbowPlot(gex.Tcells, ndims = 30) + geom_hline(yintercept=2, linetype="dashed", color = "red")
dev.off()

pdf("Graphs/Harm_pca.heatmap.pdf",width=14,height=18)
DimHeatmap(object = gex.Tcells, dims = 1:15, cells = 500, balanced = TRUE,fast = TRUE)
DimHeatmap(object = gex.Tcells, dims = 16:30, cells = 500, balanced = TRUE,fast = TRUE)
dev.off()

## -----------------------------------------------------------------------------------------------------------
# Dimensional reduction and Leiden clustering

gex.Tcells <- RunUMAP(gex.Tcells, reduction = "harmony", dims = 1:20, min.dist = 0.3, n.neighbors = 30L)
gex.Tcells <- FindNeighbors(gex.Tcells, reduction = "harmony", dims = 1:20)

resolutions <- seq(0.2, 0.8, by=0.2)
gex.Tcells <- FindClusters(gex.Tcells, resolution = resolutions, algorithm = 4, method = "igraph")


pdf("Graphs/Harm_Resolutions.pdf", height = 9, width = 16)
DimPlot(gex.Tcells, reduction = "umap", raster = F, group.by = paste("RNA_snn_res.", resolutions, sep = ""))
dev.off()

## -----------------------------------------------------------------------------------------------------------
#Define cluster resolution here!
gex.Tcells$seurat_clusters <- gex.Tcells$RNA_snn_res.0.4


## -----------------------------------------------------------------------------------------------------------
# Plotting some parameters

pdf(file = "Graphs/Integrated_UMAP.pdf", width = 8, height = 7)
DimPlot(gex.Tcells, reduction = "umap", group.by = "orig.ident", raster = F)
DimPlot(gex.Tcells, reduction = "umap", group.by = "seurat_clusters", label = T, raster = F)
DimPlot(gex.Tcells, group.by = "batchV", raster = F)
DimPlot(gex.Tcells, group.by = "Patient", raster = F)
DimPlot(gex.Tcells, group.by = "Type", raster = F)
DimPlot(gex.Tcells, group.by = "Phase", raster = F)
FeaturePlot(gex.Tcells, "MKI67", min.cutoff = 0, raster = F)
FeaturePlot(gex.Tcells, "ifn_genes1", min.cutoff = 0, raster = F)
FeaturePlot(gex.Tcells, "dissoc_genes1", min.cutoff = 0, raster = F)
FeaturePlot(gex.Tcells, c("CD3D", "CD4","CD8A",
                        "CD8B","CD14", "PDCD1"), raster = F)
dev.off()


pdf(file = "Graphs/Integrated_Violin.pdf", width = 14, height = 7)
VlnPlot(gex.Tcells, features = "CD3D", group.by = "seurat_clusters", assay= "RNA", pt.size = 0)
VlnPlot(gex.Tcells, features = "CD4", group.by = "seurat_clusters", assay= "RNA", pt.size = 0)
VlnPlot(gex.Tcells, features = "CD8A", group.by = "seurat_clusters", assay= "RNA", pt.size = 0)
VlnPlot(gex.Tcells, features = "CD8B", group.by = "seurat_clusters", assay= "RNA", pt.size = 0)
dev.off()



saveRDS(gex.Tcells, file = "saveFiles/BaseTIL_Integrated_Tcells.rds", compress = "gzip")
