## -----------------------------------------------------------------------------------------------------------
#load dependencies

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggpubr)
  library(readr)
  library(ggeasy)
  library(Seurat)
  library(UCell)
  library(harmony)

})

set.seed(12345678)


setwd("~/BaseTIL_code/CD8_Plotting/")

react.CD8 <-  readRDS(file = "~/BaseTIL_code/saveFiles/CD8_tumor_reactive.rds")


# remove files where low number tumor-reactive where identified
react.tum <- subset(react.CD8, Type %in% c("Tumor", "Rebiopsy1", "Rebiopsy2"))
react.tum <- subset(react.tum, Sample_Name %in% c( "UPN001 Rebiopsy2", "UPN008 Rebiopsy1"), invert=T)


#splitting files to perform normalization and variable feature defenition
react.tum[["RNA"]] <- split(react.tum[["RNA"]], f = react.tum$orig.ident)

react.tum <- react.tum %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) 


#remove unwanted genes from clustering
gex.features <- VariableFeatures(react.tum)
load("~/BaseTIL_code/data/Zheng_Japrin/exclude.gene.misc.human.v4.RData")


f.feat <- !(gex.features %in% all.gene.ignore.df[["seu.id"]]) &
  !(grepl("^(RP[LS]|Rp[ls])",gex.features,perl=T)) &
  !(grepl("^MT-",gex.features,perl=T)) &
  !(gex.features %in% c("MALAT1", "Malat1", "NEAT1"))
gex.features.flt <- gex.features[f.feat]
print(head(gex.features.flt))
print(paste("Length of anchors:", length(gex.features)) )
print(paste("Length of filtered anchors:", length(gex.features.flt)) )

rm(gex.features)

VariableFeatures(react.tum) <- gex.features.flt 


# scaling data and regress for unwanted effects
react.tum  <- ScaleData(react.tum, 
                      vars.to.regress = c("S.Score","G2M.Score" ,"ifn_genes1","percent_mito" ), 
                      features = rownames(react.tum),
                      assay= "RNA")



react.tum <- RunPCA(react.tum)


#integration using harmony 
react.tum <- IntegrateLayers(
  object = react.tum, 
  method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  assay = "RNA",
  verbose = T
)

react.tum[["RNA"]] <- JoinLayers(react.tum[["RNA"]])


#evaluate dimensions
ElbowPlot(react.tum, ndims = 30) + geom_hline(yintercept=2, linetype="dashed", color = "red")

DimHeatmap(object = react.tum, dims = 1:15, cells = 500, balanced = TRUE,fast = TRUE)
DimHeatmap(object = react.tum, dims = 16:20, cells = 500, balanced = TRUE,fast = TRUE)


react.tum <- RunUMAP(react.tum, reduction = "harmony", assay = "RNA", dims = 1:20)
react.tum <- FindNeighbors(react.tum, dims = 1:20, reduction = "harmony")

resolutions <- seq(0.3, 0.8, by=0.1)

react.tum <- Seurat::FindClusters(react.tum, resolution = resolutions, algorithm = 4) #, method = "igraph"


DimPlot(react.tum, reduction = "umap", raster = F, group.by = paste("RNA_snn_res.", resolutions, sep = ""))


#Define cluster resolution here!
react.tum$seurat_clusters <- react.tum$RNA_snn_res.0.6
Idents(react.tum) <- react.tum$seurat_clusters


pdf("Graphs/CD8_overallreactive_dimplot.pdf")
DimPlot(react.tum, reduction = "umap", group.by = "seurat_clusters", label = T, raster = F)
FeaturePlot(react.tum, c("PDCD1", "GNLY","CCL4", "GZMK",
                         "HAVCR2","IL7R", "CXCL13", "TCF7"), raster = F)
dev.off()


#save reclustering file
saveRDS(react.tum, file = "~/BaseTIL_code/saveFiles/CD8_tumor_reactive_clustering.rds", compress = "gzip")



#DGE to identify clusters
allMarkers <- FindAllMarkers(react.tum, only.pos = T) 
write.csv(allMarkers, file = "Graphs/CD8_TR_FindAllMarkers.csv", row.names = F, quote = F)


load("~/BaseTIL_code/data/Zheng_Japrin/exclude.gene.misc.human.v4.RData")
f.feat <- !(allMarkers$gene %in% all.gene.ignore.df[["seu.id"]]) &
  !(grepl("^(RP[LS]|Rp[ls])",allMarkers$gene,perl=T)) &
  !(grepl("^MT-",allMarkers$gene,perl=T)) &
  !(allMarkers$gene %in% c("MALAT1", "Malat1", "NEAT1"))
allMarkers <- allMarkers[f.feat,]


allMarkers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

png(paste0("Graphs/CD8_TR_heatmap.png"), type = "cairo", height = 1200, width = 700)
print(DoHeatmap(subset(react.tum, downsample = 200), features = top10$gene) + NoLegend())
dev.off()


