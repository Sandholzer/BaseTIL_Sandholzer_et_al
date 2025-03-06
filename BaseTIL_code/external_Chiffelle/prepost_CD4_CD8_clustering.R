
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(Seurat)

})

set.seed(12345678)

setwd("~/BaseTIL/Chiffelle")


prepost.sub <-  readRDS("saveFiles/prePost_Tcells_clustered.rds")

prepost.CD4 <- subset(prepost.sub, seurat_clusters %in% c(2,4,8))

prepost.CD8 <- subset(prepost.sub, seurat_clusters %in% c(1,3,6,7,9,10))




prepost.CD4$isCD8 <-
  ifelse(
    (prepost.CD4@assays[["RNA"]]$data["CD4",] < 0.1 &
       prepost.CD4@assays[["RNA"]]$data["CD8A",] > 0.1) |
      (prepost.CD4@assays[["RNA"]]$data["CD4",] < 0.1 &
         prepost.CD4@assays[["RNA"]]$data["CD8B",] > 0.1)
    ,TRUE,FALSE
  )

prepost.CD8$isCD4 <-
  ifelse(
    (prepost.CD8@assays[["RNA"]]$data["CD4",] > 0.1 &
       prepost.CD8@assays[["RNA"]]$data["CD8A",] < 0.1) |
      (prepost.CD8@assays[["RNA"]]$data["CD4",] > 0.1 &
         prepost.CD8@assays[["RNA"]]$data["CD8B",] < 0.1)
    ,TRUE, FALSE
  )


gex.CD8fromCD4 <- subset(prepost.CD4, subset = isCD8 == TRUE)
gex.CD4fromCD8 <- subset(prepost.CD8, subset = isCD4 == TRUE)


gex.trueCD4 <- subset(prepost.CD4, subset = isCD8 == FALSE)
gex.trueCD8 <- subset(prepost.CD8, subset = isCD4 == FALSE)


# intersect(colnames(gex.trueCD4), colnames(gex.CD4fromCD8))
# 
# 
# intersect(colnames(gex.trueCD4), colnames(gex.CD4fromCD8))
# 
# intersect(colnames(gex.trueCD8), colnames(gex.CD8fromCD4))
# 

prepost.CD4 <- merge(x= gex.trueCD4, y= gex.CD4fromCD8)
prepost.CD4[["RNA"]] <- JoinLayers(prepost.CD4[["RNA"]])
rm(gex.trueCD4, gex.CD4fromCD8)
gc()

prepost.CD8 <- merge(x= gex.trueCD8, y= gex.CD8fromCD4)
prepost.CD8[["RNA"]] <- JoinLayers(prepost.CD8[["RNA"]])
rm(gex.trueCD8, gex.CD8fromCD4)
gc()






# CD4 clustering 
#Integration 

prepost.CD4[["RNA"]] <- split(prepost.CD4[["RNA"]], f = prepost.CD4$orig.ident)

prepost.CD4 <- prepost.CD4 %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) 


gex.features <- VariableFeatures(prepost.CD4)
load("~/BaseTIL/Ranalysis/data/Zheng_code/exclude.gene.misc.human.v4.RData")


f.feat <- !(gex.features %in% all.gene.ignore.df[["seu.id"]]) &
  !(grepl("^(RP[LS]|Rp[ls])",gex.features,perl=T)) &
  !(grepl("^MT-",gex.features,perl=T)) &
  !(gex.features %in% c("MALAT1", "Malat1", "NEAT1"))
gex.features.flt <- gex.features[f.feat]
print(head(gex.features.flt))
print(paste("Length of anchors:", length(gex.features)) )
print(paste("Length of filtered anchors:", length(gex.features.flt)) )

rm(gex.features)

VariableFeatures(prepost.CD4) <- gex.features.flt 



prepost.CD4  <- ScaleData(prepost.CD4, 
                          vars.to.regress = c("S.Score","G2M.Score","percent_mito"), #"S.Score","G2M.Score","percent_mito" 
                          assay= "RNA") #, features = rownames(type.sub)

print("Normalization done.")

prepost.CD4 <- RunPCA(prepost.CD4, features = gex.features.flt, npcs = 50)

# DimHeatmap(prepost.CD4, dims = 1:15, cells = 500, balanced = TRUE)
# DimHeatmap(prepost.CD4, dims = 16:20, cells = 500, balanced = TRUE)


#batchV?
print("start Harmony") 


prepost.CD4 <- IntegrateLayers(
  object = prepost.CD4, 
  method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  assay = "RNA",
  #normalization.method = "SCT",
  features = gex.features.flt,
  verbose = T
)

prepost.CD4[["RNA"]] <- JoinLayers(prepost.CD4[["RNA"]])

ElbowPlot(prepost.CD4, ndims = 50) + geom_hline(yintercept=2, linetype="dashed", color = "red")


prepost.CD4 <- RunUMAP(prepost.CD4, reduction = "harmony", assay = "RNA", dims = 1:50)
prepost.CD4 <- FindNeighbors(prepost.CD4, dims = 1:50, reduction = "harmony")


prepost.CD4 <- Seurat::FindClusters(prepost.CD4, resolution = 0.5, algorithm = 4, method = "igraph") #



#Define cluster resolution here!
prepost.CD4$seurat_clusters <- prepost.CD4$RNA_snn_res.0.5
Idents(prepost.CD4) <- prepost.CD4$seurat_clusters



saveRDS(prepost.CD4, "saveFiles/prePost_CD4_clustered.rds")

pdf(paste0("Graphs/prePost_CD4_subclustering.pdf"))
# DimPlot(react.tum, reduction = "umap", group.by = "orig.ident", raster = F)
print(DimPlot(prepost.CD4, reduction = "umap", group.by = "seurat_clusters", label = T, raster = F))
print(DimPlot(prepost.CD4, reduction = "umap", group.by = "Patient", label = F, raster = F))
print(DimPlot(prepost.CD4, reduction = "umap", group.by = "Response", label = F, raster = F))
dev.off()



allMarkers <- FindAllMarkers(prepost.CD4, only.pos = T) #, test.use = "MAST"

#write.csv(allMarkers, file = paste0("Graphs/CD8_",Type_iterater,"_FindAllMarkers.csv"), row.names = F, quote = F)


load("~/BaseTIL/Ranalysis/data/Zheng_code/exclude.gene.misc.human.v4.RData")

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

png(paste0("Graphs/prePost_CD4_heatmap.png"), type = "cairo", height = 1200, width = 700)
print(DoHeatmap(subset(prepost.CD4, downsample = 200), features = top10$gene) + NoLegend())
dev.off()


png(paste0("Graphs/prePost_CD4_barPlot.png"), type = "cairo", height = 400, width = 800)
print(dittoSeq::dittoBarPlot(prepost.CD4, var = "seurat_clusters", group.by = "Patient", split.by = "Response"))
dev.off()






# CD8 clustering 
#Integration 

prepost.CD8[["RNA"]] <- split(prepost.CD8[["RNA"]], f = prepost.CD8$orig.ident)

prepost.CD8 <- prepost.CD8 %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) 


gex.features <- VariableFeatures(prepost.CD8)
load("~/BaseTIL/Ranalysis/data/Zheng_code/exclude.gene.misc.human.v4.RData")


f.feat <- !(gex.features %in% all.gene.ignore.df[["seu.id"]]) &
  !(grepl("^(RP[LS]|Rp[ls])",gex.features,perl=T)) &
  !(grepl("^MT-",gex.features,perl=T)) &
  !(gex.features %in% c("MALAT1", "Malat1", "NEAT1"))
gex.features.flt <- gex.features[f.feat]
print(head(gex.features.flt))
print(paste("Length of anchors:", length(gex.features)) )
print(paste("Length of filtered anchors:", length(gex.features.flt)) )

rm(gex.features)

VariableFeatures(prepost.CD8) <- gex.features.flt 



prepost.CD8  <- ScaleData(prepost.CD8, 
                          vars.to.regress = c("S.Score","G2M.Score","percent_mito"), #"S.Score","G2M.Score","percent_mito" 
                          assay= "RNA") #, features = rownames(type.sub)

print("Normalization done.")

prepost.CD8 <- RunPCA(prepost.CD8, features = gex.features.flt, npcs = 50)

# DimHeatmap(prepost.CD8, dims = 1:15, cells = 500, balanced = TRUE)
# DimHeatmap(prepost.CD8, dims = 16:20, cells = 500, balanced = TRUE)


#batchV?
print("start Harmony") 


prepost.CD8 <- IntegrateLayers(
  object = prepost.CD8, 
  method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  assay = "RNA",
  #normalization.method = "SCT",
  features = gex.features.flt,
  verbose = T
)

prepost.CD8[["RNA"]] <- JoinLayers(prepost.CD8[["RNA"]])

ElbowPlot(prepost.CD8, ndims = 50) + geom_hline(yintercept=2, linetype="dashed", color = "red")


prepost.CD8 <- RunUMAP(prepost.CD8, reduction = "harmony", assay = "RNA", dims = 1:50)
prepost.CD8 <- FindNeighbors(prepost.CD8, dims = 1:50, reduction = "harmony")


prepost.CD8 <- Seurat::FindClusters(prepost.CD8, resolution = 0.5, algorithm = 4, method = "igraph") #



#Define cluster resolution here!
prepost.CD8$seurat_clusters <- prepost.CD8$RNA_snn_res.0.5
Idents(prepost.CD8) <- prepost.CD8$seurat_clusters



saveRDS(prepost.CD8, "saveFiles/prePost_CD8_clustered.rds")

pdf(paste0("Graphs/prePost_CD8_subclustering.pdf"))
# DimPlot(react.tum, reduction = "umap", group.by = "orig.ident", raster = F)
print(DimPlot(prepost.CD8, reduction = "umap", group.by = "seurat_clusters", label = T, raster = F))
print(DimPlot(prepost.CD8, reduction = "umap", group.by = "Patient", label = F, raster = F))
print(DimPlot(prepost.CD8, reduction = "umap", group.by = "Response", label = F, raster = F))
dev.off()



allMarkers <- FindAllMarkers(prepost.CD8, only.pos = T) #, test.use = "MAST"

#write.csv(allMarkers, file = paste0("Graphs/CD8_",Type_iterater,"_FindAllMarkers.csv"), row.names = F, quote = F)


load("~/BaseTIL/Ranalysis/data/Zheng_code/exclude.gene.misc.human.v4.RData")

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

png(paste0("Graphs/prePost_CD8_heatmap.png"), type = "cairo", height = 1200, width = 700)
print(DoHeatmap(subset(prepost.CD8, downsample = 200), features = top10$gene) + NoLegend())
dev.off()


png(paste0("Graphs/prePost_CD8_barPlot.png"), type = "cairo", height = 400, width = 800)
print(dittoSeq::dittoBarPlot(prepost.CD8, var = "seurat_clusters", group.by = "Patient", split.by = "Response"))
dev.off()

