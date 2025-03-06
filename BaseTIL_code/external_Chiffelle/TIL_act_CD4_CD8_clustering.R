## -----------------------------------------------------------------------------------------------------------
#load dependencies

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(Seurat)

})

set.seed(12345678)

setwd("~/BaseTIL_code/external_Chiffelle")
acttil.sub <-  readRDS("saveFiles/Chiffelle_TIL_act.rds")


## -----------------------------------------------------------------------------------------------------------
#split CD4 and CD8 T cell clusters

acttil.CD4 <- subset(acttil.sub, seurat_clusters %in% c(2,5))

acttil.CD8 <- subset(acttil.sub, seurat_clusters %in% c(1,2,4,6,7))


acttil.CD4$isCD8 <-
  ifelse(
    (acttil.CD4@assays[["RNA"]]$data["CD4",] < 0.1 &
       acttil.CD4@assays[["RNA"]]$data["CD8A",] > 0.1) |
      (acttil.CD4@assays[["RNA"]]$data["CD4",] < 0.1 &
         acttil.CD4@assays[["RNA"]]$data["CD8B",] > 0.1)
    ,TRUE,FALSE
  )

acttil.CD8$isCD4 <-
  ifelse(
    (acttil.CD8@assays[["RNA"]]$data["CD4",] > 0.1 &
       acttil.CD8@assays[["RNA"]]$data["CD8A",] < 0.1) |
      (acttil.CD8@assays[["RNA"]]$data["CD4",] > 0.1 &
         acttil.CD8@assays[["RNA"]]$data["CD8B",] < 0.1)
    ,TRUE, FALSE
  )


gex.CD8fromCD4 <- subset(acttil.CD4, subset = isCD8 == TRUE)
gex.CD4fromCD8 <- subset(acttil.CD8, subset = isCD4 == TRUE)


gex.trueCD4 <- subset(acttil.CD4, subset = isCD8 == FALSE)
gex.trueCD8 <- subset(acttil.CD8, subset = isCD4 == FALSE)


acttil.CD4 <- merge(x= gex.trueCD4, y= gex.CD4fromCD8)
acttil.CD4[["RNA"]] <- JoinLayers(acttil.CD4[["RNA"]])
rm(gex.trueCD4, gex.CD4fromCD8)
gc()

acttil.CD8 <- merge(x= gex.trueCD8, y= gex.CD8fromCD4)
acttil.CD8[["RNA"]] <- JoinLayers(acttil.CD8[["RNA"]])
rm(gex.trueCD8, gex.CD8fromCD4)
gc()






## -----------------------------------------------------------------------------------------------------------
# CD4 clustering 

acttil.CD4[["RNA"]] <- split(acttil.CD4[["RNA"]], f = acttil.CD4$orig.ident)

acttil.CD4 <- acttil.CD4 %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) 


gex.features <- VariableFeatures(acttil.CD4)
load("~/BaseTIL_code/data/Zheng_Japrin/exclude.gene.misc.human.v4.RData")
f.feat <- !(gex.features %in% all.gene.ignore.df[["seu.id"]]) &
  !(grepl("^RP[LS]",gex.features,perl=T)) &
  !(grepl("^MT-",gex.features,perl=T)) &
  !(gex.features %in% c("MALAT1", "NEAT1"))
gex.features.flt <- gex.features[f.feat]
print(head(gex.features.flt))
print(paste("Length of anchors:", length(gex.features)) )
print(paste("Length of filtered anchors:", length(gex.features.flt)) )

rm(gex.features)

VariableFeatures(acttil.CD4) <- gex.features.flt 


#regress for cell cycle and percent mito
acttil.CD4  <- ScaleData(acttil.CD4, 
                          vars.to.regress = c("S.Score","G2M.Score","percent_mito"), 
                          assay= "RNA") 



acttil.CD4 <- RunPCA(acttil.CD4, features = gex.features.flt, npcs = 50)

#integration
acttil.CD4 <- IntegrateLayers(
  object = acttil.CD4, 
  method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  assay = "RNA",

  features = gex.features.flt,
  verbose = T
)

acttil.CD4[["RNA"]] <- JoinLayers(acttil.CD4[["RNA"]])

ElbowPlot(acttil.CD4, ndims = 50) + geom_hline(yintercept=2, linetype="dashed", color = "red")


acttil.CD4 <- RunUMAP(acttil.CD4, reduction = "harmony", assay = "RNA", dims = 1:50)
acttil.CD4 <- FindNeighbors(acttil.CD4, dims = 1:50, reduction = "harmony")


acttil.CD4 <- Seurat::FindClusters(acttil.CD4, resolution = 0.5, algorithm = 4, method = "igraph") #



#Define cluster resolution here!
acttil.CD4$seurat_clusters <- acttil.CD4$RNA_snn_res.0.5
Idents(acttil.CD4) <- acttil.CD4$seurat_clusters



saveRDS(acttil.CD4, "saveFiles/acttil_CD4_clustered.rds")

pdf(paste0("Graphs/acttil_CD4_subclustering.pdf"))
print(DimPlot(acttil.CD4, reduction = "umap", group.by = "seurat_clusters", label = T, raster = F))
print(DimPlot(acttil.CD4, reduction = "umap", group.by = "Patient", label = F, raster = F))
print(DimPlot(acttil.CD4, reduction = "umap", group.by = "Response", label = F, raster = F))
dev.off()



allMarkers <- FindAllMarkers(acttil.CD4, only.pos = T) 



load("~/BaseTIL_code/data/Zheng_Japrin/exclude.gene.misc.human.v4.RData")
f.feat <- !(allMarkers$gene %in% all.gene.ignore.df[["seu.id"]]) &
  !(grepl("^RP[LS]",allMarkers$gene,perl=T)) &
  !(grepl("^MT-",allMarkers$gene,perl=T)) &
  !(allMarkers$gene %in% c("MALAT1", "NEAT1"))
allMarkers <- allMarkers[f.feat,]


allMarkers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

png(paste0("Graphs/acttil_CD4_heatmap.png"), type = "cairo", height = 1200, width = 700)
print(DoHeatmap(subset(acttil.CD4, downsample = 200), features = top10$gene) + NoLegend())
dev.off()


png(paste0("Graphs/acttil_CD4_barPlot.png"), type = "cairo", height = 400, width = 800)
print(dittoSeq::dittoBarPlot(acttil.CD4, var = "seurat_clusters", group.by = "Patient", split.by = "Response"))
dev.off()






## -----------------------------------------------------------------------------------------------------------
# CD8 clustering 

acttil.CD8[["RNA"]] <- split(acttil.CD8[["RNA"]], f = acttil.CD8$orig.ident)

acttil.CD8 <- acttil.CD8 %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) 


gex.features <- VariableFeatures(acttil.CD8)
load("~/BaseTIL_code/data/Zheng_Japrin/exclude.gene.misc.human.v4.RData")
f.feat <- !(gex.features %in% all.gene.ignore.df[["seu.id"]]) &
  !(grepl("^RP[LS]",gex.features,perl=T)) &
  !(grepl("^MT-",gex.features,perl=T)) &
  !(gex.features %in% c("MALAT1", "NEAT1"))
gex.features.flt <- gex.features[f.feat]
print(head(gex.features.flt))
print(paste("Length of anchors:", length(gex.features)) )
print(paste("Length of filtered anchors:", length(gex.features.flt)) )

rm(gex.features)

VariableFeatures(acttil.CD8) <- gex.features.flt 


#regess for cell cycle and mito percent
acttil.CD8  <- ScaleData(acttil.CD8, 
                          vars.to.regress = c("S.Score","G2M.Score","percent_mito"), 
                          assay= "RNA") 


acttil.CD8 <- RunPCA(acttil.CD8, features = gex.features.flt, npcs = 50)

#integration
acttil.CD8 <- IntegrateLayers(
  object = acttil.CD8, 
  method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  assay = "RNA",
  features = gex.features.flt,
  verbose = T
)

acttil.CD8[["RNA"]] <- JoinLayers(acttil.CD8[["RNA"]])

ElbowPlot(acttil.CD8, ndims = 20) + geom_hline(yintercept=2, linetype="dashed", color = "red")


acttil.CD8 <- RunUMAP(acttil.CD8, reduction = "harmony", assay = "RNA", dims = 1:50)
acttil.CD8 <- FindNeighbors(acttil.CD8, dims = 1:50, reduction = "harmony")


acttil.CD8 <- Seurat::FindClusters(acttil.CD8, resolution = 0.5, algorithm = 4, method = "igraph") #



#Define cluster resolution here!
acttil.CD8$seurat_clusters <- acttil.CD8$RNA_snn_res.0.5
Idents(acttil.CD8) <- acttil.CD8$seurat_clusters


#save file
saveRDS(acttil.CD8, "saveFiles/acttil_CD8_clustered.rds")


#plot UMAP
pdf(paste0("Graphs/acttil_CD8_subclustering.pdf"))
print(DimPlot(acttil.CD8, reduction = "umap", group.by = "seurat_clusters", label = T, raster = F))
print(DimPlot(acttil.CD8, reduction = "umap", group.by = "Patient", label = F, raster = F))
print(DimPlot(acttil.CD8, reduction = "umap", group.by = "Response", label = F, raster = F))
dev.off()


#plot cluster specific markers
allMarkers <- FindAllMarkers(acttil.CD8, only.pos = T) 

load("~/BaseTIL_code/data/Zheng_Japrin/exclude.gene.misc.human.v4.RData")
f.feat <- !(allMarkers$gene %in% all.gene.ignore.df[["seu.id"]]) &
  !(grepl("^RP[LS]",allMarkers$gene,perl=T)) &
  !(grepl("^MT-",allMarkers$gene,perl=T)) &
  !(allMarkers$gene %in% c("MALAT1", "NEAT1"))
allMarkers <- allMarkers[f.feat,]


allMarkers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

png(paste0("Graphs/acttil_CD8_heatmap.png"), type = "cairo", height = 1200, width = 700)
print(DoHeatmap(subset(acttil.CD8, downsample = 200), features = top10$gene) + NoLegend())
dev.off()


png(paste0("Graphs/acttil_CD8_barPlot.png"), type = "cairo", height = 400, width = 800)
print(dittoSeq::dittoBarPlot(acttil.CD8, var = "seurat_clusters", group.by = "Patient", split.by = "Response"))
dev.off()

