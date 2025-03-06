## -----------------------------------------------------------------------------------------------------------
#load dependencies

suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
  library(harmony)
  library(stringr)
  
})

set.seed(12345678)

setwd("~/BaseTIL_code/external_Chiffelle")


## -----------------------------------------------------------------------------------------------------------
#load the raw data


list_files <- list.files("~/BaseTIL_code/external_Chiffelle/GSE229861_RAW_Chiffelle/", pattern = "*counts.tsv.gz")
names(list_files) <- gsub("GSM[0-9]{7}_|\\-counts\\.tsv\\.gz$", "" ,list_files)


dir_list <- list() 

for (pat_num in 1:13) {
  dir_list[[pat_num]] <- as.matrix(read.csv(paste0("~/BaseTIL_code/external_Chiffelle/GSE229861_RAW_Chiffelle/" ,list_files[[pat_num]]), sep = "\t", header = T))
  names(dir_list)[[pat_num]] <- paste0("patient",pat_num)
  
}


acttil <- CreateSeuratObject(counts = dir_list, project = "Cheffelle", min.cells = 3, min.features = 200, names.field = 2, names.delim = "\\.")
acttil


#annotate samples
acttil@meta.data$Sample_Name <- acttil@meta.data$orig.ident
acttil@meta.data[c("Patient", "Type")] <- str_split_fixed(acttil@meta.data$orig.ident, '_', 2)
acttil@meta.data$Response <- NA
acttil@meta.data[acttil@meta.data$Patient %in% c("patient2", "patient3", "patient7", "patient8", "patient9", "patient13"), "Response"] <- "R"
acttil@meta.data[acttil@meta.data$Patient %in% c("patient1", "patient4", "patient10"), "Response"] <- "SD"
acttil@meta.data[acttil@meta.data$Patient %in% c("patient5", "patient6", "patient11", "patient12"), "Response"] <- "PD"

acttil@meta.data$RNR <- NA
acttil@meta.data[acttil@meta.data$Patient %in% c("patient2", "patient3", "patient7", "patient8", "patient9", "patient13"), "RNR"] <- "R"
acttil@meta.data[acttil@meta.data$Patient %in% c("patient1", "patient4", "patient10"), "RNR"] <- "NR"
acttil@meta.data[acttil@meta.data$Patient %in% c("patient5", "patient6", "patient11", "patient12"), "RNR"] <- "NR"


head(acttil@meta.data)

## -----------------------------------------------------------------------------------------------------------
#Clustering of TIL product


acttil <- acttil %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData()

acttil[["RNA"]] <- JoinLayers(acttil[["RNA"]])


acttil <-
  PercentageFeatureSet(acttil, "MT-", col.name = "percent_mito")
acttil <-
  PercentageFeatureSet(acttil, "^RP[SL]", col.name = "percent_ribo")


acttil <-
  PercentageFeatureSet(acttil, "^HB[^(P)]", col.name = "percent_hb")

acttil <-
  PercentageFeatureSet(acttil, "PECAM1|PF4", col.name = "percent_plat")


Seurat::cc.genes.updated.2019
s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes

acttil <- Seurat::CellCycleScoring(acttil,
                                     s.features = s.genes,
                                     g2m.features = g2m.genes, search = T)





## -----------------------------------------------------------------------------------------------------------
#Plotting of Parameters for Filter selection

feats <-
  c("nFeature_RNA",
    "nCount_RNA",
    "percent_mito",
    "percent_ribo",
    "percent_hb", "percent_plat")
png(file="Graphs/QC_TIL_features_bevore_removal.png", width=1000, height=800 , type = "cairo")
VlnPlot(
  acttil,
  group.by = "Sample_Name",
  features = feats,
  pt.size = 0.1,
  ncol = 2
)+
  NoLegend()
dev.off()


png(file="Graphs/QC_TIL_featuresVSumi_bevore_removal.png", width=500, height=1000, type = "cairo")
acttil@meta.data %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent_mito)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 1250) +
  geom_hline(yintercept = 650) +
  facet_wrap(~Sample_Name)
dev.off()

acttil$log10GenesPerUMI <- log10(acttil$nFeature_RNA) / log10(acttil$nCount_RNA)

png(file="Graphs/QC_TIL_log10GenesPerUMI_bevore_removal.png", width=600, height=500, type = "cairo")
acttil@meta.data %>%
  ggplot(aes(x=log10GenesPerUMI, color = Sample_Name, fill=Sample_Name)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)+
  facet_wrap(~Sample_Name)
dev.off()


## -----------------------------------------------------------------------------------------------------------
#filtering was already done

table(acttil$Sample_Name)
#acttil <- subset(acttil, subset = nFeature_RNA > 650 & nCount_RNA > 1250 & percent_mito < 15  & percent_hb < 10 & percent_ribo > 5 )
table(acttil$Sample_Name)




## -----------------------------------------------------------------------------------------------------------  
#Clustering 


acttil[["RNA"]] <- split(acttil[["RNA"]], f = acttil$orig.ident)
  acttil <- acttil %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 3000) 
  
  
  gex.features <- VariableFeatures(acttil)
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
  
  VariableFeatures(acttil) <- gex.features.flt 
  
  
  #regress for cell cycle and mito percent
  acttil  <- ScaleData(acttil, 
                         vars.to.regress = c("S.Score","G2M.Score","percent_mito"),
                         assay= "RNA") 

  acttil <- RunPCA(acttil, features = gex.features.flt, npcs = 20)
  
  #Integration
  acttil <- IntegrateLayers(
    object = acttil, 
    method = HarmonyIntegration,
    orig.reduction = "pca", new.reduction = "harmony",
    assay = "RNA",
    features = gex.features.flt,
    verbose = T, kmeans_init_nstart=20, kmeans_init_iter_max=100
  )
  
  acttil[["RNA"]] <- JoinLayers(acttil[["RNA"]])
  
  ElbowPlot(acttil, ndims = 20) + geom_hline(yintercept=2, linetype="dashed", color = "red")
  

  acttil <- RunUMAP(acttil, reduction = "harmony", assay = "RNA", dims = 1:15)
  acttil <- FindNeighbors(acttil, dims = 1:15, reduction = "harmony")
  
  resolutions <- seq(0.3, 0.6, by=0.1)
  
  acttil <- Seurat::FindClusters(acttil, resolution = resolutions, algorithm = 4, method = "igraph") #
  
  
  DimPlot(acttil, reduction = "umap", raster = F, group.by = paste("RNA_snn_res.", resolutions, sep = ""))
  
  

  #Define cluster resolution here!
  acttil$seurat_clusters <- acttil$RNA_snn_res.0.5
  Idents(acttil) <- acttil$seurat_clusters
  
  #save file
  saveRDS(acttil, "saveFiles/Chiffelle_TIL_act.rds")
  

  
  #plot UMAP
  pdf(paste0("Graphs/TIL_subclustering.pdf"))
  print(ggarrange(plotlist = plt, ncol = 3, nrow = 2))
  print(DimPlot(acttil, reduction = "umap", group.by = "seurat_clusters", label = T, raster = F))
  print(DimPlot(acttil, reduction = "umap", group.by = "Patient", label = F, raster = F))
  print(DimPlot(acttil, reduction = "umap", group.by = "Response", label = F, raster = F))

  dev.off()
  
  
  
  #plot cluster specific markers
  allMarkers <- FindAllMarkers(acttil, only.pos = T) 

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
  
  png(paste0("Graphs/TIL_heatmap.png"), type = "cairo", height = 1200, width = 700)
  print(DoHeatmap(subset(acttil, downsample = 200), features = top10$gene) + NoLegend())
  dev.off()
  
  
  png(paste0("Graphs/TIL_barPlot.png"), type = "cairo", height = 400, width = 800)
  print(dittoSeq::dittoBarPlot(acttil, var = "seurat_clusters", group.by = "Patient", split.by = "Response"))
  dev.off()
  
