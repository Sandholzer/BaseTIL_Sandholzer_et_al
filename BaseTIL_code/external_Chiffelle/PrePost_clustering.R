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

list_files <- list.files("~/BaseTIL_code/external_Chiffelle/GSE222448_RAW_Barras/", pattern = "*CD45-counts.tsv.gz")
names(list_files) <- gsub("GSM[0-9]{7}_|\\-counts\\.tsv\\.gz$", "" ,list_files)


dir_list <- list() 

for (pat_num in 1:length(list_files)) {
  dir_list[[pat_num]] <- as.matrix(read.csv(paste0("~/BaseTIL_code/external_Chiffelle/GSE222448_RAW_Barras/" ,list_files[[pat_num]]), sep = "\t", header = T))
  names(dir_list)[[pat_num]] <- paste0("patient",pat_num)
}

prepost <- CreateSeuratObject(counts = dir_list, project = "Cheffelle", min.cells = 3, min.features = 200, names.field = 2, names.delim = "\\.")
table(prepost@meta.data$orig.ident)

#annotate samples
prepost@meta.data$Sample_Name <- prepost@meta.data$orig.ident
prepost@meta.data[c("Patient", "Type")] <- str_split_fixed(prepost@meta.data$orig.ident, '_', 2)
prepost@meta.data$Response <- NA
prepost@meta.data[prepost@meta.data$Patient %in% c("patient2", "patient3", "patient7", "patient8", "patient9", "patient13"), "Response"] <- "R"
prepost@meta.data[prepost@meta.data$Patient %in% c("patient1", "patient4", "patient10"), "Response"] <- "SD"
prepost@meta.data[prepost@meta.data$Patient %in% c("patient5", "patient6", "patient11", "patient12"), "Response"] <- "PD"

prepost@meta.data$RNR <- NA
prepost@meta.data[prepost@meta.data$Patient %in% c("patient2", "patient3", "patient7", "patient8", "patient9", "patient13"), "RNR"] <- "R"
prepost@meta.data[prepost@meta.data$Patient %in% c("patient1", "patient4", "patient10"), "RNR"] <- "NR"
prepost@meta.data[prepost@meta.data$Patient %in% c("patient5", "patient6", "patient11", "patient12"), "RNR"] <- "NR"

prepost@meta.data$prepost<- ifelse(prepost@meta.data$Type == "T30", "post", "pre")

head(prepost@meta.data)


## -----------------------------------------------------------------------------------------------------------
#Clustering of pre-ACT and post-ACT

prepost <- prepost %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData()

prepost[["RNA"]] <- JoinLayers(prepost[["RNA"]])

prepost <-
  PercentageFeatureSet(prepost, "MT-", col.name = "percent_mito")
prepost <-
  PercentageFeatureSet(prepost, "^RP[SL]", col.name = "percent_ribo")

prepost <-
  PercentageFeatureSet(prepost, "^HB[^(P)]", col.name = "percent_hb")

prepost <-
  PercentageFeatureSet(prepost, "PECAM1|PF4", col.name = "percent_plat")


Seurat::cc.genes.updated.2019
s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes

prepost <- Seurat::CellCycleScoring(prepost,
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
png(file="Graphs/QC_prepost_features_bevore_removal.png", width=1000, height=800 , type = "cairo")
VlnPlot(
  prepost,
  group.by = "Sample_Name",
  features = feats,
  pt.size = 0.1,
  ncol = 2
)+
  NoLegend()
#dev.off()


png(file="Graphs/QC_prepost_featuresVSumi_bevore_removal.png", width=500, height=1000, type = "cairo")
prepost@meta.data %>% 
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

prepost$log10GenesPerUMI <- log10(prepost$nFeature_RNA) / log10(prepost$nCount_RNA)

png(file="Graphs/QC_prepost_log10GenesPerUMI_bevore_removal.png", width=600, height=500, type = "cairo")
prepost@meta.data %>%
  ggplot(aes(x=log10GenesPerUMI, color = Sample_Name, fill=Sample_Name)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)+
  facet_wrap(~Sample_Name)
dev.off()


## -----------------------------------------------------------------------------------------------------------
#filtering was already done

table(prepost$Sample_Name)
#prepost <- subset(prepost, subset = nFeature_RNA > 650 & nCount_RNA > 1250 & percent_mito < 15  & percent_hb < 10 & percent_ribo > 5 )
table(prepost$Sample_Name)



## -----------------------------------------------------------------------------------------------------------  
#Clustering 

prepost[["RNA"]] <- split(prepost[["RNA"]], f = prepost$orig.ident)
prepost <- NormalizeData(prepost)
prepost <- FindVariableFeatures(prepost, selection.method = "vst", nfeatures = 2000)

gex.features <- VariableFeatures(prepost)
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

VariableFeatures(prepost) <- gex.features.flt

#regress for cell cycle and mito percent
prepost  <- ScaleData(prepost,
                     vars.to.regress = c("S.Score","G2M.Score" ,"ifn_genes1","percent_mito" ),
                     assay= "RNA") 

prepost <- RunPCA(prepost, features = gex.features.flt , npcs = 50)

ElbowPlot(prepost,ndims = 50)

prepost <- FindNeighbors(prepost, dims = 1:30)
prepost <- FindClusters(prepost, resolution = 0.5)
prepost <- RunUMAP(prepost, dims = 1:30)


DimPlot(prepost, group.by = "Patient")
DimPlot(prepost, split.by = "Response")

VlnPlot(prepost, "CD3E")




## -----------------------------------------------------------------------------------------------------------  
#Isolate Only T cell subset

prepost.sub <- subset(prepost, seurat_clusters %in% c(0,1,2,3,7,8,11,12,14,21))

prepost.sub <- prepost.sub %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) 

gex.features <- VariableFeatures(prepost.sub)
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

VariableFeatures(prepost.sub) <- gex.features.flt

#regress for cell cycle and mito percent
prepost.sub  <- ScaleData(prepost.sub,
                     vars.to.regress = c("S.Score","G2M.Score" ,"percent_mito" ),
                     assay= "RNA") 


prepost.sub <- RunPCA(prepost.sub, features = gex.features.flt, npcs = 20)


#Integration 
prepost.sub <- IntegrateLayers(
  object = prepost.sub, 
  method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  assay = "RNA",
  features = gex.features.flt,
  verbose = T
)
  
prepost.sub[["RNA"]] <- JoinLayers(prepost.sub[["RNA"]])

ElbowPlot(prepost.sub, ndims = 20) + geom_hline(yintercept=2, linetype="dashed", color = "red")


prepost.sub <- RunUMAP(prepost.sub, reduction = "harmony", assay = "RNA", dims = 1:15)
prepost.sub <- FindNeighbors(prepost.sub, dims = 1:15, reduction = "harmony")

resolutions <- seq(0.3, 0.6, by=0.1)

prepost.sub <- Seurat::FindClusters(prepost.sub, resolution = resolutions, algorithm = 4, method = "igraph") #



DimPlot(prepost.sub, reduction = "umap", raster = F, group.by = paste("RNA_snn_res.", resolutions, sep = ""))




#Define cluster resolution here!
prepost.sub$seurat_clusters <- prepost.sub$RNA_snn_res.0.5
Idents(prepost.sub) <- prepost.sub$seurat_clusters

#save file
saveRDS(prepost.sub, "saveFiles/prePost_Tcells_clustered.rds")

#plot UMAP

pdf(paste0("Graphs/prePost_Tcells_subclustering.pdf"))
print(ggarrange(plotlist = plt, ncol = 3, nrow = 2))
print(DimPlot(prepost.sub, reduction = "umap", group.by = "seurat_clusters", label = T, raster = F))
print(DimPlot(prepost.sub, reduction = "umap", group.by = "Patient", label = F, raster = F))
print(DimPlot(prepost.sub, reduction = "umap", group.by = "Response", label = F, raster = F))
dev.off()


#plot cluster specific markers
allMarkers <- FindAllMarkers(prepost.sub, only.pos = T)



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

png(paste0("Graphs/prePost_Tcells_heatmap.png"), type = "cairo", height = 1200, width = 700)
print(DoHeatmap(subset(prepost, downsample = 200), features = top10$gene) + NoLegend())
dev.off()


png(paste0("Graphs/prePost_Tcells_barPlot.png"), type = "cairo", height = 400, width = 800)
print(dittoSeq::dittoBarPlot(prepost, var = "seurat_clusters", group.by = "Patient", split.by = "Response"))
dev.off()

