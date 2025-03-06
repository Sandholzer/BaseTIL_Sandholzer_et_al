## -----------------------------------------------------------------------------------------------------------
#load dependencies
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggeasy)
  library(Seurat) 
  library(harmony)
  library(readr)
  library(ggpubr)
  library(viridis)
  library(purrr)
  library(tidyr)
  library(dplyr)
  library(AUCell)
  library(RColorBrewer)
  library(plotly)
  library(grid)
  library(ggsci)
  
})

set.seed(12345678)
setwd("~/BaseTIL_code/CD4_Plotting/")

Exp.CD4 <-  readRDS(file = "~/BaseTIL_code/saveFiles/CD4_3_QC_Harm_ann.rds")
Exp.CD4 <- subset(Exp.CD4, Type=="Expanded_TILs" & Patient != "UPN002")

## -----------------------------------------------------------------------------------------------------------
#sub-clustering of TIL drug product
Exp.CD4[["RNA"]] <- JoinLayers(Exp.CD4[["RNA"]])
Exp.CD4[["RNA"]] <- split(Exp.CD4[["RNA"]], f = Exp.CD4$orig.ident)


Exp.CD4 <- Exp.CD4 %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) 

#filtering out unwanted genes
gex.features <- VariableFeatures(Exp.CD4)
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

VariableFeatures(Exp.CD4) <- gex.features.flt 


#regress out cell cycle, interferon response and percent mitochondrial genes
Exp.CD4  <- ScaleData(Exp.CD4, 
                      vars.to.regress = c("S.Score","G2M.Score" ,"ifn_genes1","percent_mito" ), 
                      assay= "RNA")

Exp.CD4 <- RunPCA(Exp.CD4, features = gex.features.flt, npcs = 30)


#integration
Exp.CD4 <- IntegrateLayers(
  object = Exp.CD4, 
  method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  assay = "RNA",
  verbose = T
)

Exp.CD4[["RNA"]] <- JoinLayers(Exp.CD4[["RNA"]])


#clustering
Exp.CD4 <- RunUMAP(Exp.CD4, reduction = "harmony", assay = "RNA", dims = 2:14)
Exp.CD4 <- FindNeighbors(Exp.CD4, dims = 2:14, reduction = "harmony")
resolutions <- seq(0.2, 1.0, by=0.2)
Exp.CD4 <- Seurat::FindClusters(Exp.CD4, resolution = resolutions, algorithm = 1) #4, method = "igraph"


Exp.CD4$RNR <- ifelse(Exp.CD4$Response == "PR", "Rs", "NRs")

Exp.CD4$seurat_clusters<- Exp.CD4$RNA_snn_res.0.4
p1<- DimPlot(Exp.CD4, reduction = "umap", group.by = "Celltype", label = T, repel = T, raster = F, label.size = 9)+ ggeasy::easy_remove_axes() + ggtitle(NULL)
p2<- DimPlot(Exp.CD4, reduction = "umap", group.by = "seurat_clusters", label = T, repel = T, raster = F, label.size = 9)+ ggeasy::easy_remove_axes() + ggtitle(NULL)

p1+p2


#Annotate clusters
Exp.CD4 <- Seurat::SetIdent(Exp.CD4, value = Exp.CD4$seurat_clusters)

Exp.CD4 <-
  RenameIdents(
    Exp.CD4,
    `0` = "Teff",
    `1` = "Teff",
    `2` = "Th17",
    `3` = "Tm",
    `4` = "Th2")


Exp.CD4$Celltype_TIL <- Idents(Exp.CD4)

Exp.CD4$Celltype_TIL <- factor(Exp.CD4$Celltype_TIL, levels = c("Teff",    "Tm",  "Th2", "Th17"   )) 


#save the file
saveRDS(Exp.CD4, "../saveFiles/CD4_Expanded_TILs_subclustering.rds")


Colors <- pal_jco(alpha = 0.7)(5)


#plot UMAP
pdf(file = "Fig/CD4_TIL_Annotated_UMAP.pdf", height = 4, width = 4)
DimPlot(Exp.CD4, reduction = "umap", group.by = "Celltype_TIL", label = T, repel = F, raster = F, cols = Colors, label.size = 8,label.box = T, alpha = 0.7)+NoLegend()+ ggeasy::easy_remove_axes() + ggtitle(NULL)

DimPlot(Exp.CD4, reduction = "umap", group.by = "Celltype_TIL", label = T, cols = Colors,repel = T, raster = F, label.size = 9, alpha = 0.7)+
  NoLegend()+ 
  ggeasy::easy_remove_axes() + 
  ggtitle(NULL) 
dev.off()





#identify cluster specific markers
allMarkers <- FindAllMarkers(Exp.CD4,  only.pos = T,  assay = "RNA") #test.use = "MAST",
gex.features.flt$cluster <- factor(gex.features.flt$cluster, levels = c("Teff",    "Tm",  "Th2", "Th17"   ))

write_csv(allMarkers, "Fig/DGE_results/CD4_TIL_prod_cluster.csv")

gex.features <- allMarkers$gene
load("~/BaseTIL_code/data/Zheng_Japrin/exclude.gene.misc.human.v4.RData")
f.feat <- !(gex.features %in% all.gene.ignore.df[["seu.id"]]) &
  !(grepl("^RP[LS]",gex.features,perl=T)) &
  !(grepl("^MT-",gex.features,perl=T)) &
  !(gex.features %in% c("MALAT1", "NEAT1"))
gex.features.flt <- allMarkers[f.feat,]
print(head(gex.features.flt))
print(paste("Length of anchors:", length(gex.features)) )
print(paste("Length of filtered anchors:", length(gex.features.flt)) )



gex.features.flt %>%
  group_by(cluster) %>%
  #dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

tmp.sup <- ScaleData(Exp.CD4, features = top10$gene)
tmp.sup <- subset(tmp.sup, downsample=200)


pdf("Fig/CD4_TIL_Heatmap_Markers_An.pdf", width = 9, height = 6)
DoHeatmap(tmp.sup, features = top10$gene, group.by = "Celltype_TIL", group.colors = Colors, angle = 0, hjust = 0.5, vjust = 0.5, raster = T) #+ NoLegend()
dev.off()





## -----------------------------------------------------------------------------------------------------------
# Calculate Th17 frequency

Exp.CD4@meta.data$isTh17 <- FALSE
Exp.CD4@meta.data[Exp.CD4@meta.data$Celltype_TIL == "Th17","isTh17"] <- TRUE

dittoSeq::dittoBarPlot(Exp.CD4, var = "isTh17", group.by = "Patient", split.by = "Response")

data_Th17<- dittoSeq::dittoBarPlot(Exp.CD4, var = "isTh17", group.by = "Patient", split.by = "RNR", data.out = T)

result_Th17<- data_Th17$data %>% filter(label==TRUE) %>% mutate(percent = percent*100)

# Calculate mean and standard error for each group
summary_data <- result_Th17 %>%
  group_by(RNR) %>%
  summarise(
    mean_value = mean(percent),
    sd_value = sd(percent),
    se_value = sd_value / sqrt(length(percent))
  )

bar_colors <- c("NRs" = "#CD2626", "Rs" = "#009ACD")


# Create the plot
pdf("Fig/CD4_ExpTIL_Th17_percent_response.pdf", height = 3, width = 2)
p <- ggplot(result_Th17, aes(x = RNR, y = percent, fill = RNR)) +
  geom_bar(stat = "summary", fun = "mean", alpha = 0.4, width = 0.7, position = position_dodge(width = 0.8), color = NA, show.legend = FALSE) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, aes(color = RNR), position = position_dodge(width = 0.8), show.legend = FALSE) +
  geom_jitter(aes(color = RNR), position = position_jitterdodge(jitter.width = 0.35, dodge.width = 0.8, seed = 1), size = 1.5, show.legend = FALSE) +
  scale_fill_manual(name = "Response", values = bar_colors) +
  scale_color_manual(values = bar_colors) +
  theme_minimal() +
  labs(x = "", y = "Percent of CD4 T cells") +
  ggtitle("Th17")+
  theme_classic()+
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = 16),
    axis.line = element_line(color = "black"),
    plot.title = element_text(size = 16,hjust = 0.5, face = "bold")
  ) 

p
p+stat_compare_means(aes(group = RNR), label = "p.signif", method = "t.test")
dev.off()

