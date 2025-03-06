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

setwd("~/BaseTIL_code/CD8_Plotting/")
Exp.CD8 <-  readRDS(file = "../saveFiles/CD8_annotated.rds")
Exp.CD8 <- subset(Exp.CD8, Type=="Expanded_TILs" & Patient != "UPN002") #excluded to high gdT impair clustering

## -----------------------------------------------------------------------------------------------------------
#sub-clustering of TIL drug product
Exp.CD8[["RNA"]] <- split(Exp.CD8[["RNA"]], f = Exp.CD8$orig.ident)


Exp.CD8 <- Exp.CD8 %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) 

#filtering out unwanted genes
gex.features <- VariableFeatures(Exp.CD8)
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

VariableFeatures(Exp.CD8) <- gex.features.flt 

#regress out cell cycle, interferon response and percent mitochondrial genes
Exp.CD8  <- ScaleData(Exp.CD8, 
                       vars.to.regress = c("S.Score","G2M.Score" ,"ifn_genes1","percent_mito" ), 
                       assay= "RNA") 


Exp.CD8 <- RunPCA(Exp.CD8, features = gex.features.flt, npcs = 30)


#integration
Exp.CD8 <- IntegrateLayers(
  object = Exp.CD8, 
  method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  assay = "RNA",
  verbose = T
)

Exp.CD8[["RNA"]] <- JoinLayers(Exp.CD8[["RNA"]])

#clustering
Exp.CD8 <- RunUMAP(Exp.CD8, reduction = "harmony", assay = "RNA", dims = 2:14)
Exp.CD8 <- FindNeighbors(Exp.CD8, dims = 2:14, reduction = "harmony")
resolutions <- seq(0.2, 0.6, by=0.2)
Exp.CD8 <- Seurat::FindClusters(Exp.CD8, resolution = resolutions, algorithm = 4, method = "igraph") #


Exp.CD8$RNR <- ifelse(Exp.CD8$Response == "PR", "Rs", "NRs")

Exp.CD8$seurat_clusters<- Exp.CD8$RNA_snn_res.0.4
p1<- DimPlot(Exp.CD8, reduction = "umap", group.by = "Celltype", label = T, repel = T, raster = F, label.size = 9)+ ggeasy::easy_remove_axes() + ggtitle(NULL)
p2<- DimPlot(Exp.CD8, reduction = "umap", group.by = "seurat_clusters", label = T, repel = T, raster = F, label.size = 9)+ ggeasy::easy_remove_axes() + ggtitle(NULL)

p1+p2



#Annotate clusters
Exp.CD8 <- Seurat::SetIdent(Exp.CD8, value = Exp.CD8$seurat_clusters)

Exp.CD8 <-
  RenameIdents(
    Exp.CD8,
    `1` = "Tm",
    `2` = "Teff",
    `3` = "Teff",
    `4` = "Teff",
    `5` = "Tc17",
    `6` = "Tgd")


Exp.CD8$Celltype_TIL <- Idents(Exp.CD8)


Exp.CD8$Celltype_TIL <- factor(Exp.CD8$Celltype_TIL, levels = c("Teff",    "Tm",  "Tgd", "Tc17"  )) 

#save the file
saveRDS(Exp.CD8, "../saveFiles/CD8_Expanded_TILs_subclustering.rds")



Colors <- pal_jco(alpha = 0.7)(4)


#plot UMAP
pdf(file = "Fig/CD8_TIL_Annotated_UMAP.pdf", height = 4, width = 5)
DimPlot(Exp.CD8, reduction = "umap", group.by = "Celltype_TIL", label = T, repel = T, raster = F, cols = Colors, label.size = 8,label.box = T, alpha = 0.7)+NoLegend()+ 
  ggeasy::easy_remove_axes() + 
  ggtitle(NULL)

DimPlot(Exp.CD8, reduction = "umap", group.by = "Celltype_TIL", label = T, cols = Colors, repel = F, raster = F, label.size = 9, alpha = 0.7)+
  NoLegend()+ 
  ggeasy::easy_remove_axes() + 
  ggtitle(NULL) 
dev.off()


#identify cluster specific markers
allMarkers <- FindAllMarkers(Exp.CD8,  only.pos = T,  assay = "RNA") 
allMarkers$cluster <- factor(allMarkers$cluster, levels = c("Teff",    "Tm",  "Tgd", "Tc17" ))

write_csv(allMarkers, "Fig/DGE_results/CD8_TIL_prod_cluster.csv")


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

tmp.sup <- ScaleData(Exp.CD8, features = top10$gene)
tmp.sup <- subset(tmp.sup, downsample=200)


pdf("Fig/CD8_TIL_Heatmap_Markers_An.pdf", width = 9, height = 7)
DoHeatmap(tmp.sup, features = top10$gene, group.by = "Celltype_TIL",
          group.colors = Colors, angle = 0, hjust = 0.5, vjust = 0.5, disp.max = 1.7, disp.min = -1.6) #+ NoLegend()
dev.off()







## -----------------------------------------------------------------------------------------------------------
# Calculate Tc17 frequency

Exp.CD8@meta.data$isTc17 <- FALSE
Exp.CD8@meta.data[Exp.CD8@meta.data$Celltype_TIL == "Tc17","isTc17"] <- TRUE

dittoSeq::dittoBarPlot(Exp.CD8, var = "isTc17", group.by = "Patient", split.by = "Response")

data_Tc17<- dittoSeq::dittoBarPlot(Exp.CD8, var = "isTc17", group.by = "Patient", split.by = "RNR", data.out = T)

result_Tc17<- data_Tc17$data %>% filter(label==TRUE) %>% mutate(percent = percent*100)

# Calculate mean and standard error for each group
summary_data <- result_Tc17 %>%
  group_by(RNR) %>%
  summarise(
    mean_value = mean(percent),
    sd_value = sd(percent),
    se_value = sd_value / sqrt(length(percent))
  )

bar_colors <- c("NRs" = "#CD2626", "Rs" = "#009ACD")


# Create the plot
pdf("Fig/CD8_ExpTIL_Tc17_percent_response.pdf", height = 3, width = 2)
p <- ggplot(result_Tc17, aes(x = RNR, y = percent, fill = RNR)) +
  geom_bar(stat = "summary", fun = "mean", alpha = 0.4, width = 0.7, position = position_dodge(width = 0.8), color = NA, show.legend = FALSE) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, aes(color = RNR), position = position_dodge(width = 0.8), show.legend = FALSE) +
  geom_jitter(aes(color = RNR), position = position_jitterdodge(jitter.width = 0.35, dodge.width = 0.8, seed = 2), size = 1.5, show.legend = FALSE) +
  scale_fill_manual(name = "Response", values = bar_colors) +
  scale_color_manual(values = bar_colors) +
  theme_minimal() +
  labs(x = "", y = "Percent of CD8 T cells") +
  ggtitle("Tc17")+
  theme_classic()+
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = 16),
    axis.line = element_line(color = "black"),
    plot.title = element_text(size = 16,hjust = 0.5, face = "bold")
  ) 

p

dev.off()

