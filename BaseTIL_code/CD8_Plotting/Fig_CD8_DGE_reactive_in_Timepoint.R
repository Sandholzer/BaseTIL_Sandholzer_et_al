## -----------------------------------------------------------------------------------------------------------
#load dependencies


suppressPackageStartupMessages({
  library(ggplot2)
  library(Seurat)
  library(stringr)
  library(ggrepel)
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(ggnewscale)
})

set.seed(1234567)
setwd("~/BaseTIL_code/CD8_Plotting")

gex.CD8 <- readRDS(file = "../saveFiles/CD8_annotated.rds")

#remove control PBMC sample from analysis
gex.CD8 <- subset(gex.CD8, Type != "PBMC_Ctrl")

#prevents problems with NAs
gex.CD8$reactive_clones <- ifelse(gex.CD8$overall_reactive == TRUE & !is.na(gex.CD8$overall_reactive), TRUE, FALSE)


for (arg_type in c( "Tumor" , "PreREP" , "Expanded_TILs", "PBMC_7dpt","Post")) {
  


print(paste("DGE of:", arg_type))

if (arg_type %in% c("post", "Post")) {
  arg_type2 <- c("Rebiopsy1", "Rebiopsy2")
} else{
  arg_type2 <- arg_type
}



print(paste("Starting:", arg_type))

# Filtering for differential expression analysis
tum.CD8 <- subset(gex.CD8, Type %in% arg_type2 & Patient %in% c("UPN001", "UPN003", "UPN006", "UPN008", "UPN009", "UPN011")) #,"UPN002" removed due to low cell number and high gdT cells frequency
tum.CD8 <- subset(tum.CD8, Sample_Name %in% c("UPN008 Rebiopsy1", "UPN001 Rebiopsy2", "UPN006 Rebiopsy1", "UPN006 Rebiopsy2"), invert = T) # samples with low numbers of detected tumor-reactive clones were removed


Idents(tum.CD8) <- tum.CD8$reactive_clones


#DGE using single cells
markers <- FindMarkers(tum.CD8, ident.1 = "TRUE",
                       ident.2 = "FALSE", logfc.threshold = 0) 
markers$gene <- rownames(markers)
write.csv(markers, paste0("DGE_results/dge_", arg_type, ".csv"), row.names = F)


#Filtering of non relevant genes
load("../data/Zheng_Japrin/exclude.gene.misc.human.v4.RData")


f.feat <- !(markers$gene %in% all.gene.ignore.df[["seu.id"]]) &
  !(grepl("^(RP[LS]|Rp[ls])",markers$gene,perl=T)) &
  !(grepl("^MT-|^MTRN|^MTAT|^MTND|^MRP",markers$gene,perl=T)) &
  !(grepl("^RP[0-9]+-",markers$gene,perl=T)) &
  !(markers$gene %in% c("MALAT1", 'XIST', 'NEAT1'))
markers.flt <- markers[f.feat,]
print(head(markers.flt))
print(paste("Length of anchors:", nrow(markers)) )
print(paste("Length of filtered anchors:", nrow(markers.flt)) )



#aggregate based on patient for DGE as patient replicate
tum.CD8_aggregate <- AggregateExpression(tum.CD8, group.by = c("Sample_Name", "reactive_clones"), return.seurat = T)

counts <- tum.CD8@meta.data %>% group_by(Sample_Name, reactive_clones) %>% summarise(cell_count = n())
counts$reactive_clones <- factor(counts$reactive_clones)
tum.CD8_aggregate@meta.data$reactive_clones <- factor(tum.CD8_aggregate@meta.data$reactive_clones)
colData_meta <- left_join(tum.CD8_aggregate@meta.data, counts, by = c("Sample_Name", "reactive_clones"))



dds <- DESeqDataSetFromMatrix(tum.CD8_aggregate[["RNA"]]$counts, 
                              colData = colData_meta, 
                              design = ~ reactive_clones + Sample_Name)
rld <- rlog(dds, blind=TRUE)

dds <- DESeq(dds)
plotDispEsts(dds)


resultsNames(dds)

res <- results(dds, 
               name = "reactive_clones_TRUE_vs_FALSE",
               alpha = 0.1)


library(tibble)
res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  arrange(padj)


res_table_thres <- res_tbl[!is.na(res_tbl$pvalue), ] 

write.csv(res_tbl, paste0("DGE_results/dge_", arg_type, "_bulk_new.csv"), row.names = F)


#Filtering of non relevant genes
load("../data/Zheng_Japrin/exclude.gene.misc.human.v4.RData")

f.feat <- !(res_table_thres$gene %in% all.gene.ignore.df[["seu.id"]]) &
  !(grepl("^(RP[LS]|Rp[ls])",res_table_thres$gene,perl=T)) &
  !(grepl("^MT-|^MTRN|^MTAT|^MTND|^MRP",res_table_thres$gene,perl=T)) &
  !(grepl("^RP[0-9]+-",res_table_thres$gene,perl=T)) &
  !(res_table_thres$gene %in% c("MALAT1", 'XIST', 'NEAT1'))
markers_bulk.flt <- res_table_thres[f.feat,]
print(head(markers_bulk.flt))
print(paste("Length of anchors:", nrow(res_table_thres)) )
print(paste("Length of filtered anchors:", nrow(markers_bulk.flt)) )


## -----------------------------------------------------------------------------------------------------------
#general plotting for data exploration

pdf(paste0("DGE_results/DGE_react_", arg_type, ".pdf"), width = 7, height = 9)
markers.flt %>%
  dplyr::rename(FC = avg_log2FC) %>%
  dplyr::mutate(name = ifelse(gene %in% c(markers.flt$gene[markers.flt$p_val_adj < 0.01 & abs(markers.flt$avg_log2FC) > 0.5]), gene, NA)) %>%
  dplyr::mutate(sig = ifelse(p_val_adj < 0.01 & abs(FC) > 0.15, "yes", "no")) %>%
  ggplot(aes(FC, -log10(p_val_adj))) +
  geom_point(aes(color = sig), size = 0.9) +
  geom_text_repel(aes(label = name),
                  parse = F,
                  size = 3.5,
                  box.padding = 0.52,
                  segment.color = "black",
                  show.legend = FALSE) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 15,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  scale_color_manual(values = c("#D6D6D6","#CD2626"))+ ggtitle("single-Cell")

EnhancedVolcano::EnhancedVolcano(markers.flt, x = "avg_log2FC", y = "p_val_adj", lab = markers.flt$gene)

markers_bulk.flt %>%
  dplyr::rename(FC = log2FoldChange) %>%
  dplyr::mutate(name = ifelse(gene %in% c(markers_bulk.flt$gene[markers_bulk.flt$padj < 0.01 & abs(markers_bulk.flt$log2FoldChange) > 0.5]), gene, NA)) %>%
  dplyr::mutate(sig = ifelse(padj < 0.01 & abs(FC) > 0.15, "yes", "no")) %>%
  ggplot(aes(FC, -log10(padj))) +
  geom_point(aes(color = sig), size = 0.9) +
  geom_text_repel(aes(label = name),
                  parse = F,
                  size = 3.5,
                  box.padding = 0.52,
                  segment.color = "black",
                  show.legend = FALSE) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 15,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  scale_color_manual(values = c("#D6D6D6","#CD2626")) + ggtitle("bulk")

EnhancedVolcano::EnhancedVolcano(markers_bulk.flt, x = "log2FoldChange", y = "padj", lab = markers_bulk.flt$gene, pCutoff = 0.3, FCcutoff = 0.25)


dev.off()

}


## -----------------------------------------------------------------------------------------------------------
#Plotting the graphs used in manuscript for baseline Tumor


# loading the single cell data for baseline Tumor
markers  <- read.csv("DGE_results/dge_Tumor.csv")


#now for pseudo bulk results
markers_bulk  <- read.csv("DGE_results/dge_Tumor_bulk_new.csv")


markers_bulk$p_val <- markers_bulk$pvalue
markers_bulk$avg_log2FC <- markers_bulk$log2FoldChange
markers_bulk$p_val_adj <- markers_bulk$padj

merged <- left_join(markers, markers_bulk, by = "gene", suffix = c(".sc", ".bulk") )


#filter the merged table
load("../data/Zheng_Japrin/exclude.gene.misc.human.v4.RData")
f.feat <- !(merged$gene %in% all.gene.ignore.df[["seu.id"]]) &
  !(grepl("^(RP[LS]|Rp[ls])",merged$gene,perl=T)) &
  !(grepl("^MT-",merged$gene ,perl=T)) &
  !(merged$gene %in% c("MALAT1", "NEAT1"))
print(paste("Length of anchors:", length(merged$gene)) )
merged <- merged[f.feat,]
print(head(merged))
print(paste("Length of filtered anchors:", length(merged$gene)) )



#for annotation based on threshold
merged$name <- ifelse(merged$p_val_adj.sc < 1e-50 & merged$p_val_adj.bulk < 0.05 ,merged$gene, NA )


#select markers to annotate
own.gene <- c("KIR2DL4", "ETV1") #"CD109",,  "ANXA9", "TNS3""TNFRSF9","KLRC2",,"BTAF""LAYN",,"LAG3""VCAM1",
exhaustion.gene <- c("SNX9","CXCR6","TOX","CD27",  "DUSP4","ITGAE","CTLA4","PDCD1","ENTPD1","TIGIT","HAVCR2", "CXCL13", "TNFRSF9")
mem.gene <- c("IL7R","SELL", "CCR7", "KLF2", "S1PR1", "CCR4", "TCF7", "KLF3")
selection <- c(own.gene, exhaustion.gene, mem.gene)

merged$name <- ifelse(merged$gene %in% selection ,merged$gene, NA )


#define if dots should be shown based on threshold
merged$show <-  ifelse(merged$p_val_adj.sc < 1e-25 , "yes", "no") 
merged$show <- factor(merged$show, levels = c("no", "yes"))

merged$sig <-  ifelse(merged$p_val_adj.sc < 1e-50 & merged$p_val_adj.bulk < 0.05, TRUE, FALSE) 



#define value for color scaling 
merged$color_scale <- (merged$avg_log2FC.sc + merged$avg_log2FC.bulk)
merged[merged$color_scale>1.5,"color_scale"] <- 1.5
merged[merged$color_scale<= -1.5,"color_scale"] <- -1.5

write.csv(merged, paste0("Fig/DGE_results/DGE_Reactive_correlation_pre-ACT_Tumor.csv"), row.names = F)


pdf("Fig/DGE_results/DGE_Reactive_correlation_Tumor.pdf", width = 5, height = 5)
plt1<- ggplot(merged) + 
  geom_point(data = merged[merged$sig == FALSE &merged$show == "yes",],aes(x=avg_log2FC.sc, y=avg_log2FC.bulk, color = color_scale), size = 1)+ 
  scale_colour_gradient2(
    low = "#7E8CBE",
    mid = "gray90",
    high = "#F8766D",
    midpoint = 0,
    na.value = "grey50",
    limits=c(-1.5,1.5)
  ) +
  new_scale_colour() +
  geom_point(data = merged[merged$sig == TRUE & merged$show == "yes",], aes(x=avg_log2FC.sc, y=avg_log2FC.bulk, color = color_scale),  size = 1) + 
  scale_colour_gradient2(
    low = "#3C4DA0",
    mid = "gray90",
    high = "#CD1F2D",
    midpoint = 0,
    na.value = "grey50",
    limits=c(-1.5,1.5))+
  theme_classic()+
  geom_hline(yintercept = 0, linetype = "dotted", color = "black")+
  geom_vline(xintercept = 0, linetype = "dotted", color = "black")+
  geom_label_repel(aes(label = name,x=avg_log2FC.sc, y=avg_log2FC.bulk),
                   parse = F,
                   size = 3,
                   box.padding = 0.25,
                   segment.color = "black",
                   min.segment.length = 0,max.time = 1,
                   show.legend = FALSE, max.overlaps = 20, point.padding = 0, label.padding = 0.15) + NoLegend() + 
  xlab("single-cell log2FC") + 
  ylab("patient averaged log2FC") + xlim(-5, 4) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 15,color="black"),
    axis.text = element_text(size = 12,color="black")
  )

plt1
dev.off()









## -----------------------------------------------------------------------------------------------------------
#Plotting the graphs used in manuscript for preREP (intermediate product)


# loading the single cell data for PreREP
markers  <- read.csv("DGE_results/dge_PreREP.csv")


#now for pseudo bulk results
markers_bulk  <- read.csv("DGE_results/dge_PreREP_bulk_new.csv")

markers_bulk$p_val <- markers_bulk$pvalue
markers_bulk$avg_log2FC <- markers_bulk$log2FoldChange
markers_bulk$p_val_adj <- markers_bulk$padj


merged <- left_join(markers, markers_bulk, by = "gene", suffix = c(".sc", ".bulk") )

#filter the merged table
load("../data/Zheng_Japrin/exclude.gene.misc.human.v4.RData")
f.feat <- !(merged$gene %in% all.gene.ignore.df[["seu.id"]]) &
  !(grepl("^(RP[LS]|Rp[ls])",merged$gene,perl=T)) &
  !(grepl("^MT-",merged$gene ,perl=T)) &
  !(merged$gene %in% c("MALAT1", "NEAT1"))
print(paste("Length of anchors:", length(merged$gene)) )
merged <- merged[f.feat,]
print(head(merged))
print(paste("Length of filtered anchors:", length(merged$gene)) )



merged$name <- ifelse(merged$p_val_adj.sc < 10^-50 & merged$p_val_adj.bulk < 0.05 ,merged$gene, NA )

#select markers to annotate
own.gene <- c("HLA-DPA1", "HLA-DRB5","ETV1", "CD80") #"CD109","CD99",,"ITGB2" "CD8A","CD27","CD8B", "IL13",,"GZMH", "ITM2A""GZMH","IRF5", "KLRC2"
exhaustion.gene <- c("HLA-DOA", "HLA-DQB1", "HLA-DQA1","HLA-DRB1", "CD74",
                     "CD70","PDCD1", "CD86", "CTLA4", "IRF5", "CTSH", "NPDC1") #
mem.gene <- c("LEF1","IL12RB2", "CCR7",  "IL7R", "TNFRSF18", "TNFRSF4", "TRDC", "TCF7","TRGC1")
selection <- c(own.gene, exhaustion.gene, mem.gene)


selection[!(selection %in% merged[merged$sig, "gene"])]

merged$name <- ifelse(merged$gene %in% selection ,merged$gene, NA )


#define if dots should be shown based on threshold
merged$show <-  ifelse(merged$p_val_adj.sc < 1e-25 , "yes", "no") #
merged$show <- factor(merged$show, levels = c("no", "yes"))

selection[!(selection %in% merged[merged$show == "yes", "gene"])]


merged$sig <-  ifelse(merged$p_val_adj.sc < 1e-50 & merged$p_val_adj.bulk < 0.05, TRUE, FALSE) #


#define value for color scaling 
merged$color_scale <- (merged$avg_log2FC.sc + merged$avg_log2FC.bulk)
merged[merged$color_scale>1.5,"color_scale"] <- 1.5
merged[merged$color_scale<= -1.5,"color_scale"] <- -1.5

write.csv(merged, paste0("Fig/DGE_results/DGE_Reactive_correlation_preREP_IntermediateProduct.csv"), row.names = F)

pdf("Fig/DGE_results/DGE_Reactive_correlation_PreREP.pdf", width = 5, height = 5)
plt2<- ggplot(merged) + 
  geom_point(data = merged[merged$sig == FALSE &merged$show == "yes",],aes(x=avg_log2FC.sc, y=avg_log2FC.bulk, color = color_scale), size = 1)+ 
  scale_colour_gradient2(
    low = "#7E8CBE",
    mid = "gray90",
    high = "#F8766D",
    midpoint = 0,
    na.value = "grey50",
    limits=c(-1.5,1.5)
  ) +
  new_scale_colour() +
  geom_point(data = merged[merged$sig == TRUE & merged$show == "yes",], aes(x=avg_log2FC.sc, y=avg_log2FC.bulk, color = color_scale),  size = 1) + 
  scale_colour_gradient2(
    low = "#3C4DA0",
    mid = "gray90",
    high = "#CD1F2D",
    midpoint = 0,
    na.value = "grey50",
    limits=c(-1.5,1.5))+
  theme_classic()+
  geom_hline(yintercept = 0, linetype = "dotted", color = "black")+
  geom_vline(xintercept = 0, linetype = "dotted", color = "black")+
  geom_label_repel(aes(label = name,x=avg_log2FC.sc, y=avg_log2FC.bulk),
                   parse = F,
                   size = 3,
                   box.padding = 0.25,
                   segment.color = "black",
                   min.segment.length = 0, max.time = 5,
                   show.legend = FALSE, max.overlaps = 20, point.padding = 0, label.padding = 0.15) + NoLegend() + 
  xlab("single-cell log2FC") + 
  ylab("patient averaged log2FC") + 
  xlim(-3.2, 2.5)+
  ylim(-1.5, 1.2)+
  theme(
    legend.position = "none",
    axis.title = element_text(size = 15,color="black"),
    axis.text = element_text(size = 12,color="black")
  )

plt2
dev.off()





## -----------------------------------------------------------------------------------------------------------
#Plotting the graphs used in manuscript for ExpandedTILs (TIL product)


# loading the single cell data for ExpandedTILs
markers  <- read.csv("DGE_results/dge_Expanded_TILs.csv")

#now for pseudo bulk results
markers_bulk  <- read.csv("DGE_results/dge_Expanded_TILs_bulk_new.csv")


markers_bulk$p_val <- markers_bulk$pvalue
markers_bulk$avg_log2FC <- markers_bulk$log2FoldChange
markers_bulk$p_val_adj <- markers_bulk$padj

merged <- left_join(markers, markers_bulk, by = "gene", suffix = c(".sc", ".bulk") )

#filter the merged table
load("../data/Zheng_Japrin/exclude.gene.misc.human.v4.RData")
f.feat <- !(merged$gene %in% all.gene.ignore.df[["seu.id"]]) &
  !(grepl("^(RP[LS]|Rp[ls])",merged$gene,perl=T)) &
  !(grepl("^MT-",merged$gene ,perl=T)) &
  !(merged$gene %in% c("MALAT1", "NEAT1"))
print(paste("Length of anchors:", length(merged$gene)) )
merged <- merged[f.feat,]
print(head(merged))
print(paste("Length of filtered anchors:", length(merged$gene)) )




merged$name <- ifelse(merged$p_val_adj.sc < 10^-50 & merged$p_val_adj.bulk < 0.05 ,merged$gene, NA )

#select markers to annotate
own.gene <- c("ETV1", "GPR141", "IRF5", "HSF4", "TNFRSF9") #, "NCR1"
exhaustion.gene <- c("TNFSF4","HLA-DQA1","HLA-DPA1","CD70","TOX", "TNFRSF9", "PDCD1", "EPDR1", "RGS1")
mem.gene <- c("TRDC","IL12RB2", "TNFRSF18", "IL4",  "CD300A",    "IL7R", "GPR55")
selection <- c(own.gene, exhaustion.gene, mem.gene)

selection[!(selection %in% merged[merged$sig, "gene"])]
merged$name <- ifelse(merged$gene %in% selection ,merged$gene, NA )


#define if dots should be shown based on threshold
merged$show <-  ifelse(merged$p_val_adj.sc < 1e-25 , "yes", "no") #
merged$show <- factor(merged$show, levels = c("no", "yes"))

selection[!(selection %in% merged[merged$show == "yes", "gene"])]


merged$sig <-  ifelse(merged$p_val_adj.sc < 1e-50 & merged$p_val_adj.bulk < 0.05, TRUE, FALSE) #


#define value for color scaling 
merged$color_scale <- (merged$avg_log2FC.sc + merged$avg_log2FC.bulk)
merged[merged$color_scale>1.5,"color_scale"] <- 1.5
merged[merged$color_scale<= -1.5,"color_scale"] <- -1.5

write.csv(merged, paste0("Fig/DGE_results/DGE_Reactive_correlation_TIL_Product.csv"), row.names = F)

pdf("Fig/DGE_results/DGE_Reactive_correlation_ExpandedTILs.pdf", width = 5, height = 5)
plt3<- ggplot(merged) + 
  geom_point(data = merged[merged$sig == FALSE &merged$show == "yes",],aes(x=avg_log2FC.sc, y=avg_log2FC.bulk, color = color_scale), size = 1)+ 
  scale_colour_gradient2(
    low = "#7E8CBE",
    mid = "gray90",
    high = "#F8766D",
    midpoint = 0,
    na.value = "grey50",
    limits=c(-1.5,1.5)
  ) +
  new_scale_colour() +
  geom_point(data = merged[merged$sig == TRUE & merged$show == "yes",], aes(x=avg_log2FC.sc, y=avg_log2FC.bulk, color = color_scale),  size = 1) + 
  scale_colour_gradient2(
    low = "#3C4DA0",
    mid = "gray90",
    high = "#CD1F2D",
    midpoint = 0,
    na.value = "grey50",
    limits=c(-1.5,1.5))+
  theme_classic()+
  geom_hline(yintercept = 0, linetype = "dotted", color = "black")+
  geom_vline(xintercept = 0, linetype = "dotted", color = "black")+
  geom_label_repel(mapping = aes(label = name, x=avg_log2FC.sc, y=avg_log2FC.bulk), 
                   parse = F,
                   size = 3,
                   box.padding = 0.25,
                   segment.color = "black",
                   min.segment.length = 0,max.time = 1,
                   show.legend = FALSE, max.overlaps = 10, point.padding = 0, label.padding = 0.15) + NoLegend() + 
  xlab("single-cell log2FC") + 
  ylab("patient averaged log2FC") +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 15,color="black"),
    axis.text = element_text(size = 12,color="black")
  )+ 
  xlim(-3, 3)+
  ylim(-1.8, 1.8)
plt3
dev.off()





## -----------------------------------------------------------------------------------------------------------
#Plotting the graphs used in manuscript for PBMC_7dpt (PBMC 7 days post transfer)


# loading the single cell data for PBMC
markers  <- read.csv("DGE_results/dge_PBMC_7dpt.csv")


#now for pseudo bulk results
markers_bulk  <- read.csv("DGE_results/dge_PBMC_7dpt_bulk_new.csv")

markers_bulk$p_val <- markers_bulk$pvalue
markers_bulk$avg_log2FC <- markers_bulk$log2FoldChange
markers_bulk$p_val_adj <- markers_bulk$padj

merged <- left_join(markers, markers_bulk, by = "gene", suffix = c(".sc", ".bulk") )

#filter the merged table
load("../data/Zheng_Japrin/exclude.gene.misc.human.v4.RData")
f.feat <- !(merged$gene %in% all.gene.ignore.df[["seu.id"]]) &
  !(grepl("^(RP[LS]|Rp[ls])",merged$gene,perl=T)) &
  !(grepl("^MT-",merged$gene ,perl=T)) &
  !(merged$gene %in% c("MALAT1", "NEAT1"))
print(paste("Length of anchors:", length(merged$gene)) )
merged <- merged[f.feat,]
print(head(merged))
print(paste("Length of filtered anchors:", length(merged$gene)) )


merged$name <- ifelse(merged$p_val_adj.sc < 10^-10 & merged$p_val_adj.bulk < 0.05 ,merged$gene, NA )

#select markers to annotate
own.gene <- c("CCL23",  "RGS18")
exhaustion.gene <- c("ETV1","CD33","IL7","XCL2","EOMES","CTSG","IL18RAP", "HLA-DQA1", "ITM2A", "ETV7", "TOX", "RBPMS", "CXCR6")
mem.gene <- c("ITM2C", "GPR183", "ZNF683", "TNFRSF9", "IRF5", "KLF3")
selection <- c(own.gene, exhaustion.gene, mem.gene)
merged$name <- ifelse(merged$gene %in% selection ,merged$gene, NA )


selection[!(selection %in% merged[merged$sig, "gene"])]
merged$name <- ifelse(merged$gene %in% selection ,merged$gene, NA )


#define if dots should be shown based on threshold
merged$show <-  ifelse(merged$p_val_adj.sc < 1e-10 , "yes", "no") #
merged$show <- factor(merged$show, levels = c("no", "yes"))

selection[!(selection %in% merged[merged$show == "yes", "gene"])]


merged$sig <-  ifelse(merged$p_val_adj.sc < 1e-50 & merged$p_val_adj.bulk < 0.05, TRUE, FALSE) #


#define value for color scaling 
merged$color_scale <- (merged$avg_log2FC.sc + merged$avg_log2FC.bulk)
merged[merged$color_scale>1.5,"color_scale"] <- 1.5
merged[merged$color_scale<= -1.5,"color_scale"] <- -1.5

write.csv(merged, paste0("Fig/DGE_results/DGE_Reactive_correlation_PBMC_7dpt.csv"), row.names = F)

pdf("Fig/DGE_results/DGE_Reactive_correlation_PBMC.pdf", width = 5, height = 5)
plt4<-   ggplot(merged) + 
  geom_point(data = merged[merged$sig == FALSE &merged$show == "yes",],aes(x=avg_log2FC.sc, y=avg_log2FC.bulk, color = color_scale), size = 1)+ 
  scale_colour_gradient2(
    low = "#7E8CBE",
    mid = "gray90",
    high = "#F8766D",
    midpoint = 0,
    na.value = "grey50",
    limits=c(-1.5,1.5)
  ) +
  new_scale_colour() +
  geom_point(data = merged[merged$sig == TRUE & merged$show == "yes",], aes(x=avg_log2FC.sc, y=avg_log2FC.bulk, color = color_scale),  size = 1) + 
  scale_colour_gradient2(
    low = "#3C4DA0",
    mid = "gray90",
    high = "#CD1F2D",
    midpoint = 0,
    na.value = "grey50",
    limits=c(-1.5,1.5))+
  theme_classic()+
  geom_hline(yintercept = 0, linetype = "dotted", color = "black")+
  geom_vline(xintercept = 0, linetype = "dotted", color = "black")+
  geom_label_repel(aes(label = name,x=avg_log2FC.sc, y=avg_log2FC.bulk),
                   parse = F,
                   size = 3,
                   box.padding = 0.25,
                   segment.color = "black",
                   min.segment.length = 0,max.time = 1,
                   show.legend = FALSE, max.overlaps = 20, point.padding = 0, label.padding = 0.15) + NoLegend() + 
  xlab("single-cell log2FC") + 
  ylab("patient averaged log2FC") +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 15,color="black"),
    axis.text = element_text(size = 12,color="black")
  )+ ylim(-1.8, 2.3)
# xlim(-5, 3.5) #+

plt4
dev.off()










## -----------------------------------------------------------------------------------------------------------
#Plotting the graphs used in manuscript for post-ACT 


# loading the single cell data for post-ACT
markers  <- read.csv("DGE_results/dge_Post.csv")


#now for pseudo bulk results
markers_bulk  <- read.csv("DGE_results/dge_Post_bulk_new.csv")
markers_bulk$p_val <- markers_bulk$pvalue
markers_bulk$avg_log2FC <- markers_bulk$log2FoldChange
markers_bulk$p_val_adj <- markers_bulk$padj


merged <- left_join(markers, markers_bulk, by = "gene", suffix = c(".sc", ".bulk") )

#filter the merged table
load("../data/Zheng_Japrin/exclude.gene.misc.human.v4.RData")
f.feat <- !(merged$gene %in% all.gene.ignore.df[["seu.id"]]) &
  !(grepl("^(RP[LS]|Rp[ls])",merged$gene,perl=T)) &
  !(grepl("^MT-",merged$gene ,perl=T)) &
  !(merged$gene %in% c("MALAT1", "NEAT1"))
print(paste("Length of anchors:", length(merged$gene)) )
merged <- merged[f.feat,]
print(head(merged))
print(paste("Length of filtered anchors:", length(merged$gene)) )

merged$name <- ifelse(merged$p_val_adj.sc < 10^-50 & merged$p_val_adj.bulk < 0.05 ,merged$gene, NA )


#select markers to annotate
own.gene <- c( "KIR2DL3", "KIR2DL4", "HAVCR2", "KLRC2", "ADGRG1") # "KLRC3", "KLRC4", "LAG3", "TNFSF9", "ITGA2","KLRC2","TNFSF4",, "TGIF1", "TOX2"
exhaustion.gene <- c("CXCL13","ETV1","ENTPD1","TOX", "SNX9", "DUSP4", "CD27", "PDCD1",  "CTLA4")
mem.gene <- c("SELL","IL7R","TCF7", "LEF1","S1PR1", "ZNF683", "KLF3", "TNFSF14")
selection <- c(own.gene, exhaustion.gene, mem.gene)
merged$name <- ifelse(merged$gene %in% selection ,merged$gene, NA )

#define if dots should be shown based on threshold
merged$show <-  ifelse(merged$p_val_adj.sc < 1e-25 , "yes", "no") #
merged$show <- factor(merged$show, levels = c("no", "yes"))

selection[!(selection %in% merged[merged$show == "yes", "gene"])]


merged$sig <-  ifelse(merged$p_val_adj.sc < 1e-50 & merged$p_val_adj.bulk < 0.05, TRUE, FALSE) #
selection[!(selection %in% merged[merged$sig, "gene"])]

#define value for color scaling 
merged$color_scale <- (merged$avg_log2FC.sc + merged$avg_log2FC.bulk)
merged[merged$color_scale>1.5,"color_scale"] <- 1.5
merged[merged$color_scale<= -1.5,"color_scale"] <- -1.5

write.csv(merged, paste0("Fig/DGE_results/DGE_Reactive_correlation_postACT.csv"), row.names = F)


pdf("Fig/DGE_results/DGE_Reactive_correlation_Post.pdf", width = 5, height = 5)
plt5<- ggplot(merged) + 
  geom_point(data = merged[merged$sig == FALSE &merged$show == "yes",],aes(x=avg_log2FC.sc, y=avg_log2FC.bulk, color = color_scale), size = 1)+ 
  scale_colour_gradient2(
    low = "#7E8CBE",
    mid = "gray90",
    high = "#F8766D",
    midpoint = 0,
    na.value = "grey50",
    #limits=c(-1.5,1.5)
  ) +
  new_scale_colour() +
  geom_point(data = merged[merged$sig == TRUE & merged$show == "yes",], aes(x=avg_log2FC.sc, y=avg_log2FC.bulk, color = color_scale),  size = 1) + 
  scale_colour_gradient2(
    low = "#3C4DA0",
    mid = "gray90",
    high = "#CD1F2D",
    midpoint = 0,
    na.value = "grey50",
    limits=c(-1.5,1.5))+
  theme_classic()+
  geom_hline(yintercept = 0, linetype = "dotted", color = "black")+
  geom_vline(xintercept = 0, linetype = "dotted", color = "black")+
  geom_label_repel(aes(label = name,x=avg_log2FC.sc, y=avg_log2FC.bulk),
                   parse = F,
                   size = 3,
                   box.padding = 0.25,
                   segment.color = "black",
                   min.segment.length = 0, max.time = 1,
                   show.legend = FALSE, max.overlaps = 20, point.padding = 0, label.padding = 0.15) + NoLegend() + 
  xlab("single-cell log2FC") + 
  ylab("patient averaged log2FC") + 
  xlim(-3.5, 3.5) +
  ylim(-3, 3) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 15,color="black"),
    axis.text = element_text(size = 12,color="black")
  )

plt5
dev.off()



pdf("Fig/DGE_results/DGE_Reactive_correlation_all.pdf", width = 25, height = 5)
ggpubr::ggarrange(plt1,plt2,plt3,plt4,plt5, ncol = 5, nrow = 1)
dev.off()







