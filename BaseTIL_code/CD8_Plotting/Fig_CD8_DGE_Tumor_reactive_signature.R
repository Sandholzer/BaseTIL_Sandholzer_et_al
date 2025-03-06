## -----------------------------------------------------------------------------------------------------------
#load dependencies

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(stringr)
  library(ggplot2)
  library(ggalluvial)
  library(ggrepel)
  library(ggpubr)
  library(readr)
  library(ggeasy)
  library(purrr)
  library(ComplexHeatmap)
  library(DESeq2)
  library(EnhancedVolcano)
  library(tibble)
})

set.seed(1234567)
setwd("~/BaseTIL_code/CD8_Plotting")

gex.CD8 <- readRDS(file = "../saveFiles/CD8_annotated.rds")



## -----------------------------------------------------------------------------------------------------------
# Differential expression analysis comparing reactive to bystander T cells in all tumor derived samples



reactive_clones <-subset(gex.CD8@meta.data, subset = Type %in% c("Tumor", "Rebiopsy1", "Rebiopsy2") & selfreactive == TRUE)
reactive_clones <-  reactive_clones[!is.na(reactive_clones$ntb),]


gex.CD8$reactive_clones <- "None"
gex.CD8$reactive_clones <- gex.CD8$ntb %in%  reactive_clones$ntb


#remove samples were too low TR were identified
tum.CD8 <- subset(gex.CD8, Type %in% c("Tumor", "Rebiopsy1", "Rebiopsy2") & Patient %in% c("UPN001", "UPN003", "UPN006", "UPN008", "UPN009", "UPN011")) #"UPN002",
tum.CD8@meta.data[tum.CD8@meta.data$Sample_Name == "UPN006 Rebiopsy2", "Type"] <- "Rebiopsy1"
tum.CD8@meta.data[tum.CD8@meta.data$Sample_Name == "UPN006 Rebiopsy2", "Sample_Name"] <- "UPN006 Rebiopsy1" # consider UPN006 rebiopsies as one sample due to low cell numbers and reduction of resulting variance
tum.CD8 <- subset(tum.CD8, Sample_Name %in% c("UPN008 Tumor","UPN008 Rebiopsy1", "UPN001 Rebiopsy2", "UPN011 Rebiopsy1", "UPN001 Rebiopsy1"), invert = T) 

#only use T cells defined as "selfreactive" = tumor-reactive against the tumor lesions they where detected in 
Idents(tum.CD8) <- tum.CD8$selfreactive


## -----------------------------------------------------------------------------------------------------------
# DGE considering single cells


markers <- FindMarkers(tum.CD8, ident.1 = "TRUE",
                       ident.2 = "FALSE", logfc.threshold = 0) 
markers$gene <- rownames(markers)
write.csv(markers, paste0("DGE_results/dge_tumorReactive_signature.csv"), row.names = F)

# Filtering for non relevant genes
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




## -----------------------------------------------------------------------------------------------------------
# DGE considering patient replicates

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


# Check the coefficients for the comparison
resultsNames(dds)

# Generate results object
res <- results(dds, 
               name = "reactive_clones_TRUE_vs_FALSE",
               alpha = 0.1)


markers_bulk <- res %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  arrange(padj)


markers_bulk <- markers_bulk[!is.na(markers_bulk$pvalue), ] 

write.csv(markers_bulk, paste0("DGE_results/dge_tumorReactive_signature_bulk.csv"), row.names = F)

# Filtering for non relevant genes
load("../data/Zheng_Japrin/exclude.gene.misc.human.v4.RData")

f.feat <- !(markers_bulk$gene %in% all.gene.ignore.df[["seu.id"]]) &
  !(grepl("^(RP[LS]|Rp[ls])",markers_bulk$gene,perl=T)) &
  !(grepl("^MT-|^MTRN|^MTAT|^MTND|^MRP",markers_bulk$gene,perl=T)) &
  !(grepl("^RP[0-9]+-",markers_bulk$gene,perl=T)) &
  !(markers_bulk$gene %in% c("MALAT1", 'XIST', 'NEAT1'))
markers_bulk.flt <- markers_bulk[f.feat,]
print(head(markers_bulk.flt))
print(paste("Length of anchors:", nrow(markers_bulk)) )
print(paste("Length of filtered anchors:", nrow(markers_bulk.flt)) )




## -----------------------------------------------------------------------------------------------------------
# Plotting results in vulcano plots

pdf("DGE_results/DGE_react_tumorReactive_signature.pdf", width = 7, height = 9)
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

EnhancedVolcano(markers.flt, x = "avg_log2FC", y = "p_val_adj", lab = markers.flt$gene)

markers_bulk.flt %>%
  dplyr::rename(FC = log2FoldChange) %>%
  dplyr::mutate(name = ifelse(gene %in% c(markers_bulk.flt$gene[markers_bulk.flt$pvalue < 0.01 & abs(markers_bulk.flt$log2FoldChange) > 0.5]), gene, NA)) %>%
  dplyr::mutate(sig = ifelse(pvalue < 0.01 & abs(FC) > 0.15, "yes", "no")) %>%
  ggplot(aes(FC, -log10(pvalue))) +
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

EnhancedVolcano(markers_bulk.flt, x = "log2FoldChange", y = "pvalue", lab = markers_bulk.flt$gene, pCutoff = 0.05, FCcutoff = 0.25)


dev.off()



## -----------------------------------------------------------------------------------------------------------
# merge both outputs to compare in correlation plot

markers_bulk$p_val <- markers_bulk$pvalue
markers_bulk$avg_log2FC <- markers_bulk$log2FoldChange
markers_bulk$p_val_adj <- markers_bulk$padj

merged <- left_join(markers, markers_bulk, by = "gene", suffix = c(".sc", ".bulk") )

load("../data/Zheng_Japrin/exclude.gene.misc.human.v4.RData")
f.feat <- !(merged$gene %in% all.gene.ignore.df[["seu.id"]]) &
  !(grepl("^(RP[LS]|Rp[ls])",merged$gene,perl=T)) &
  !(grepl("^MT-",merged$gene ,perl=T)) &
  !(merged$gene %in% c("MALAT1", "NEAT1"))
print(paste("Length of anchors:", length(merged$gene)) )
merged <- merged[f.feat,]
print(head(merged))
print(paste("Length of filtered anchors:", length(merged$gene)) )

write.csv(merged, paste0("Fig/DGE_results/dge_tumorReactive_signature_merged.csv"), row.names = F)


merged$name <- ifelse(merged$p_val_adj.sc < 1e-100 & merged$p_val.bulk < 0.05 ,merged$gene, NA )


merged$sig <-  ifelse(merged$p_val_adj.sc < 1e-100 & merged$p_val.bulk < 0.05, "yes", "no") #
merged$sig <- factor(merged$sig, levels = c("no", "yes"))

merged$color_scale <- (merged$avg_log2FC.sc + merged$avg_log2FC.bulk)
merged[merged$color_scale>1.5,"color_scale"] <- 1.5
merged[merged$color_scale<= -1.5,"color_scale"] <- -1.5

pdf("DGE_results/DGE_Reactive_signature_correlation.pdf", width = 5, height = 5)
plt1<- ggplot(data = merged[merged$sig == "yes",], aes(x=avg_log2FC.sc, y=avg_log2FC.bulk)) + 
  geom_point(aes(color = color_scale), size = 0.5) + 
  theme_classic()+
  geom_hline(yintercept = 0, linetype = "dotted", color = "black")+
  geom_vline(xintercept = 0, linetype = "dotted", color = "black")+
  geom_text_repel(aes(label = name),
                  parse = F,
                  size = 3,
                  box.padding = 0.25,
                  segment.color = "black",
                  show.legend = FALSE, max.overlaps = 50)+ scale_colour_gradient2(
                    low = "#3C4DA0",
                    mid = "gray90",
                    high = "#CD1F2D",
                    midpoint = 0,
                    na.value = "grey50",
                    limits=c(-1.5,1.5)
                  ) + NoLegend() + 
  xlab("single-cell log2FC") + 
  ylab("patient averaged log2FC") + xlim(-5, 5)
plt1
dev.off()




