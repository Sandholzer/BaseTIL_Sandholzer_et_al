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
  library(escape)
  library(ggrepel)
  library(scater)
  library(Seurat)
  library(cowplot)
  library(reticulate)
  library(SeuratWrappers)
  library(SeuratDisk)
  library(RColorBrewer)
  library(circlize)
  library(ComplexHeatmap)
  library(tidyverse)
  library(escape)
})

set.seed(12345678)


setwd("~/BaseTIL_code/CD8_Plotting/")

react.tum <-  readRDS(file = "~/BaseTIL_code/saveFiles/CD8_tumor_reactive_clustering.rds")



react.tum$prePost <- ifelse(react.tum$Type %in% c("Rebiopsy1", "Rebiopsy2"), "post-ACT", "pre-ACT") 
react.tum$prePost <- factor(react.tum$prePost, levels = c("pre-ACT", "post-ACT"))

#abbreviate response definition
react.tum$RNR <- ifelse(react.tum$Response %in% c("PR", "CR"), "Rs", "NRs") 
react.tum$RNR <- factor(react.tum$RNR, levels = c("NRs", "Rs"))


#annotate clusters
react.tum$seurat_clusters <- react.tum$RNA_snn_res.0.4
Idents(react.tum) <- react.tum$seurat_clusters
react.tum <-
  RenameIdents(
    react.tum,
    `1` = "Ttex",
    `2` = "GZMK+ Tex",
    `3` = "ITGAE+ Tex",
    `4` = "CXCL13+ Tex",
    `5` = "Tn",
    `6` = "Tpex",
    `7` = "Tpex")


DimPlot(react.tum, label = T)

react.tum$TR_anno <- Idents(react.tum)


react.tum$TR_anno <- factor(
  react.tum$TR_anno,
  levels = rev(c("Ttex",
             "GZMK+ Tex",
             "CXCL13+ Tex",
             "ITGAE+ Tex",
             "Tpex",
             "Tn")
))


#save annotated file
saveRDS(react.tum, "~/BaseTIL_code/saveFiles/CD8_tumor_reactive_annotated.rds")

## -----------------------------------------------------------------------------------------------------------
#save as h5ad for RNA velocity analysis


#creates an v4 assay since v5 is not supported
react.tum[["RNA3"]] <- as(object = react.tum[["RNA"]], Class = "Assay")
DefaultAssay(react.tum) <- "RNA3"
react.tum[["RNA"]] <- NULL
react.tum <- RenameAssays(object = react.tum, RNA3 = 'RNA')

VariableFeatures(react.tum)

SaveH5Seurat(react.tum, filename = "~/BaseTIL_code/saveFiles/CD8_tumor_reactive.h5Seurat", overwrite = T)
Convert("~/BaseTIL_code/saveFiles/CD8_tumor_reactive.h5Seurat", dest = "h5ad", overwrite = T)





## -----------------------------------------------------------------------------------------------------------
#plotting

Colors <- RColorBrewer::brewer.pal(n=12, "Paired")

pdf("Fig/CD8_TR_barplot_frequency.pdf", height = 3, width = 3.5)
dittoSeq::dittoBarPlot(react.tum, group.by = "prePost", var = "TR_anno", split.by = "RNR", retain.factor.levels = T, color.panel = Colors) &
  ggtitle(NULL) & xlab(NULL) & theme(axis.text.x = element_text(angle = 45,  hjust=1))
dev.off()


pdf("Fig/CD8_TR_Umap.pdf", height = 4, width = 5)
DimPlot(react.tum, reduction = "umap", group.by = "TR_anno", label = T, raster = F, repel = T, pt.size = 0.5, cols = Colors, order = F,
        label.size = 5, label.box = T) & NoLegend()&
  ggeasy::easy_remove_axes() & ggtitle(NULL)
dev.off()


pdf("Fig/CD8_TR_Umap_split.pdf", height = 4, width = 10)
DimPlot(react.tum, reduction = "umap", group.by = "TR_anno", label = F, raster = F, split.by = "prePost", pt.size = 0.5, cols = Colors,
        label.size = 4, label.box = F)  & 
  ggeasy::easy_remove_axes()& ggtitle(NULL)
dev.off()

pdf("Fig/TR_Umap_marker.pdf", height = 6.5, width = 10)
Nebulosa::plot_density(react.tum, c("IL7R",  "SELL","TCF7",
                                    "GZMK","ITGAE", "CXCL13",
                                    "PDCD1","HAVCR2","KIR2DL4")) & 
                                      ggeasy::easy_remove_axes() & 
                                      ggeasy::easy_remove_legend_title() &
                                      theme(axis.text.x=element_blank(),
                                            axis.text.y=element_blank(),
                                            plot.title = element_text(hjust = 0.5))


Nebulosa::plot_density(react.tum, c("ADGRG1", "ZNF683",  "SNX9",
                                    "TNFRSF9", "XCL1", "XCL2",
                                    "ENTPD1","ITGA1","ITGB1")) & 
                                      ggeasy::easy_remove_axes() & 
                                      ggeasy::easy_remove_legend_title() &
                                      theme(axis.text.x=element_blank(),
                                            axis.text.y=element_blank(),
                                            plot.title = element_text(hjust = 0.5))


Nebulosa::plot_density(react.tum, c("NFKB1", "ZNF683",  "NFKB2",
                                         "NR4A2", "TNFRSF18", "NFKBID",
                                         "LAIR2","HOPX","NCR3")) & 
  ggeasy::easy_remove_axes() & 
  ggeasy::easy_remove_legend_title() &
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        plot.title = element_text(hjust = 0.5))


dev.off()




## -----------------------------------------------------------------------------------------------------------
#DGE for heatmap to verify cluster annotation
aggro <- AggregateExpression(react.tum, return.seurat = T, group.by = "TR_anno")


aggro$TR_anno <- factor(
  aggro$TR_anno,
  levels = rev(c("Ttex",
                 "GZMK+ Tex",
                 "CXCL13+ Tex",
                 "ITGAE+ Tex",
                 "Tpex",
                 "Tn")
  ))



markers_clusters <- FindAllMarkers(react.tum, only.pos = T)

write.csv(markers_clusters, "Graphs/TR_clustering_markers.csv")


load("~/BaseTIL_code/data/Zheng_Japrin/exclude.gene.misc.human.v4.RData")

f.feat <- !(markers_clusters$gene %in% all.gene.ignore.df[["seu.id"]]) &
  !(grepl("^(RP[LS]|Rp[ls])",markers_clusters$gene,perl=T)) &
  !(grepl("^MT-",markers_clusters$gene,perl=T)) &
  !(markers_clusters$gene %in% c("MALAT1", "Malat1", "NEAT1"))
gex.features.flt <- markers_clusters[f.feat,]
print(head(gex.features.flt))
print(paste("Length of anchors:", length(markers_clusters$gene)) )
print(paste("Length of filtered anchors:", length(gex.features.flt$gene)) )


gex.features.flt$cluster <- factor(
  gex.features.flt$cluster,
  levels = rev(c("Ttex",
                 "GZMK+ Tex",
                 "CXCL13+ Tex",
                 "ITGAE+ Tex",
                 "Tpex",
                 "Tn")
  ))

gex.features.flt %>%
  group_by(cluster) %>%
  #dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

top_genes <- top10$gene
top_genes[top_genes == "LINC02341"] <- "CXCL13"

mat<- aggro@assays[["RNA"]]$scale.data[c(top_genes), ] %>% as.matrix()

quantile(mat, c(0.05, 0.95))


pdf("Fig/TR_Cluster_annotation_heatmap.pdf", height = 10, width = 3.5)
Heatmap(mat, name = "Z score",  
        column_split = aggro@meta.data$TR_anno,
        show_column_dend = T,
        cluster_columns = F,
        cluster_column_slices = TRUE,
        cluster_rows = F,
        column_title_gp = gpar(fontsize = 12),
        column_gap = unit(0.75, "mm"),
        show_row_dend = FALSE,
        col=colorRamp2(c(-1.5, 0, 1.7), c("#3C4DA0", "white", "#CD1F2D")),
        row_names_gp = gpar(fontsize = 12),
        column_title_rot = 45,
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = rev(c(Colors[1:6]))))),
        show_column_names = F,
        use_raster = TRUE,
        raster_quality = 4,
        row_names_side = "left",
)
dev.off()


## -----------------------------------------------------------------------------------------------------------
#highlight expansion state on UMAP

pdf("Fig/CD8_TR_Expansionstate_split_Type.pdf",  height = 4, width = 7)
DimPlot(react.tum, group.by = "cloneType_percent", raster = F,  pt.size = 0.5) &  scale_color_manual(
  "Clonetype Group",
  values = c(
    "Small (0 < X <= 0.001)" =  "#0D0887FF",
    "Medium (0.001 < X <= 0.01)" = "#9C179EFF",
    "Large (0.01 < X <= 0.1)" = "#ED7953FF",
    "Hyperexpanded (0.1 < X <= 1)" = "#F0F921FF"
  ), na.value = "gray90") & easy_remove_axes()
dev.off()


## -----------------------------------------------------------------------------------------------------------
#Signature scoring for matched clones of pre-ACT and post-ACT


#load signatures
published_signatures <- read_csv("~/BaseTIL_code/data/Published_signatures.csv")
gene.sets <- published_signatures %>%
  pivot_longer(cols = everything(), names_to = "Column") %>%
  group_by(Column) %>%
  summarize(Values = list(value)) %>%
  deframe()

gene.sets<- lapply(gene.sets, function(x) x <- x[!is.na(x)])

selection <- names(gene.sets)
selection <- c("Oliveira.TTE", "Li.CD8.DYS",
                "Jansen_Stem like", "Jansen_Term diff" 
               )
gene.sets <- gene.sets[selection]


preTCR <- react.tum@meta.data[react.tum$Type=="Tumor", "ntb"]
postTCR <- react.tum@meta.data[react.tum$Type %in% c("Rebiopsy1", "Rebiopsy2"), "ntb"]
react.shared <- subset(react.tum,  cells =rownames(react.tum@meta.data[react.tum$ntb %in% preTCR & react.tum$ntb %in% postTCR,]))


react.shared$Type <- as.character(react.shared$Type)
react.shared$Type <- ifelse(react.shared$Type %in% c("Rebiopsy1", "Rebiopsy2"), "post-ACT", react.shared$Type)


react.shared$reactivity <- NA
react.shared@meta.data[react.shared$Tumor_reactive==TRUE & react.shared$Post_reactive==FALSE,"reactivity"] <- "pre-ACT"
react.shared@meta.data[react.shared$Tumor_reactive==FALSE & react.shared$Post_reactive==TRUE,"reactivity"] <- "post-ACT"
react.shared@meta.data[react.shared$Tumor_reactive==TRUE & react.shared$Post_reactive==TRUE,"reactivity"] <- "both"


aggro_seur <- AggregateExpression(react.shared, return.seurat = T, group.by = c("ntb", "Type", "Patient", "reactivity"))
aggro_seur$Type_post <- ifelse(aggro_seur$Type %in% c("Rebiopsy1", "Rebiopsy2"), "post-ACT", aggro_seur$Type)

# Calculate signature scores
aggro_seur_new <- AddModuleScore_UCell(aggro_seur, gene.sets, name = "")


ES2 <- data.frame(aggro_seur_new[[]])
ES2$orig.ident <- NULL
ES2$Type_post <- NULL
ES2$Type <- ifelse(ES2$Type %in% c("Rebiopsy1", "Rebiopsy2"), "post-ACT", ES2$Type)
plt <- list()
vector <- ES2 %>% select (-c(ntb, Type, Patient, reactivity)) %>% colnames()


# Plot it separately for reponder and non-responder
pdf("Fig/Signature_prePost_reactive.pdf", height = 4, width = 3)
for (iterater in vector) {
  plot<-  ES2[ES2$Patient %in% c("UPN008", "UPN009", "UPN011", "UPN006"),] %>%
  ggplot()+
    aes(x = factor(Type,levels = c("Tumor", "post-ACT"), labels = c("pre-ACT", "post-ACT")),
        y = .data[[iterater]],
        #color = Patient,
        group = ntb) +
    geom_point(size=0.75) +
    geom_line(linewidth=0.25) +
    xlab(NULL)+
    ylab("Signature Score") +
    ggtitle(iterater) +
    theme_classic() + theme(text = element_text(size=20),
                            axis.text.x = element_text(angle = 45, hjust=1))+
    stat_compare_means(comparisons = list(c("pre-ACT", "post-ACT")), label = "p.signif", method = "t.test", paired = T, size = 8)+ 
    scale_y_continuous(expand = expansion(mult = 0.1))
  print(plot)


  plot<-  ES2[ES2$Patient %in% c("UPN001"),] %>%
  ggplot()+
    aes(x = factor(Type,levels = c("Tumor", "post-ACT"), labels = c("pre-ACT", "post-ACT")),
        y = .data[[iterater]],
        #color = Patient,
        group = ntb) +
    geom_point(size=0.75) +
    geom_line(linewidth=0.25) +
    xlab(NULL)+
    ylab("Signature Score") +
    ggtitle(iterater) +
    theme_classic() + theme(text = element_text(size=20),
                            axis.text.x = element_text(angle = 45, hjust=1))+
    stat_compare_means(comparisons = list(c("pre-ACT", "post-ACT")), label = "p.signif", method = "t.test", paired = T, size = 8) +
    scale_y_continuous(expand = expansion(mult = 0.1))
  print(plot)


}
dev.off()


## -----------------------------------------------------------------------------------------------------------
#Differential expression analysis between pre-ACT and post-ACT tumor reactive clones

Idents(react.tum) <- react.tum$prePost


markers_R <- FindMarkers(subset(react.tum, Patient %in% c("UPN001", "UPN003") & Frequency>10), 
                         ident.1 = "post-ACT",
                         ident.2 = "pre-ACT" 
) 


markers_R$gene <- rownames(markers_R)


#filter irrelevant genes
load("~/BaseTIL_code/data/Zheng_Japrin/exclude.gene.misc.human.v4.RData")
f.feat <- !(markers_R$gene %in% all.gene.ignore.df[["seu.id"]]) &
  !(grepl("^(RP[LS]|Rp[ls])",markers_R$gene,perl=T)) &
  !(grepl("^MT-|^MTRN|^MTAT|^MTND|^MRP",markers_R$gene,perl=T)) &
  !(grepl("^RP[0-9]+-",markers_R$gene,perl=T)) &
  !(markers_R$gene %in% c("MALAT1", 'XIST', 'NEAT1', "IGHV7-4-1"))
markers_R.flt <- markers_R[f.feat,]
print(head(markers_R.flt))
print(paste("Length of anchors:", nrow(markers_R)) )
print(paste("Length of filtered anchors:", nrow(markers_R.flt)) )


EnhancedVolcano::EnhancedVolcano(markers_R.flt, x = "avg_log2FC", y = "p_val_adj", lab = markers_R.flt$gene)






markers_NR <- FindMarkers(subset(react.tum, Patient %in% c("UPN009", "UPN008", "UPN011") & Frequency > 10), 
                          ident.1 = "post-ACT",
                          ident.2 = "pre-ACT" 
)


markers_NR$gene <- rownames(markers_NR)


#filter irrelevant genes
load("~/BaseTIL_code/data/Zheng_Japrin/exclude.gene.misc.human.v4.RData")
f.feat <- !(markers_NR$gene %in% all.gene.ignore.df[["seu.id"]]) &
  !(grepl("^(RP[LS]|Rp[ls])",markers_NR$gene,perl=T)) &
  !(grepl("^MT-|^MTRN|^MTAT|^MTND|^MRP",markers_NR$gene,perl=T)) &
  !(grepl("^RP[0-9]+-",markers_NR$gene,perl=T)) &
  !(markers_NR$gene %in% c("MALAT1", 'XIST', 'NEAT1', "IGHV7-4-1"))
markers_NR.flt <- markers_NR[f.feat,]
print(head(markers_NR.flt))
print(paste("Length of anchors:", nrow(markers_NR)) )
print(paste("Length of filtered anchors:", nrow(markers_NR.flt)) )


EnhancedVolcano::EnhancedVolcano(markers_NR.flt, x = "avg_log2FC", y = "p_val_adj", lab = markers_NR.flt$gene, parseLabels = T, drawConnectors = T)



markers_R.flt$abs_pct_diff <- abs(markers_R.flt$pct.1 - markers_R.flt$pct.2)
heatmap_markers_R <- markers_R.flt %>% arrange(p_val_adj, desc(abs_pct_diff)) %>% head(n=20)%>% arrange(avg_log2FC)

markers_NR.flt$abs_pct_diff <- abs(markers_NR.flt$pct.1 - markers_NR.flt$pct.2)
heatmap_markers_NR <- markers_NR.flt %>% arrange(p_val_adj, desc(abs_pct_diff)) %>% head(n=20) %>% arrange(avg_log2FC)


#Heatmap plotting
aggro <- AggregateExpression(subset(react.tum, Frequency>10), return.seurat = T, group.by = c("Patient", "prePost", "Response", "ntb"))
#aggro <- subset(react.tum, Frequency>5)


preR_genes <- c("VCAM1", "ID3", "MS4A6A", "NR4A3", "ITGAE", "NFKB1", "ZNF331", "IRF1")

postR_genes <- c("GZMA", "CXCR4")

preNR_genes <- c("KIR2DL4","KLRD1",  "KLRC4","KLRC2")

postNR_genes <- c("CXCL13", "SNX9","PDCD1", "CTLA4","HAVCR2","ADGRG1", "THEMIS" , "HLA-DRB5","GZMB", "GNLY","GZMH", "PRF1", "TNFSF4",  "NKG7") #, "CCL3", "CCL4", "CCL4L2", "CCL5"




all_markers <- c(preR_genes, postR_genes, preNR_genes, postNR_genes)


length(all_markers)
all_markers <- all_markers[all_markers %in% rownames(aggro@assays[["RNA"]]$data)]
length(all_markers)

mat<- aggro@assays[["RNA"]]$scale.data[all_markers, ] %>% as.matrix()
#mat<- t(scale(t(mat)))

quantile(mat, c(0.05, 0.95))


aggro@meta.data$Response2 <- ifelse(aggro@meta.data$Patient %in% c("UPN001", "UPN003"), "R", "NR")
aggro@meta.data$groups <- paste(aggro@meta.data$Response2, aggro@meta.data$prePost, sep = " ")
aggro@meta.data$groups <- factor(aggro@meta.data$groups, levels = c("R pre-ACT", "R post-ACT", "NR pre-ACT", "NR post-ACT"))



library(circlize)
library(ComplexHeatmap)

pdf("Fig/TR_prePost_Heatmap_genes.pdf", height = 6, width = 12)
Heatmap(mat, name = "Z score",  
        column_split = aggro@meta.data$groups,
        show_column_dend = T,
        cluster_columns = F,
        cluster_column_slices = TRUE,
        cluster_rows = F,
        column_title_gp = gpar(fontsize = 12),
        column_gap = unit(0.75, "mm"),
        show_row_dend = FALSE,
        col=colorRamp2(c(-1.5, 0, 1.5), c("#3C4DA0", "white", "#CD1F2D")),
        row_names_gp = gpar(fontsize = 12),
        column_title_rot = 0,
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("#A6CEE3", "#1F78B4", "#FB9A99", "#E31A1C")))),
        show_column_names = F,
        use_raster = TRUE,
        raster_quality = 4,
        row_names_side = "left",
        heatmap_legend_param = list(at = c(-1.5, -0.75, 0, 0.75, 1.5))
)
dev.off()


## -----------------------------------------------------------------------------------------------------------
#Signature scoring for selected clones of R and NR post-ACT

# looking for clones that underwend expansion
table(react.tum@meta.data[react.tum$Patient == "UPN001" & react.tum@meta.data$cloneType == "Hyperexpanded (200 < X <= 1000)", "ntb"])
table(react.tum@meta.data[react.tum$Patient == "UPN009" & react.tum@meta.data$cloneType == "Hyperexpanded (200 < X <= 1000)", "ntb"])


react.tum$highlight<- ifelse(react.tum$ntb == "TGTGCCAGCAGTGAAACTCTAGCGATAGCCGGGGAGCTGTTTTTT", "R", NA)

react.tum$highlight<- ifelse(react.tum$ntb == "TGTGCCAGCGGATTACCGGAGACCCAGTACTTC", "NR", react.tum$highlight)

# highlighting the clones in UAMP
pdf("Fig/Example_Clones_response.pdf", height = 3.5, width = 6)
DimPlot(react.tum, group.by = "highlight", split.by = "prePost", pt.size = 0.8, order = T) &
  ggtitle("") &
  scale_color_manual(na.value = "gray80", values = c("#E31A1C", "#1F78B4")) &
  ggeasy::easy_remove_axes() &
  theme(legend.position="bottom" )
dev.off()



# identifies differential expressed genes between the clones post-ACT
Idents(react.tum) <- react.tum$RNR
markers <- FindMarkers(subset(react.tum, prePost == "post-ACT"), ident.1 ="NRs", ident.2 ="Rs" )
markers$gene <- rownames(markers)



#filter irrelevant genes
load("~/BaseTIL_code/data/Zheng_Japrin/exclude.gene.misc.human.v4.RData")
f.feat <- !(markers$gene %in% all.gene.ignore.df[["seu.id"]]) &
  !(grepl("^(RP[LS]|Rp[ls])",markers$gene,perl=T)) &
  !(grepl("^MT-|^MTRN|^MTAT|^MTND|^MRP",markers$gene,perl=T)) &
  !(grepl("^RP[0-9]+-",markers$gene,perl=T)) &
  !(markers$gene %in% c("MALAT1", 'XIST', 'NEAT1'))
markers.flt <- markers[f.feat,]
print(head(markers.flt))
print(paste("Length of anchors:", nrow(markers)) )
print(paste("Length of filtered anchors:", nrow(markers.flt)) )

markers.flt$gene[markers.flt$p_val_adj < 0.01 & markers.flt$avg_log2FC > 0.5]

#select markers to plot
show.genes_pos <- c("HAVCR2", "PDCD1","CTLA4","LAG3",  "TIGIT","CXCL13",  "CXCR6", "KLRD1", "KIR2DL4",  "CCL3","CCL4", "CCL4L2")
show.genes_neg <- c("NR4A2", "NR4A3", "GZMK", "GZMM", "DUSP4","JUND", "FOSL2", "TNFSF9")
markers.flt$selection <- ifelse(markers.flt$gene %in% c(show.genes_pos, show.genes_neg),markers.flt$gene, NA)


sub_clones<- subset(react.tum, prePost == "post-ACT" & Type == "Rebiopsy1")
mat <- GetAssayData(sub_clones, layer = "scale.data")

mat<- mat[c(show.genes_pos,show.genes_neg),]
quantile(mat, c(0.05, 0.95))


#plot the genes in a heamap per cell
pdf("Fig/TR_example_clones_Heatmap_genes.pdf", height = 5, width = 10)
draw(Heatmap(mat, name = "Z score",
             column_split = factor(sub_clones@meta.data$RNR, levels = c("NRs", "Rs"), labels = c("NR", "R")),
             show_column_dend = T,
             cluster_columns = F,
             cluster_column_slices = TRUE,
             cluster_rows = F,
             column_title_gp = gpar(fontsize = 12),
             column_gap = unit(0.75, "mm"),
             show_row_dend = FALSE,
             col= colorRamp2(seq(from=-1.5,to=1.5,length=11),rev(brewer.pal(11, "RdBu"))),
             row_names_gp = gpar(fontsize = 12),
             column_title_rot = 0,
             top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("#E31A1C", "#1F78B4" )))),
             show_column_names = F,
             use_raster = TRUE,
             raster_quality = 4,
             row_names_side = "left",
             heatmap_legend_param = list(at = c(-1.5, -0.75, 0, 0.75, 1.5),
                                         legend_direction = "horizontal"
             )
             
             
), heatmap_legend_side = "bottom")
dev.off()


