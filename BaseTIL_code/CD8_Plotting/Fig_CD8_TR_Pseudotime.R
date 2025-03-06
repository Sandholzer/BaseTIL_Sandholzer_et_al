## -----------------------------------------------------------------------------------------------------------
#load dependencies

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratWrappers)
  library(monocle3)
  library(ggplot2)
  library(patchwork)
  library(ComplexHeatmap)
  library(dplyr)
  library(RColorBrewer)
  library(circlize)
})

set.seed(1234567)
setwd("~/BaseTIL_code/CD8_Plotting/")
react.tum <-  readRDS(file = "~/BaseTIL_code/saveFiles/CD8_tumor_reactive_Tumor_ann.rds")


## -----------------------------------------------------------------------------------------------------------


#create cell data set for monocle3
cds <- as.cell_data_set(react.tum, assay = "RNA") 
cds <- estimate_size_factors(cds)

#combine all cells in one partition for just one pseudotime trajectory
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
recreate.partitions

cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

#transfer seurat clustering
list.cluster <- react.tum$RNA_snn_res.0.4
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster

cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- react.tum@reductions$umap@cell.embeddings



#learn cell trajectories
cds <- learn_graph(cds, use_partition = T, close_loop = T, learn_graph_control = NULL)


plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F,
           label_branch_points = T, label_roots = T, label_leaves = F,
           group_label_size = 5)

#order cells and use cluster 5 (Tn) as starting point 
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds[, clusters(cds) %in% c(5)]))


plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = F,
           label_branch_points = F, label_roots = F, label_leaves = F, show_trajectory_graph = T, cell_size = 0.7)

## -----------------------------------------------------------------------------------------------------------
#plot pseudotime highlighted on umap
pdf("Fig/TR_Pseudotime.pdf", height = 4, width = 6)
plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = F,
           label_branch_points = F, label_roots = F, label_leaves = F, show_trajectory_graph = F, cell_size = 0.7) + 
            ggeasy::easy_remove_axes() +
            theme(axis.line.y.left = element_blank(),
                  axis.line.x.bottom = element_blank())
dev.off()


#transfer pseudotime to seurat object
react.tum$pseudotime <- pseudotime(cds)
FeaturePlot(object = react.tum, feature= "pseudotime") + viridis::scale_color_viridis(option = "plasma")




## -----------------------------------------------------------------------------------------------------------
#plot boxplot to show pseudotime per cluster

library(RColorBrewer)
Colors <- RColorBrewer::brewer.pal(n=12, "Paired")


head(pseudotime(cds), 10)

cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

df("Fig/Boxplot_Pseudotime.pdf", height = 4, width = 4)
ggplot(data.pseudo, aes(reorder(TR_anno, monocle3_pseudotime), monocle3_pseudotime, fill = TR_anno)) + 
  geom_boxplot() + theme_classic() + scale_fill_manual(values=Colors) + 
  xlab("")+
  ylab("pseudotime")+ theme(axis.text.x = element_text(angle = 45, hjust=1)) + NoLegend()
dev.off()

library(dplyr)


## -----------------------------------------------------------------------------------------------------------
#DGE on trajectory path


#test for genes changing along the trajectory path
modulated_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 10)
genes <- row.names(subset(modulated_genes, q_value == 0 & morans_I > 0.2)) 

#remove genes for plotting
genes <- genes[!(genes %in% c(all.gene.ignore.df[["seu.id"]], "FOXP3", "CCR8", "CD8A", "B2M", "ACTB", "GAPDH", "TNFRSF18"))]

#include additional genes for plotting
genes <- c(genes, "GZMK","IFNG", "CXCL13", "PDCD1", "CCL4", "CCL4L2", "CCL3", "CXCR6", "NR4A2")


#create matrix and order along pseudotime
pt.matrix <- exprs(cds)[match(genes,rownames(rowData(cds))),order(pseudotime(cds))]

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes;


#add heatmap annotation for pseudotime scale
ha <-  HeatmapAnnotation(
  pseudotime = seq(
    length.out = 9062,
    from = 0,
    to = 8
  ),
  col = list(
    pseudotime = colorRamp2(seq(
      length.out = 100,
      from = 0,
      to = 8
    ), 
    plasma(100))),
  annotation_legend_param = list(pseudotime = list(direction = "horizontal"
  )))


#create heatmap using K means Clustering
htkm <- Heatmap(
  pt.matrix,
  name                         = "Z score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "RdYlBu"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  #clustering_method_rows = "complete",
  km = 4,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE, 
  top_annotation = ha,
   
  
  heatmap_legend_param = list(legend_direction = "horizontal" 
                              #legend_width = unit(3, "cm")
  )
  
  )


#create heatmap using Hierarchical Clustering
hthc <- Heatmap(
  pt.matrix,
  name                         = "Z score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "RdYlBu"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation = ha,
  heatmap_legend_param = list(legend_direction = "horizontal" 
                              #legend_width = unit(3, "cm")
  ))


#plot heatmaps
pdf("Graphs/TR_trajectory_heatmap.pdf", height = 8, width = 6)
draw(htkm, heatmap_legend_side="top", annotation_legend_side="top", merge_legend = TRUE)
draw(hthc, heatmap_legend_side="top", annotation_legend_side="top", merge_legend = TRUE)
dev.off()



