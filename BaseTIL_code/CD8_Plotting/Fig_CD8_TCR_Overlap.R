## -----------------------------------------------------------------------------------------------------------
#load dependencies

suppressPackageStartupMessages({
  library(readr)
  library(Seurat)
  library(ggeasy)
  library(ggpubr)
  library(immunarch)
  library(viridis)
  library(purrr)
  library(tidyr)
  library(scRepertoire)
  library(stringr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
})

set.seed(1234567)
setwd("~/BaseTIL_code/CD8_Plotting")

gex.CD8 <- readRDS(file = "~/BaseTIL_code/saveFiles/CD8_annotated.rds")

## -----------------------------------------------------------------------------------------------------------
#calculates the TCR overlap between samples and then evaluates the mean over the patients


gex.CD8$Type_new <- factor(gex.CD8$Type, 
                           levels = c("Tumor","PreREP","Expanded_TILs", "PBMC_7dpt", "Rebiopsy1", "Rebiopsy2"),
                           labels= c("pre-ACT","Intermediate product","TIL product", "PBMC 7dpt", "post-ACT.1", "post-ACT.2"))
gex.CD8$Type <- as.character(gex.CD8$Type)

gex.CD8$Sample_Name_new <- paste(gex.CD8$Patient, gex.CD8$Type_new, sep = " ")

#calculate overlap
overlap_matrix <- clonalOverlap(subset(gex.CD8, Type != "PBMC_Ctrl"), 
                                cloneCall="nt", 
                                chain = "TRB", 
                                method = "overlap", group.by = "Sample_Name_new", exportTable = T)

overlap_matrix[lower.tri(overlap_matrix, diag = TRUE)] <- t(overlap_matrix)[lower.tri(overlap_matrix, diag = TRUE)]


#create empty matrix and fill it with values of mean
empty_matrix <- matrix(nrow = 5, ncol = 5)
Type_vector <- c("pre-ACT", "Intermediate product", "TIL product", "PBMC 7dpt", "post-ACT")
pat_vector <- c("UPN001", "UPN002", "UPN003", "UPN006", "UPN008", "UPN009", "UPN011")
colnames(empty_matrix) <- Type_vector
rownames(empty_matrix) <- Type_vector


for (first_type in c("pre-ACT", "Intermediate product", "TIL product", "PBMC 7dpt")) {
  for (second_type in c("pre-ACT", "Intermediate product", "TIL product", "PBMC 7dpt")) {
    empty_matrix[first_type, second_type] <- mean(diag(as.matrix(overlap_matrix[paste(pat_vector, first_type), paste(pat_vector, second_type)])))
  }
  
}


for (first_type in c("pre-ACT", "Intermediate product", "TIL product", "PBMC 7dpt")) {
  result_vector <- c()
  for (patient_vec in c("UPN001", "UPN006", "UPN008", "UPN009", "UPN011")) {
    post_vector <- list(UPN001=c("UPN001 post-ACT.1", "UPN001 post-ACT.2"),
                        UPN006=c("UPN006 post-ACT.1", "UPN006 post-ACT.2"),
                        UPN008=c("UPN008 post-ACT.2"),
                        UPN009=c("UPN009 post-ACT.1", "UPN009 post-ACT.2"),
                        UPN011=c("UPN011 post-ACT.1"))
    result_vector <- c(result_vector, as.vector(as.matrix(overlap_matrix[paste(patient_vec, first_type), post_vector[[patient_vec]]])))
  }
  empty_matrix[first_type, "post-ACT"] <- mean(result_vector)
}


for (second_type in c("pre-ACT", "Intermediate product", "TIL product", "PBMC 7dpt")) {
  result_vector <- c()
  for (patient_vec in c("UPN001", "UPN006", "UPN008", "UPN009", "UPN011")) {
    post_vector <- list(UPN001=c("UPN001 post-ACT.1", "UPN001 post-ACT.2"),
                        UPN006=c("UPN006 post-ACT.1", "UPN006 post-ACT.2"),
                        UPN008=c("UPN008 post-ACT.2"),
                        UPN009=c("UPN009 post-ACT.1", "UPN009 post-ACT.2"),
                        UPN011=c("UPN011 post-ACT.1"))
    result_vector <- c(result_vector, as.vector(as.matrix(overlap_matrix[post_vector[[patient_vec]], paste(patient_vec, second_type)])))
  }
  empty_matrix["post-ACT", second_type] <- mean(result_vector)
}

mat <- as.matrix(empty_matrix)

#plot the heatmap
pdf("Fig/CD8_Mean_overlap.pdf", height = 3.5, width = 5.5)
ht <- ComplexHeatmap::Heatmap(mat, 
                              name = "Mean Overlap",
                              column_order = paste(Type_vector),
                              row_order = paste(Type_vector),
                              cluster_columns = F,
                              cluster_rows = F,
                              column_names_rot = 45,
                              col = viridis(256), 
                              cell_fun = function(j, i, x, y, width, height, fill) {
                                if(mat[i, j] < 0.21 & !is.na(mat[i, j])){
                                  grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 10, col="white"))
                                } else {
                                  grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 10, col="black"))
                                }
                              }) 

draw(ht, padding = unit(c(2, 25, 2, 2), "mm"))
dev.off()





## -----------------------------------------------------------------------------------------------------------
#plot single TCR overlap values per patient

#patient UPN001
overlap_matrix <- clonalOverlap(subset(gex.CD8, Patient == "UPN001"), 
                                cloneCall="nt", 
                                chain = "TRB", 
                                method = "overlap", group.by = "Sample_Name", exportTable = T)


overlap_matrix[lower.tri(overlap_matrix, diag = TRUE)] <- t(overlap_matrix)[lower.tri(overlap_matrix, diag = TRUE)]
mat <- as.matrix(overlap_matrix)
Type_vector <- c("Tumor", "PreREP", "Expanded_TILs", "PBMC_7dpt", "Rebiopsy1", "Rebiopsy2")

pdf("Graphs/CD8_UPN001_overlap.pdf", height = 5, width = 7)
ht <- ComplexHeatmap::Heatmap(mat, 
                              name = "Overlap",
                              column_order = paste("UPN001", Type_vector),
                              row_order = paste("UPN001", Type_vector),
                              cluster_columns = F,
                              cluster_rows = F,
                              column_names_rot = 45,
                              heatmap_legend_param = list(at = seq(0, 0.5, 0.1)),
                              col = colorRamp2(seq(0, 0.5, 0.01), viridis(51)),
                              cell_fun = function(j, i, x, y, width, height, fill) {
                                grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 10))
                              })
draw(ht, padding = unit(c(2, 25, 2, 2), "mm"))
dev.off()

#patient UPN002
overlap_matrix <- clonalOverlap(subset(gex.CD8, Patient == "UPN002"), 
                                cloneCall="nt", 
                                chain = "TRB", 
                                method = "overlap", group.by = "Sample_Name", exportTable = T)


overlap_matrix[lower.tri(overlap_matrix, diag = TRUE)] <- t(overlap_matrix)[lower.tri(overlap_matrix, diag = TRUE)]
mat <- as.matrix(overlap_matrix)
Type_vector <- c("Tumor", "PreREP", "Expanded_TILs", "PBMC_7dpt")

pdf("Graphs/CD8_UPN002_overlap.pdf", height = 5, width = 7)
ht <- ComplexHeatmap::Heatmap(mat, 
                              name = "Overlap",
                              column_order = paste("UPN002", Type_vector),
                              row_order = paste("UPN002", Type_vector),
                              cluster_columns = F,
                              cluster_rows = F,
                              column_names_rot = 45,
                              heatmap_legend_param = list(at = seq(0, 0.5, 0.1)),
                              col = colorRamp2(seq(0, 0.5, 0.01), viridis(51)),
                              cell_fun = function(j, i, x, y, width, height, fill) {
                                grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 10))
                              })
draw(ht, padding = unit(c(2, 25, 2, 2), "mm"))
dev.off()


#patient UPN003
overlap_matrix <- clonalOverlap(subset(gex.CD8, Patient == "UPN003"), 
                                cloneCall="nt", 
                                chain = "TRB", 
                                method = "overlap", group.by = "Sample_Name", exportTable = T)

overlap_matrix[lower.tri(overlap_matrix, diag = TRUE)] <- t(overlap_matrix)[lower.tri(overlap_matrix, diag = TRUE)]
mat <- as.matrix(overlap_matrix)
Type_vector <- c("Tumor", "PreREP", "Expanded_TILs", "PBMC_7dpt")

pdf("Graphs/CD8_UPN003_overlap.pdf", height = 5, width = 7)
ht <- ComplexHeatmap::Heatmap(mat, 
                              name = "Overlap",
                              column_order = paste("UPN003", Type_vector),
                              row_order = paste("UPN003", Type_vector),
                              cluster_columns = F,
                              cluster_rows = F,
                              column_names_rot = 45,
                              heatmap_legend_param = list(at = seq(0, 0.5, 0.1)),
                              col = colorRamp2(seq(0, 0.5, 0.01), viridis(51)),
                              cell_fun = function(j, i, x, y, width, height, fill) {
                                grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 10))
                              })
draw(ht, padding = unit(c(2, 25, 2, 2), "mm"))
dev.off()


#patient UPN006
overlap_matrix <- clonalOverlap(subset(gex.CD8, Patient == "UPN006"), 
                                cloneCall="nt", 
                                chain = "TRB", 
                                method = "overlap", group.by = "Sample_Name", exportTable = T)

overlap_matrix[lower.tri(overlap_matrix, diag = TRUE)] <- t(overlap_matrix)[lower.tri(overlap_matrix, diag = TRUE)]
mat <- as.matrix(overlap_matrix)
Type_vector <- c("Tumor", "PreREP", "Expanded_TILs", "PBMC_7dpt", "Rebiopsy1", "Rebiopsy2")

pdf("Graphs/CD8_UPN006_overlap.pdf", height = 5, width = 7)
ht <- ComplexHeatmap::Heatmap(mat, 
                              name = "Overlap",
                              column_order = paste("UPN006", Type_vector),
                              row_order = paste("UPN006", Type_vector),
                              cluster_columns = F,
                              cluster_rows = F,
                              column_names_rot = 45,
                              heatmap_legend_param = list(at = seq(0, 0.5, 0.1)),
                              col = colorRamp2(seq(0, 0.5, 0.01), viridis(51)),
                              cell_fun = function(j, i, x, y, width, height, fill) {
                                grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 10))
                              })
draw(ht, padding = unit(c(2, 25, 2, 2), "mm"))
dev.off()


#patient UPN008
overlap_matrix <- clonalOverlap(subset(gex.CD8, Patient == "UPN008"), 
                                cloneCall="nt", 
                                chain = "TRB", 
                                method = "overlap", group.by = "Sample_Name", exportTable = T)

overlap_matrix[lower.tri(overlap_matrix, diag = TRUE)] <- t(overlap_matrix)[lower.tri(overlap_matrix, diag = TRUE)]
mat <- as.matrix(overlap_matrix)
Type_vector <- c("Tumor", "PreREP", "Expanded_TILs", "PBMC_7dpt", "Rebiopsy1", "Rebiopsy2")

pdf("Graphs/CD8_UPN008_overlap.pdf", height = 5, width = 7)
ht <- ComplexHeatmap::Heatmap(mat, 
                              name = "Overlap",
                              column_order = paste("UPN008", Type_vector),
                              row_order = paste("UPN008", Type_vector),
                              cluster_columns = F,
                              cluster_rows = F,
                              column_names_rot = 45,
                              heatmap_legend_param = list(at = seq(0, 0.5, 0.1)),
                              col = colorRamp2(seq(0, 0.5, 0.01), viridis(51)),
                              cell_fun = function(j, i, x, y, width, height, fill) {
                                grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 10))
                              })
draw(ht, padding = unit(c(2, 25, 2, 2), "mm"))
dev.off()

#patient UPN009
overlap_matrix <- clonalOverlap(subset(gex.CD8, Patient == "UPN009"), 
                                cloneCall="nt", 
                                chain = "TRB", 
                                method = "overlap", group.by = "Sample_Name", exportTable = T)

overlap_matrix[lower.tri(overlap_matrix, diag = TRUE)] <- t(overlap_matrix)[lower.tri(overlap_matrix, diag = TRUE)]
mat <- as.matrix(overlap_matrix)
Type_vector <- c("Tumor", "PreREP", "Expanded_TILs", "PBMC_7dpt", "Rebiopsy1", "Rebiopsy2")

pdf("Graphs/CD8_UPN009_overlap.pdf", height = 5, width = 7)
ht <- ComplexHeatmap::Heatmap(mat, 
                              name = "Overlap",
                              column_order = paste("UPN009", Type_vector),
                              row_order = paste("UPN009", Type_vector),
                              cluster_columns = F,
                              cluster_rows = F,
                              column_names_rot = 45,
                              heatmap_legend_param = list(at = seq(0, 0.5, 0.1)),
                              col = colorRamp2(seq(0, 0.5, 0.01), viridis(51)),
                              cell_fun = function(j, i, x, y, width, height, fill) {
                                grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 10))
                              })
draw(ht, padding = unit(c(2, 25, 2, 2), "mm"))
dev.off()


#patient UPN011
overlap_matrix <- clonalOverlap(subset(gex.CD8, Patient == "UPN011"), 
                                cloneCall="nt", 
                                chain = "TRB", 
                                method = "overlap", group.by = "Sample_Name", exportTable = T)

overlap_matrix[lower.tri(overlap_matrix, diag = TRUE)] <- t(overlap_matrix)[lower.tri(overlap_matrix, diag = TRUE)]
mat <- as.matrix(overlap_matrix)
Type_vector <- c("Tumor", "PreREP", "Expanded_TILs", "PBMC_7dpt", "Rebiopsy1")

pdf("Graphs/CD8_UPN011_overlap.pdf", height = 5, width = 7)
ht <- ComplexHeatmap::Heatmap(mat, 
                              name = "Overlap",
                              column_order = paste("UPN011", Type_vector),
                              row_order = paste("UPN011", Type_vector),
                              cluster_columns = F,
                              cluster_rows = F,
                              column_names_rot = 45,
                              heatmap_legend_param = list(at = seq(0, 0.5, 0.1)),
                              col = colorRamp2(seq(0, 0.5, 0.01), viridis(51)),
                              cell_fun = function(j, i, x, y, width, height, fill) {
                                grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 10))
                              })
draw(ht, padding = unit(c(2, 25, 2, 2), "mm"))
dev.off()




