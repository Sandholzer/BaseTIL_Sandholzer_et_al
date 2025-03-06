## -----------------------------------------------------------------------------------------------------------
#load dependencies

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(stringr)
  library(ggplot2)
  library(ggalluvial)
  library(ggpubr)
  library(readr)
  library(ggeasy)
  library(purrr)
  library(ComplexHeatmap)
  library(UCell)
})

set.seed(1234567)
setwd("~/BaseTIL_code/CD8_Plotting")

gex.CD8 <- readRDS(file = "../saveFiles/CD8_annotated.rds")

#remove control PBMC sample from analysis
gex.CD8 <- subset(gex.CD8, Type != "PBMC_Ctrl")

#Remove sample UPN008 Rebiopsy1, non-tumor material was resected
gex.CD8 <- subset(gex.CD8, Sample_Name != "UPN008 Rebiopsy1")



## -----------------------------------------------------------------------------------------------------------
# Calculate expansion and contraction dynamics based on logFC between consecutive samples

gex.CD8@meta.data$cellnames <- Cells(gex.CD8)
meta <- gex.CD8@meta.data[!(is.na(gex.CD8$ntb) | gex.CD8$ntb == "NA"),]
meta$pat_barcode <- paste(meta$Patient , meta$ntb, sep = "-")


res_data <- meta %>%
  group_by(Sample_Name, ntb) %>%
  summarize(clone_count = n()) %>%
  group_by(Sample_Name) %>%
  mutate(clone_percent = clone_count/sum(clone_count)*100)%>% 
  tidyr::separate(Sample_Name, into = c("Patient", "Type"), sep=" ", remove = F) %>% 
  unite(col = "pat_barcode", Patient, ntb, sep = "-") %>% 
  group_by(pat_barcode)%>% select(-Sample_Name,)  %>% 
  pivot_wider(names_from = Type, values_from = c(clone_percent, clone_count))


res_data <- res_data %>% mutate(Tumor_to_PreREP = log(clone_percent_PreREP/clone_percent_Tumor), 
                                Tumor_to_Exp = log(clone_percent_Expanded_TILs/clone_percent_Tumor),
                                PreREP_to_Exp = log(clone_percent_Expanded_TILs/clone_percent_PreREP),
                                Exp_to_PBMC = log(clone_percent_PBMC_7dpt/clone_percent_Expanded_TILs),
                                PBMC_to_R1 = log(clone_percent_Rebiopsy1/clone_percent_PBMC_7dpt),
                                PBMC_to_R2 = log(clone_percent_Rebiopsy2/clone_percent_PBMC_7dpt),
                                PBMC_to_Post = log(mean(c(clone_percent_Rebiopsy1,clone_percent_Rebiopsy2), na.rm = T) /clone_percent_PBMC_7dpt),
                                Tumor_to_R1 = log(clone_percent_Rebiopsy1/clone_percent_Tumor),
                                Tumor_to_R2 = log(clone_percent_Rebiopsy2/clone_percent_Tumor),
                                Tumor_to_Post = log(mean(c(clone_percent_Rebiopsy1,clone_percent_Rebiopsy2), na.rm = T)/clone_percent_Tumor),
                                Exp_to_R1 = log(clone_percent_Rebiopsy1/clone_percent_Expanded_TILs),
                                Exp_to_R2 = log(clone_percent_Rebiopsy2/clone_percent_Expanded_TILs)
                                ) %>% as.data.frame()


merged <- left_join(meta, res_data[,c("pat_barcode","Tumor_to_PreREP", "Tumor_to_Exp", "PreREP_to_Exp", "Exp_to_PBMC", "PBMC_to_R1", "PBMC_to_R2", "PBMC_to_Post","Tumor_to_R1", "Tumor_to_R2", "Tumor_to_Post", "Exp_to_R1", "Exp_to_R2")], by="pat_barcode")
rownames(merged) <- merged$cellnames
merged <- merged %>% select(Tumor_to_PreREP,Tumor_to_Exp, PreREP_to_Exp, Exp_to_PBMC, PBMC_to_R1, PBMC_to_R2,PBMC_to_Post, Tumor_to_R1, Tumor_to_R2,Tumor_to_Post, Exp_to_R1, Exp_to_R2)


gex.CD8 <- AddMetaData(gex.CD8, metadata = merged)
merged.gex_noNA<- subset(gex.CD8, cells = rownames(gex.CD8@meta.data[!is.na(gex.CD8$ntb),]))


reactive <- subset(merged.gex_noNA, overall_reactive == TRUE)

## Save object of tumor-reactive CD8 T cells for further analysis
saveRDS(reactive, "../saveFiles/CD8_tumor_reactive.rds")




## -----------------------------------------------------------------------------------------------------------
# Plotting of expansion and contraction dynamics in UMAP and bar plot

# Subset baseline tumor samples to highlight dynamics in the UMAP of pre-ACT
tum <- subset(merged.gex_noNA, subset = Type %in% c("Tumor"))


# Plot UMAPS with highlighted expansion and contraction values per cell

plt <- list()


# Reorder the cells to position NA cells in the back and randomize all others to prevent sample bias
ordered_indices <- c(which(is.na(tum@meta.data$Tumor_to_PreREP)), sample(which(!is.na(tum@meta.data$Tumor_to_PreREP))))
new.cell.order <- rownames(tum@meta.data[ordered_indices,])

plt[["Expansion_T_to_Pre_inTumor"]]<- FeaturePlot(
  tum,
  features = "Tumor_to_PreREP",
  cells =new.cell.order, 
  #pt.size = 0.5,
  order = F,
  min.cutoff = -1.6,
  max.cutoff = 1.6
) & ggtitle("Expansion in preREP") & 
  theme(
  plot.title = element_text(size=14, face = "bold", hjust = 0.5),
  legend.title = element_text(size = 12)) &
  scale_color_gradient2(name = "Expansion \nScore",
                        low = "blue",
                        high = "red",
                        mid = "white",
                        na.value = "gray80",
                        guide = "colourbar"
  ) & ggeasy::easy_remove_axes()


# Reorder the cells to position NA cells in the back and randomize all others to prevent sample bias
ordered_indices <- c(which(is.na(tum@meta.data$PreREP_to_Exp)), sample(which(!is.na(tum@meta.data$PreREP_to_Exp))))
new.cell.order <- rownames(tum@meta.data[ordered_indices,])

plt[["Expansion_Pre_to_Exp_inTumor"]]<- FeaturePlot(
  tum,
  features = "PreREP_to_Exp",
  cells =new.cell.order, 
  #pt.size = 0.5,
  order = F,
  min.cutoff = -1.6,
  max.cutoff = 1.6
) & ggtitle("Expansion in REP")& 
  theme(
    plot.title = element_text(size=14, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 12)) &
  scale_color_gradient2(name = "Expansion \nScore",
                        low = "blue",
                        high = "red",
                        mid = "white",
                        na.value = "gray80",
                        guide = "colourbar"
  ) & ggeasy::easy_remove_axes()





# Reorder the cells to position NA cells in the back and randomize all others to prevent sample bias
ordered_indices <- c(which(is.na(tum@meta.data$Exp_to_PBMC)), sample(which(!is.na(tum@meta.data$Exp_to_PBMC))))
new.cell.order <- rownames(tum@meta.data[ordered_indices,])


plt[["Expansion_Exp_to_PBMC_inTumor"]]<- FeaturePlot(
  tum,
  features = "Exp_to_PBMC",
  cells =new.cell.order, 
  #pt.size = 0.5,
  order = F,
  min.cutoff = -1.6,
  max.cutoff = 1.6
) & ggtitle("Expansion in blood")& 
  theme(
    plot.title = element_text(size=14, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 12)) &
  scale_color_gradient2(name = "Expansion \nScore",
                        low = "blue",
                        high = "red",
                        mid = "white",
                        na.value = "gray80",
                        guide = "colourbar"
  )& ggeasy::easy_remove_axes()




# Subset post treatment samples to highlight dynamics in the UMAP of post-ACT
post <- subset(merged.gex_noNA, subset = Type %in% c("Rebiopsy1", "Rebiopsy2"))


# Reorder the cells to position NA cells in the back and randomize all others to prevent sample bias
ordered_indices <- c(which(is.na(post@meta.data$PBMC_to_Post)), sample(which(!is.na(post@meta.data$PBMC_to_Post))))
new.cell.order <- rownames(post@meta.data[ordered_indices,])


plt[["Expansion_PBMC_to_post_inPost"]]<- FeaturePlot(
  post,
  #split.by = "overall_reactive",
  features = "PBMC_to_Post",
  cells =new.cell.order, 
  #pt.size = 0.5,
  order = F,
  min.cutoff = -1.6,
  max.cutoff = 1.6
) & ggtitle("Expansion post therapy")& 
  theme(
    plot.title = element_text(size=14, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 12)) &
  scale_color_gradient2(name = "Expansion \nScore",
                        low = "blue",
                        high = "red",
                        mid = "white",
                        na.value = "gray80",
                        guide = "colourbar"
  )& ggeasy::easy_remove_axes()



# Plotting bar plots for quantification of cells per comparison within pre-ACT time point


for (comparison in c("Tumor_to_PreREP", "PreREP_to_Exp", "Exp_to_PBMC")) {
  
  
  tum@meta.data[,"Expansion_Score"] <- "undetected"
  tum@meta.data[!is.na(tum@meta.data[,comparison]) & tum@meta.data[,comparison] >= -0.25 & tum@meta.data[,comparison] < 0.25, "Expansion_Score"] <- "-0.25 \u2264 x < 0.25"
  tum@meta.data[!is.na(tum@meta.data[,comparison]) & tum@meta.data[,comparison] < -0.25, "Expansion_Score"] <- "x < -0.25"
  tum@meta.data[!is.na(tum@meta.data[,comparison]) & tum@meta.data[,comparison] >= 0.25, "Expansion_Score"] <- "x \u2265 0.25"
  
  tum@meta.data[,"Expansion_Score"]  <- factor(tum@meta.data[,"Expansion_Score"] , levels = rev(c( "undetected","x < -0.25", "-0.25 \u2264 x < 0.25", "x \u2265 0.25")))
  
  
  plt[[comparison]]<- dittoSeq::dittoBarPlot(subset(tum, Expansion_Score != "undetected"), var = "Expansion_Score", group.by = "Celltype",retain.factor.levels = T, scale = "count" 
  ) + xlab("")+ggtitle(NULL)+
    #ggtitle(comparison) + 
    theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size = 12)
    ) +
    scale_fill_manual(name = "Expansion \nScore",
                      values = c(
                        "x \u2265 0.25" = "#FF3219",
                        "-0.25 \u2264 x < 0.25" = "gray70",
                        "x < -0.25" = "#4423FF"
                      ), 
                      labels = c("x \u2265 0.25" ,
                                 "-0.25 \u2264 x < 0.25",
                                 "x < -0.25"))
  
  
}



# Plotting bar plots for quantification of cells per comparison within post-ACT time point


for (comparison in c("PBMC_to_Post")) {
  
  
  post@meta.data[,"Expansion_Score"] <- "undetected"
  post@meta.data[!is.na(post@meta.data[,comparison]) & post@meta.data[,comparison] >= -0.25 & post@meta.data[,comparison] < 0.25, "Expansion_Score"] <- "-0.25 \u2264 x < 0.25"
  post@meta.data[!is.na(post@meta.data[,comparison]) & post@meta.data[,comparison] < -0.25, "Expansion_Score"] <- "x < -0.25"
  post@meta.data[!is.na(post@meta.data[,comparison]) & post@meta.data[,comparison] >= 0.25, "Expansion_Score"] <- "x \u2265 0.25"
  
  post@meta.data[,"Expansion_Score"]  <- factor(post@meta.data[,"Expansion_Score"] , levels = rev(c( "undetected","x < -0.25", "-0.25 \u2264 x < 0.25", "x \u2265 0.25")))
  
  plt[[comparison]]<- dittoSeq::dittoBarPlot(subset(post, Expansion_Score != "undetected"), var = "Expansion_Score", group.by = "Celltype",retain.factor.levels = T, scale = "count"
  ) + xlab("")+ ggtitle(NULL) +
    theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size = 12)
          ) +
    #ggtitle(comparison) + 
    scale_fill_manual(name = "Expansion \nScore",
                      values = c(
                        "x \u2265 0.25" = "#FF3219",
                        "-0.25 \u2264 x < 0.25" = "gray70",
                        "x < -0.25" = "#4423FF"
                      ), 
                      labels = c("x \u2265 0.25" ,
                                 "-0.25 \u2264 x < 0.25",
                                 "x < -0.25"))
  
}

p1<- ggarrange(plotlist = plt[1:4], nrow = 1, ncol = 4, common.legend = T, legend = "right" )+
  theme(plot.margin = margin(0,1,0.0,0.5, "cm")) 
p2<- ggarrange(plotlist = plt[5:8], nrow = 1, ncol = 4, common.legend = T, legend = "right")


grDevices::cairo_pdf("Fig/Expansion_together.pdf", height = 6, width = 14)
ggarrange(p1, p2, ncol = 1)
dev.off()





## -----------------------------------------------------------------------------------------------------------
# Feature plot of tumor reactive and virus reactive signatures in lesion derived samples (pre-ACT and post-ACT)

prepostTumor <- subset(merged.gex_noNA, subset = Type %in% c("Tumor", "Rebiopsy1", "Rebiopsy2"))

markers <- list()

#Using own defined signature of tumor reactive T cells and filter it
own_signature_preflt  <- read.csv("Fig/DGE_results/dge_tumorReactive_signature_merged.csv")
own_signature_preflt <- own_signature_preflt[own_signature_preflt$avg_log2FC.sc > 0.5 & own_signature_preflt$avg_log2FC.bulk > 0.5,]
own_signature_preflt$name <- ifelse(own_signature_preflt$p_val_adj.sc < 1e-200 & own_signature_preflt$p_val.bulk < 0.05 ,own_signature_preflt$gene, NA )
reactive_signature<- own_signature_preflt$name
reactive_signature <- reactive_signature[!is.na(reactive_signature)]
write_csv(own_signature_preflt, "Fig/DGE_results/dge_tumorReactive_signature_merged_filtered.csv")

markers$own_signature <- reactive_signature

NeoTCR_sig <- read.csv(file = "../data/Lowery_NeoTCR.csv")
markers$NeoTCR8 <- NeoTCR_sig$NeoTCR8
markers$HanadaCD8 <- c("PDCD1", "ITGAE", "CXCL13", "ENTPD1", "BATF", "GZMB", "CD27", "TIGIT", "PHLDA1", "CD74", "HLA-DMA", "HLA-DRA", "HLA-DRB1", "HLA-DPB1", "CD3D", "CD82", "ARL3", "HMOX1", "ALOX5AP", "DUSP4", "CARS", "LSP1", "CCND2", "TPI1", "GAPDH", "ITM2A", "HMGN3", "CHST12", "NAP1L4")
published_sig <- read.csv(file = "../data/Oliveira_Tumor-spec_Virus-spec.csv")
markers$Oliveira_Tumor_spec <- published_sig$Oliveira.Tumor.specific
markers$Oliveira_Virs_spec <- published_sig$Oliveira.Virus.specific

prepostTumor <- AddModuleScore_UCell(prepostTumor, features = markers, name = "")


FeaturePlot(prepostTumor, features = names(markers), min.cutoff = 0.25,  order = T) +
  ggeasy::easy_remove_axes()


pdf("Fig/CD8_NeoTCR_sig_featureplot.pdf",  height = 3, width = 4)
FeaturePlot(prepostTumor, features = "own_signature", min.cutoff = 0.15,  order = T) +
    ggtitle("Tumor-reactive signature")+
  ggeasy::easy_remove_axes()
FeaturePlot(prepostTumor, features = "NeoTCR8", min.cutoff = 0.25,  order = T) +
  ggeasy::easy_remove_axes()

FeaturePlot(prepostTumor, features = "Oliveira_Tumor_spec", min.cutoff = 0.15,  order = T) +
  ggeasy::easy_remove_axes()

FeaturePlot(prepostTumor, features = "Oliveira_Virs_spec", min.cutoff = 0.2,  order = T) +
  ggeasy::easy_remove_axes()

Nebulosa::plot_density(prepostTumor, c("NeoTCR8")) & 
  ggeasy::easy_remove_axes() & 
  ggeasy::easy_remove_legend_title() &
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        plot.title = element_text(hjust = 0.5))

Nebulosa::plot_density(prepostTumor, c("Oliveira_Tumor_spec")) & 
  ggeasy::easy_remove_axes() & 
  ggeasy::easy_remove_legend_title() &
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        plot.title = element_text(hjust = 0.5))

Nebulosa::plot_density(prepostTumor, c("Oliveira_Virs_spec")) & 
  ggeasy::easy_remove_axes() & 
  ggeasy::easy_remove_legend_title() &
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        plot.title = element_text(hjust = 0.5))
dev.off()






