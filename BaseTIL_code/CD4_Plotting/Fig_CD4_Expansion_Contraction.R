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
})

set.seed(1234567)
setwd("~/BaseTIL_code/CD4_Plotting")

gex.CD4 <- readRDS(file = "../saveFiles/CD4_annotated.rds")

#remove control PBMC sample from analysis
gex.CD4 <- subset(gex.CD4, Type != "PBMC_Ctrl")

#Remove sample UPN008 Rebiopsy1, non-tumor material was removed
gex.CD4 <- subset(gex.CD4, Sample_Name != "UPN008 Rebiopsy1")



## -----------------------------------------------------------------------------------------------------------
# Calculate expansion and contraction dynamics based on logFC between consecutive samples

meta <- gex.CD4@meta.data[!(is.na(gex.CD4$ntb) | gex.CD4$ntb == "NA"),]
meta$pat_barcode <- paste(meta$Patient , meta$ntb, sep = "-")



res_data <- meta %>%
  group_by(Sample_Name, ntb) %>%
  summarize(clone_count = n()) %>%
  group_by(Sample_Name) %>%
  mutate(clone_percent = clone_count/sum(clone_count)*100)%>% 
  tidyr::separate(Sample_Name, into = c("Patient", "Type"), sep=" ", remove = F) %>% 
  unite(col = "pat_barcode", Patient, ntb, sep = "-") %>% 
  group_by(pat_barcode)%>% dplyr::select(-Sample_Name,)  %>% 
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
rownames(merged) <- merged$barcode
merged <- merged %>% dplyr::select(Tumor_to_PreREP,Tumor_to_Exp, PreREP_to_Exp, Exp_to_PBMC, PBMC_to_R1, PBMC_to_R2,PBMC_to_Post, Tumor_to_R1, Tumor_to_R2,Tumor_to_Post, Exp_to_R1, Exp_to_R2)


merged.gex_noNA <- AddMetaData(gex.CD4, metadata = merged)
merged.gex_noNA<- subset(merged.gex_noNA, cells = rownames(gex.CD4@meta.data[!is.na(gex.CD4$ntb),]))


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
# Feature plot of tumor reactive and virus reactive signatures in pre-ACT

markers <- list()


NeoTCR_sig <- read.csv(file = "../data/Lowery_NeoTCR.csv")
markers$NeoTCR4 <- NeoTCR_sig$NeoTCR4
markers$T_EX <- c("TOX", "LAG3", "CASP8", "TIGIT", "CTLA4", "ITGAE", "PDCD1", "CXCL13","CD27")
markers$Tfh <- c( "CD40LG" ,"TOX2"   ,"MAF"   , "CD200",  "BATF"  )
markers$HanadaCD4 <- c("CXCL13", "NR3C1", "ADGRG1", "NMG", "ITM2A", "ETV7", "COTL1", "B2M", "IGFL2", "CD4")


tum <- AddModuleScore_UCell(tum, features = markers, name = "")


pdf("Fig/CD4_NeoTCR_sig_featureplot.pdf",  height = 3, width = 4)
FeaturePlot(tum, features = c("HanadaCD4"), min.cutoff = 0.4, ncol = 1, order = T) +
    ggeasy::easy_remove_axes()
FeaturePlot(tum, features = c("NeoTCR4"), min.cutoff = 0.25, ncol = 1, order = T) + ggeasy::easy_remove_axes()
Nebulosa::plot_density(tum, c("HanadaCD4")) & 
  ggeasy::easy_remove_axes() & 
  ggeasy::easy_remove_legend_title() &
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        plot.title = element_text(hjust = 0.5))
Nebulosa::plot_density(tum, c("NeoTCR4")) & 
  ggeasy::easy_remove_axes() & 
  ggeasy::easy_remove_legend_title() &
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        plot.title = element_text(hjust = 0.5))
dev.off()



## -----------------------------------------------------------------------------------------------------------
# Highlighting of overall reactive clones in UMAP

tum$reactive_temp <- ifelse(tum$overall_reactive == TRUE, tum$Patient, NA)
tum$reactive_temp <- factor(tum$reactive_temp, levels = rev(c("UPN001", "UPN002", "UPN003", "UPN006", "UPN011", "UPN008", "UPN009"))) #


pdf("Fig/CD4_DimPlot_Reactive_Overall_baseTumor.pdf", height = 3, width = 4.5)
DimPlot(subset(tum, Patient != "UPN011"), 
        group.by = "reactive_temp", raster = F,  order = T, pt.size = 0.5
)& easy_remove_axes()  & 
  scale_color_manual(name = "Patients", values = c(
    "UPN001" =  "#BC3C29B2",
    "UPN002" = "#0072B5FF",
    "UPN003" = "#E18727B2",
    "UPN006" = "#20854EB2",
    "UPN008" = "#7876B1B2",
    "UPN009" = "#6F99ADB2",
    "UPN011" = "#EE4C97B2"
  ), na.value = "gray80")  & 
  ggtitle("Tumor-reactive CD4 T cells") & theme(strip.text.x = element_text(size=20, face = "bold"),
                        legend.text = element_text(size=20),
                        legend.title = element_text(size=20, face = "bold"))
dev.off()

