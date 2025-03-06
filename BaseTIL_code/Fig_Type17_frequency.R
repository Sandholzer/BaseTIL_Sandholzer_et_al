## -----------------------------------------------------------------------------------------------------------
#load dependencies

suppressPackageStartupMessages({
library(Seurat)
library(stringr)
library(ggplot2)
library(dittoSeq)
library(ggeasy)
library(ggpubr)
library(purrr)
library(tidyr)
library(dplyr)
})
setwd("~/BaseTIL_code")

## -----------------------------------------------------------------------------------------------------------
#load both datasets

gex.CD8 <- readRDS(file = "saveFiles/CD8_annotated.rds")
gex.CD8 <- subset(gex.CD8, Type != "PBMC_Ctrl" & Patient != "UPN002")
gex.CD8$Type_temp <- factor(gex.CD8$prepost,
                            levels = rev(c("pre-ACT","Intermediate product","TIL product", "PBMC 7dpt", "post-ACT")))
gex.CD8$RNR <- ifelse(gex.CD8$Response == "PR", "Rs", "NRs")


gex.CD4 <- readRDS(file = "saveFiles/CD4_annotated.rds")
gex.CD4 <- subset(gex.CD4, Type != "PBMC_Ctrl" & Patient != "UPN002")

gex.CD4$Type_temp <- factor(gex.CD4$prepost,
                            levels = rev(c("pre-ACT","Intermediate product","TIL product", "PBMC 7dpt", "post-ACT")))
gex.CD4$RNR <- ifelse(gex.CD4$Response == "PR", "Rs", "NRs")


## -----------------------------------------------------------------------------------------------------------
# Calculate Type17 percent

#Calculate Tc17 
gex.CD8@meta.data$isTc17 <- FALSE
gex.CD8@meta.data[gex.CD8@meta.data$Celltype == "Tc17","isTc17"] <- TRUE
data_Tc17<- dittoSeq::dittoBarPlot(gex.CD8, var = "isTc17", group.by = "Sample_Name_new", split.by = "RNR", data.out = T)
result_Tc17<- data_Tc17$data %>% filter(label==TRUE) %>% mutate(percent = percent*100)

#Calculate Th17 
gex.CD4@meta.data$isTh17 <- FALSE
gex.CD4@meta.data[gex.CD4@meta.data$Celltype == "Th17","isTh17"] <- TRUE
data_Th17<- dittoSeq::dittoBarPlot(gex.CD4, var = "isTh17", group.by = "Sample_Name_new", split.by = "RNR", data.out = T)
result_Th17<- data_Th17$data %>% filter(label==TRUE) %>% mutate(percent = percent*100)

#Calculate total Type 17
result <- merge(result_Tc17, result_Th17, by = c("label", "grouping", "RNR"), suffixes = c(".CD8", ".CD4"))
result$percent <- (result$count.CD8 + result$count.CD4) / (result$label.count.total.per.facet.CD8 + result$label.count.total.per.facet.CD4) * 100
result[,c("Patient", "Type")]<- str_split_fixed(result$grouping, pattern = " ", n = 2)
result <- subset(result, Type %in% c("pre-ACT",
                                     "Intermediate product",
                                     "TIL product"))

#Plot frequencies of Type17 during expansion
pdf("Fig/Type17_freq_expansion.pdf", height = 4, width = 5)
result %>% ggplot(aes(
  factor(
    Type,
    levels = c(
      "pre-ACT",
      "Intermediate product",
      "TIL product"
    )
  ),
  percent,
  color = RNR,
  group = Patient
)) + ggtitle("Type 17 T cells") +
  geom_point() + geom_line() +
  theme_classic() +
  theme( 
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = 16),
    axis.line = element_line(color = "black"),
    plot.title = element_text(size = 16,hjust = 0.5, face = "bold")) +
  labs(
    x = "",
    y = "Percent of T cells"
  )+ 
  scale_color_manual(name = "Response", values = c("NRs" = "#CD2626", "Rs" = "#009ACD")) +
  theme(axis.text.x = element_text(angle = 45,  hjust=1)) 
dev.off()


#Plot frequencies of Type17 in pre-ACT lesions
pdf("Fig/Type17_freq_Tumor.pdf", height = 3, width = 2)
bar_colors <- c("NRs" = "#CD2626", "Rs" = "#009ACD")

ggplot(subset(result, Type == "pre-ACT"), aes(x = RNR, y = percent, fill = RNR)) +
  geom_bar(stat = "summary", fun = "mean", alpha = 0.4, width = 0.7, position = position_dodge(width = 0.8), color = NA, show.legend = FALSE) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, aes(color = RNR), position = position_dodge(width = 0.8), show.legend = FALSE) +
  geom_jitter(aes(color = RNR), position = position_jitterdodge(jitter.width = 0.8, dodge.width = 0.8, seed = 3), size = 1.5, show.legend = FALSE) +
  scale_fill_manual(name = "Response", values = bar_colors) +
  scale_color_manual(values = bar_colors) +
  theme_minimal() +
  labs(x = "", y = "Percent of T cells") +
  ggtitle("Type 17")+
  theme_classic()+
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = 16),
    axis.line = element_line(color = "black"),
    plot.title = element_text(size = 15,hjust = 0.5, face = "bold")
  ) + 
  ylim(c(0,14.5))
dev.off()







## -----------------------------------------------------------------------------------------------------------
# Highlighting Tc17 clones at different time points


Exp.CD8 <-  readRDS(file = "../saveFiles/CD8_Expanded_TILs_subclustering.rds")

gex.CD8@meta.data$Patient_ntb <- paste(gex.CD8@meta.data$Patient, gex.CD8@meta.data$ntb, sep = "_")

Exp.CD8@meta.data$Patient_ntb <- paste(Exp.CD8@meta.data$Patient, Exp.CD8@meta.data$ntb, sep = "_")
Tc17_ntb<- Exp.CD8@meta.data %>% filter(Celltype_TIL == "Tc17" & Type == "Expanded_TILs" & ntb != "NA" & !is.na(ntb) & nta != "NA" & !is.na(nta)) %>% pull(Patient_ntb)
non_Tc17_ntb<- Exp.CD8@meta.data %>% filter(Celltype_TIL != "Tc17" & Type == "Expanded_TILs" & ntb != "NA" & !is.na(ntb) & nta != "NA" & !is.na(nta)) %>% pull(Patient_ntb)

Tc17_ntb<- Tc17_ntb[!(Tc17_ntb %in% non_Tc17_ntb)]

gex.CD8$Type_new_temp <- ifelse(gex.CD8$Type_new %in% c("post-ACT.1", "post-ACT.2"), "post-ACT", as.character(gex.CD8$Type_new))
gex.CD8$Type_new_temp <- factor(gex.CD8$Type_new_temp, levels = c("pre-ACT", "Intermediate product", "TIL product", "PBMC 7dpt", "post-ACT"))


gex.CD8$is_Th17_ntb <- ifelse(gex.CD8$Patient_ntb %in% Tc17_ntb, "Tc17", NA)

pdf("Fig/CD8_Tc17_clones.pdf", height = 4, width = 13)
DimPlot(subset(gex.CD8, subset = Type_new_temp %in% c("pre-ACT", "Intermediate product", "TIL product")), group.by="is_Th17_ntb", split.by = "Type_new_temp", raster=F, pt.size = 0.5, order=T) & 
  scale_color_manual(values = c("#009ACD"), na.value = "gray80") &
  ggeasy::easy_remove_axes()&
  theme(
    title = element_text(size=24, face = "bold"),
    strip.text.x = element_text(size=20, face = "bold"),
    legend.text = element_text(size=20)
  ) & 
  ggtitle("Tc17 clones in TIL product")
dev.off()


pdf("Graphs/CD8_Tc17_clones_full.pdf", height = 4, width = 18)
DimPlot(subset(gex.CD8, subset = Type_new_temp %in% c("pre-ACT", "Intermediate product", "TIL product", "PBMC 7dpt", "post-ACT")), group.by="is_Th17_ntb", split.by = "Type_new_temp", raster=F, pt.size = 0.5, order=T) & 
  scale_color_manual(values = c("#009ACD"), na.value = "gray80") &
  ggeasy::easy_remove_axes()&
  theme(
    title = element_text(size=24, face = "bold"),
    strip.text.x = element_text(size=20, face = "bold"),
    legend.text = element_text(size=20)
  ) & 
  ggtitle("Tc17 clones in TIL product")
dev.off()


## -----------------------------------------------------------------------------------------------------------
# Highlighting Th17 clones at different time points


Exp.CD4 <-  readRDS(file = "../saveFiles/CD4_Expanded_TILs_subclustering.rds")

gex.CD4@meta.data$Patient_ntb <- paste(gex.CD4@meta.data$Patient, gex.CD4@meta.data$ntb, sep = "_")

Exp.CD4@meta.data$Patient_ntb <- paste(Exp.CD4@meta.data$Patient, Exp.CD4@meta.data$ntb, sep = "_")
Tc17_ntb<- Exp.CD4@meta.data %>% filter(seurat_clusters == 2 & Type == "Expanded_TILs" & ntb != "NA" & !is.na(ntb) & nta != "NA" & !is.na(nta)) %>% pull(Patient_ntb)
non_Tc17_ntb<- Exp.CD4@meta.data %>% filter(seurat_clusters != 2 & Type == "Expanded_TILs" & ntb != "NA" & !is.na(ntb) & nta != "NA" & !is.na(nta)) %>% pull(Patient_ntb)

Tc17_ntb<- Tc17_ntb[!(Tc17_ntb %in% non_Tc17_ntb)]

gex.CD4$Type_new_temp <- ifelse(gex.CD4$Type_new %in% c("post-ACT.1", "post-ACT.2"), "post-ACT", as.character(gex.CD4$Type_new))
gex.CD4$Type_new_temp <- factor(gex.CD4$Type_new_temp, levels = c("pre-ACT", "Intermediate product", "TIL product", "PBMC 7dpt", "post-ACT"))


gex.CD4$is_Th17_ntb <- ifelse(gex.CD4$Patient_ntb %in% Tc17_ntb, "Th17", NA)

pdf("Fig/CD4_Th17_clones.pdf", height = 4, width = 13)
DimPlot(subset(gex.CD4, subset = Type_new_temp %in% c("pre-ACT", "Intermediate product", "TIL product")), group.by="is_Th17_ntb", split.by = "Type_new_temp", raster=F, pt.size = 0.5, order=T) & 
  scale_color_manual(values = c("#009ACD"), na.value = "gray80") &
  ggeasy::easy_remove_axes()&
  theme(
    title = element_text(size=24, face = "bold"),
    strip.text.x = element_text(size=20, face = "bold"),
    legend.text = element_text(size=20)
  ) & 
  ggtitle("Th17 clones in TIL product")
dev.off()


pdf("Graphs/CD4_Th17_clones_full.pdf", height = 4, width = 18)
DimPlot(subset(gex.CD4, subset = Type_new_temp %in% c("pre-ACT", "Intermediate product", "TIL product", "PBMC 7dpt", "post-ACT")), group.by="is_Th17_ntb", split.by = "Type_new_temp", raster=F, pt.size = 0.5, order=T) & 
  scale_color_manual(values = c("#009ACD"), na.value = "gray80") &
  ggeasy::easy_remove_axes()&
  theme(
    title = element_text(size=24, face = "bold"),
    strip.text.x = element_text(size=20, face = "bold"),
    legend.text = element_text(size=20)
  ) & 
  ggtitle("Th17 clones in TIL product")
dev.off()


