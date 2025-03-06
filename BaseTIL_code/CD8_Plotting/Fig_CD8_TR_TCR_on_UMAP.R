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
setwd("~/BaseTIL_code/CD8_Plotting")

gex.CD8 <- readRDS(file = "../saveFiles/CD8_annotated.rds")


## -----------------------------------------------------------------------------------------------------------
# Highlight overall reactive clones in UAMP 


gex.CD8$reactive_temp <- ifelse(gex.CD8$Tumor_reactive == TRUE, gex.CD8$Patient, NA)
gex.CD8$reactive_temp <- factor(gex.CD8$reactive_temp, levels = rev(c("UPN001", "UPN002", "UPN003", "UPN006", "UPN011", "UPN009")))

gex.CD8$prepost <-  ifelse(as.character(gex.CD8$Type_new) %in% c("post-ACT.1", "post-ACT.2"), "post-ACT", as.character(gex.CD8$Type_new))
gex.CD8$prepost <- factor(gex.CD8$prepost, level = c("pre-ACT", "Intermediate product", "TIL product", "PBMC 7dpt", "post-ACT"))



gex.CD8$reactive_temp <- ifelse(gex.CD8$overall_reactive == TRUE, gex.CD8$Patient, NA)
gex.CD8$reactive_temp <- factor(gex.CD8$reactive_temp, levels = rev(c("UPN001", "UPN002", "UPN003", "UPN006", "UPN011", "UPN008", "UPN009")))


pdf("Fig/CD8_DimPlot_Reactive_Overall.pdf", height = 4, width = 22)
DimPlot(subset(gex.CD8, Type != "PBMC_Ctrl") , 
        group.by = "reactive_temp", 
        split.by = "prepost", raster = F,  order = T, pt.size = 0.5
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
  ggtitle(NULL) & theme(strip.text.x = element_text(size=20, face = "bold"),
                        legend.text = element_text(size=20),
                        legend.title = element_text(size=20, face = "bold"))
dev.off()

