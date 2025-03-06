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
  library(scRepertoire)
})

set.seed(1234567)
setwd("~/BaseTIL_code")

gex.all <- readRDS(file = "saveFiles/BaseTIL_Integrated_Tcells.rds")



## -----------------------------------------------------------------------------------------------------------
# Plot percent Richness of TCRs 


gex.all$Type <- as.character(gex.all$Type)
gex.all$Sample_Name <- as.character(gex.all$Sample_Name)



quant <- clonalQuant(gex.all, 
                     cloneCall="nt", 
                     chain = "TRB", 
                     scale = TRUE, group.by = "Sample_Name", exportTable = T)


quant[c('Patient', 'Type')] <- str_split_fixed(quant$Sample_Name, ' ', 2)

pdf("Fig/percRichness_Expansion.pdf", height = 4, width = 3.2)
quant[quant$Type %in% c("Tumor", "PreREP", "Expanded_TILs"),] %>%
  ggplot(aes(factor(Type, 
                    levels = c("Tumor","PreREP","Expanded_TILs"),
                    labels= c("pre-ACT","Intermediate product","TIL product")
                    ), scaled, color = Patient, group = Patient)) +
  geom_point(size=2) + geom_line() +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 13, colour = "black")) +
  labs(
    x = "",
    y = "Percent richness"
  ) +  ylim(c(0, 100))+
  scale_color_manual(name = "Patients", values = c(
    "UPN001" =  "#BC3C29B2",
    "UPN002" = "#0072B5FF",
    "UPN003" = "#E18727B2",
    "UPN006" = "#20854EB2",
    "UPN008" = "#7876B1B2",
    "UPN009" = "#6F99ADB2",
    "UPN011" = "#EE4C97B2"
  )) +theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



pdf("Fig/percRichness_prePost.pdf", height = 4, width = 3.2)
quant[quant$Type %in% c("Tumor", "Rebiopsy1", "Rebiopsy2") & quant$Patient %in% c("UPN001", "UPN006", "UPN008", "UPN009", "UPN011") & !(quant$Sample_Name %in% c("UPN008 Rebiopsy1")),] %>%
  ggplot(aes(factor(Type, 
                    levels = c("Tumor", "Rebiopsy1", "Rebiopsy2"),
                    labels= c("pre-ACT", "post-ACT.1", "post-ACT.2")
                    ), scaled, color = Patient, group = Patient)) +
  geom_point(size=2) + geom_line() +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 13, colour = "black")) +
  labs(
    x = "",
    y = "Percent richness"
  ) +  ylim(c(0, 100))+
  scale_color_manual(name = "Patients", values = c(
    "UPN001" =  "#BC3C29B2",
    "UPN002" = "#0072B5FF",
    "UPN003" = "#E18727B2",
    "UPN006" = "#20854EB2",
    "UPN008" = "#7876B1B2",
    "UPN009" = "#6F99ADB2",
    "UPN011" = "#EE4C97B2"
  )) +theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

