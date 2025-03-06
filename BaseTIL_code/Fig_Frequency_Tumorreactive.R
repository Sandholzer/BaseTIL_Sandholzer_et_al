## -----------------------------------------------------------------------------------------------------------
#load dependencies

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(stringr)
  library(ggplot2)
  library(ggalluvial)
  library(ggpubr)
  
})

set.seed(1234567)
setwd("~/BaseTIL_code")



gex.CD8 <- readRDS(file = paste0("saveFiles/CD8_annotated.rds"))

gex.CD4 <- readRDS(file = paste0("saveFiles/CD4_annotated.rds"))


metadata <- list()

metadata[["CD4_metadata"]] <- gex.CD4@meta.data
metadata[["CD8_metadata"]] <- gex.CD8@meta.data


saveRDS(metadata ,"saveFiles/Tcells_metadata.Rds")



## -----------------------------------------------------------------------------------------------------------
#Percentage of Tumor reactive T cells in therapy 

result_CD4 <- metadata[["CD4_metadata"]] %>%
  group_by(Sample_Name, overall_reactive, Response, Type, Patient) %>%
  summarize(count_CD4 = n()) %>%
  ungroup() %>%
  group_by(Sample_Name) %>%
  mutate(total_count_CD4 = sum(count_CD4), frequency_TR_CD4 = count_CD4/total_count_CD4)


result_CD8 <- metadata[["CD8_metadata"]] %>%
  group_by(Sample_Name, overall_reactive, Response, Type, Patient) %>%
  summarize(count_CD8 = n()) %>%
  ungroup() %>%
  group_by(Sample_Name) %>%
  mutate(total_count_CD8 = sum(count_CD8), frequency_TR_CD8 = count_CD8/total_count_CD8)

result_CD8 <- subset(result_CD8, overall_reactive == TRUE)
result_CD4 <- subset(result_CD4, overall_reactive == TRUE)

results <-  full_join(result_CD4,result_CD8, by = c("Sample_Name", "Response", "overall_reactive", "Type", "Patient"))

results[is.na(results)] <- 0
results$total_count <- results$total_count_CD4 + results$total_count_CD8
results$count <- results$count_CD4 + results$count_CD8
results$frequency_TR <- results$count/results$total_count

pdf("Graphs/Freq_reactive.pdf", height = 4, width = 5)
results %>%ggplot(aes(factor(Type, 
                             levels = c("Tumor","PreREP","Expanded_TILs", "PBMC_7dpt", "Rebiopsy1", "Rebiopsy2"),
                             labels= c("pre-ACT","Intermediate product","TIL product", "PBMC 7dpt", "post-ACT.1", "post-ACT.2")
                             ), frequency_TR, color = Patient, group = Patient)) +
  geom_point(size=2) + geom_line() +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 13, colour = "black")) +
  labs(
    x = "",
    y = "Percentage of \ntumor-reactive T cells"
  ) + scale_color_manual(name = "Patients", values = c(
    "UPN001" =  "#BC3C29B2",
    "UPN002" = "#0072B5FF",
    "UPN003" = "#E18727B2",
    "UPN006" = "#20854EB2",
    "UPN008" = "#7876B1B2",
    "UPN009" = "#6F99ADB2",
    "UPN011" = "#EE4C97B2"
  ))+ theme(axis.text.x = element_text(angle = 45,  hjust=1))  #+ scale_color_manual(values = c("#009ACD","#EE799F", "#CD2626")) 

results %>%ggplot(aes(factor(Type, 
                             levels = c("Tumor","PreREP","Expanded_TILs", "PBMC_7dpt", "Rebiopsy1", "Rebiopsy2"),
                             labels= c("pre-ACT","Intermediate product","TIL product", "PBMC 7dpt", "post-ACT.1", "post-ACT.2")
                             ), frequency_TR_CD8, color = Patient, group = Patient)) +
  geom_point(size=2) + geom_line() +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 13, colour = "black")) +
  labs(
    x = "",
    y = "Percentage of \nCD8+ tumor-reactive T cells"
  )+ scale_color_manual(name = "Patients", values = c(
    "UPN001" =  "#BC3C29B2",
    "UPN002" = "#0072B5FF",
    "UPN003" = "#E18727B2",
    "UPN006" = "#20854EB2",
    "UPN008" = "#7876B1B2",
    "UPN009" = "#6F99ADB2",
    "UPN011" = "#EE4C97B2"
  )) + theme(axis.text.x = element_text(angle = 45,  hjust=1))  #+ scale_color_manual(values = c("#009ACD","#EE799F", "#CD2626")) 

results %>%ggplot(aes(factor(Type, 
                             levels = c("Tumor","PreREP","Expanded_TILs", "PBMC_7dpt", "Rebiopsy1", "Rebiopsy2"),
                             labels= c("pre-ACT","Intermediate product","TIL product", "PBMC 7dpt", "post-ACT.1", "post-ACT.2")
                             ), frequency_TR_CD4, color = Patient, group = Patient)) +
  geom_point(size=2) + geom_line() +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 13, colour = "black")) +
  labs(
    x = "",
    y = "Percentage of \nCD4+ tumor-reactive T cells"
  ) + scale_color_manual(name = "Patients", values = c(
    "UPN001" =  "#BC3C29B2",
    "UPN002" = "#0072B5FF",
    "UPN003" = "#E18727B2",
    "UPN006" = "#20854EB2",
    "UPN008" = "#7876B1B2",
    "UPN009" = "#6F99ADB2",
    "UPN011" = "#EE4C97B2"
  ))+ theme(axis.text.x = element_text(angle = 45,  hjust=1))  #+ scale_color_manual(values = c("#009ACD","#EE799F", "#CD2626")) 
dev.off()








pdf("Fig/Freq_reactive_in_expansion.pdf", height = 4, width = 3.2)
results[results$Type %in% c("Tumor","PreREP","Expanded_TILs"),] %>%
  ggplot(aes(factor(Type, 
                    levels = c("Tumor","PreREP","Expanded_TILs"),
                    labels= c("pre-ACT","Intermediate product","TIL product")
                    ), frequency_TR, color = Patient, group = Patient)) +
  geom_point(size=2) + geom_line() +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 13, colour = "black")) +
  labs(
    x = "",
    y = "Percentage of \ntumor-reactive T cells"
  ) + scale_color_manual(name = "Patients", values = c(
    "UPN001" =  "#BC3C29B2",
    "UPN002" = "#0072B5FF",
    "UPN003" = "#E18727B2",
    "UPN006" = "#20854EB2",
    "UPN008" = "#7876B1B2",
    "UPN009" = "#6F99ADB2",
    "UPN011" = "#EE4C97B2"
  ))+ theme(axis.text.x = element_text(angle = 45,  hjust=1))  #+ scale_color_manual(values = c("#009ACD","#EE799F", "#CD2626")) 

results[results$Type %in% c("Tumor","PreREP","Expanded_TILs"),] %>%
  ggplot(aes(factor(Type, 
                    levels = c("Tumor","PreREP","Expanded_TILs"),
                    labels= c("pre-ACT","Intermediate product","TIL product")
                    ), frequency_TR_CD8, color = Patient, group = Patient)) +
  geom_point(size=2) + geom_line() +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 13, colour = "black")) +
  labs(
    x = "",
    y = "Percentage of \nCD8+ tumor-reactive T cells"
  )+ scale_color_manual(name = "Patients", values = c(
    "UPN001" =  "#BC3C29B2",
    "UPN002" = "#0072B5FF",
    "UPN003" = "#E18727B2",
    "UPN006" = "#20854EB2",
    "UPN008" = "#7876B1B2",
    "UPN009" = "#6F99ADB2",
    "UPN011" = "#EE4C97B2"
  )) + theme(axis.text.x = element_text(angle = 45,  hjust=1))  #+ scale_color_manual(values = c("#009ACD","#EE799F", "#CD2626")) 

results[results$Type %in% c("Tumor","PreREP","Expanded_TILs"),] %>%
  ggplot(aes(factor(Type, 
                    levels = c("Tumor","PreREP","Expanded_TILs"),
                    labels= c("pre-ACT","Intermediate product","TIL product")
                    ), frequency_TR_CD4, color = Patient, group = Patient)) +
  geom_point(size=2) + geom_line() +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 13, colour = "black")) +
  labs(
    x = "",
    y = "Percentage of \nCD4+ tumor-reactive T cells"
  ) + scale_color_manual(name = "Patients", values = c(
    "UPN001" =  "#BC3C29B2",
    "UPN002" = "#0072B5FF",
    "UPN003" = "#E18727B2",
    "UPN006" = "#20854EB2",
    "UPN008" = "#7876B1B2",
    "UPN009" = "#6F99ADB2",
    "UPN011" = "#EE4C97B2"
  ))+ theme(axis.text.x = element_text(angle = 45,  hjust=1))  #+ scale_color_manual(values = c("#009ACD","#EE799F", "#CD2626")) 
dev.off()



