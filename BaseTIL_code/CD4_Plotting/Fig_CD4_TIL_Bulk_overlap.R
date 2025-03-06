## -----------------------------------------------------------------------------------------------------------
#load dependencies

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(scRepertoire)
  library(readr)
  library(dplyr)
  library(AUCell)
  library(UCell)
  library(ggeasy)
  library(ggpubr)
  library(EnhancedVolcano)
  library(mixtools)
  library(dittoSeq)
  library(immunarch)
  library(stringr)
})

set.seed(1234567)
setwd("~/BaseTIL_code/CD4_Plotting/")

gex.CD4 <- readRDS(file = "../saveFiles/CD4_annotated.rds")

## -----------------------------------------------------------------------------------------------------------
#load TCR data

file_path = paste0("~/BaseTIL_code/TCR_activation_assay/", c("UPN001_3", "UPN002_2", "UPN003", "UPN006", "UPN008", "UPN009", "UPN011_2"), "/results")

immdata <- repLoad(file_path, .mode = "single")

plt <- list()

Sample_Name_vec <- c("UPN001 Tumor", "UPN001 PreREP", "UPN001 Expanded_TILs", "UPN001 PBMC_7dpt",  "UPN001 Rebiopsy1", "UPN001 Rebiopsy2",
                     "UPN002 Tumor", "UPN002 PreREP", "UPN002 Expanded_TILs", "UPN002 PBMC_7dpt",  
                     "UPN003 Tumor", "UPN003 PreREP", "UPN003 Expanded_TILs", "UPN003 PBMC_7dpt",  
                     "UPN006 Tumor", "UPN006 PreREP", "UPN006 Expanded_TILs", "UPN006 PBMC_7dpt",  "UPN006 Rebiopsy1", "UPN006 Rebiopsy2",
                     "UPN008 Tumor", "UPN008 PreREP", "UPN008 Expanded_TILs", "UPN008 PBMC_7dpt",  "UPN008 Rebiopsy1","UPN008 Rebiopsy2",
                     "UPN009 Tumor", "UPN009 PreREP", "UPN009 Expanded_TILs", "UPN009 PBMC_7dpt",  "UPN009 Rebiopsy1", "UPN009 Rebiopsy2",
                     "UPN011 Tumor", "UPN011 PreREP", "UPN011 Expanded_TILs", "UPN011 PBMC_7dpt",  "UPN011 Rebiopsy1")

for (Sample_Name_iterater in Sample_Name_vec) {

  patient <-  str_split_fixed(Sample_Name_iterater, ' ', 2)
  if(patient[1,1] == "UPN001") patient[1,1] <- "UPN001_3"
  if(patient[1,1] == "UPN002") patient[1,1] <- "UPN002_2"
  if(patient[1,1] == "UPN011") patient[1,1] <- "UPN011_2"
  
  concat<- repFilter(immdata, .method = "by.meta", .query = list(Patient_Type = include(paste0(patient[1,1], "_TIL"))))
  
  concat <- bind_rows(concat[["data"]], .id = "column_label")
  
  
  gex.CD4$Tracked_clones <- "None"
  gex.CD4$Tracked_clones <- gex.CD4$ntb %in%  concat$CDR3.nt
    

  plt[[paste(Sample_Name_iterater)]] <-
    DimPlot(
      subset(gex.CD4, Sample_Name == Sample_Name_iterater),
      label = F,
      cells.highlight = WhichCells(gex.CD4, expression = Tracked_clones == TRUE),
      sizes.highlight = 0.5, raster=FALSE
    ) +
    scale_color_manual(
      labels = c("Unselected", "NeoTCRs"),
      values = c("grey90", "#F8766D")
    ) + NoLegend() + easy_remove_axes() + ggtitle(paste(Sample_Name_iterater)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  
}

plt2 <- list(ggarrange(plotlist = plt[c("UPN001 Tumor", "UPN001 PreREP", "UPN001 Expanded_TILs", "UPN001 PBMC_7dpt",  "UPN001 Rebiopsy1", "UPN001 Rebiopsy2")] , ncol = 6, nrow = 1),
             ggarrange(plotlist = plt[c("UPN002 Tumor", "UPN002 PreREP", "UPN002 Expanded_TILs", "UPN002 PBMC_7dpt")] , ncol = 6, nrow = 1),
             ggarrange(plotlist = plt[c("UPN003 Tumor", "UPN003 PreREP", "UPN003 Expanded_TILs", "UPN003 PBMC_7dpt")] , ncol = 6, nrow = 1),
             ggarrange(plotlist = plt[c("UPN006 Tumor", "UPN006 PreREP", "UPN006 Expanded_TILs", "UPN006 PBMC_7dpt",  "UPN006 Rebiopsy1", "UPN006 Rebiopsy2")] , ncol = 6, nrow = 1),
             ggarrange(plotlist = plt[c("UPN008 Tumor", "UPN008 PreREP", "UPN008 Expanded_TILs", "UPN008 PBMC_7dpt",  "UPN008 Rebiopsy1",  "UPN008 Rebiopsy2")] , ncol = 6, nrow = 1),
             ggarrange(plotlist = plt[c("UPN009 Tumor", "UPN009 PreREP", "UPN009 Expanded_TILs", "UPN009 PBMC_7dpt",  "UPN009 Rebiopsy1", "UPN009 Rebiopsy2")] , ncol = 6, nrow = 1),
             ggarrange(plotlist = plt[c("UPN011 Tumor", "UPN011 PreREP", "UPN011 Expanded_TILs", "UPN011 PBMC_7dpt",  "UPN011 Rebiopsy1")] , ncol = 6, nrow = 1))


png("Graphs/CD4_inBulkTIL.png", width = 1200, height = 1400, type = "cairo")
ggarrange(plotlist = plt2 , ncol = 1, nrow = 7)
dev.off()







plt <- list()
gex.CD4$Tracked_clones <- FALSE

Patient_vec <- c("UPN001","UPN002", "UPN003", "UPN006", "UPN008", "UPN009", "UPN011")

for (Patient_iterater in Patient_vec) {
  if(Patient_iterater == "UPN001") {pat_2 <- "UPN001_3"}
  else if(Patient_iterater == "UPN002") {pat_2 <- "UPN002_2"}
  else if(Patient_iterater == "UPN011") {pat_2 <- "UPN011_2"} else {pat_2 <- Patient_iterater}
  
  x.df <-
    as.data.frame(repFilter(
      immdata,
      .method = "by.meta",
      .query = list(Patient_Type = include(paste0(pat_2, "_TIL")))
    )$data[[1]])
  x.df <- x.df[, c("CDR3.nt", "Clones", "Proportion")]
  
  gex.CD4@meta.data$Tracked_clones <- ifelse(gex.CD4@meta.data$ntb %in%  x.df$CDR3.nt & gex.CD4@meta.data$Patient == Patient_iterater, Patient_iterater, gex.CD4$Tracked_clones)
} 

table(gex.CD4@meta.data$Tracked_clones)




pdf("Fig/TIL_PBMC_overlay.pdf", width = 4, height = 4)

DimPlot(
  subset(gex.CD4, Type == "PBMC_7dpt"),
  label = F,
  cells.highlight = WhichCells(gex.CD4, expression = Tracked_clones != FALSE),
  sizes.highlight = 0.75, raster=FALSE
) +
  scale_color_manual(
    labels = c("Unselected", "NeoTCRs"),
    values = c("grey80",  "#009ACD")
  ) + NoLegend() + easy_remove_axes() + ggtitle("Shared clones with \n TIL product (bulk)") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()



