## -----------------------------------------------------------------------------------------------------------
#load dependencies

suppressPackageStartupMessages({
  library(Seurat)
  library(scRepertoire)
  library(tidyverse)
  library(Matrix)
  library(immunarch)
  library(viridis)
  
})

set.seed(1234567)
setwd("~/BaseTIL_code")



## -----------------------------------------------------------------------------------------------------------
# Function to To add information from clonotype.csv
add_clonotype <- function(tcr_folder, seurat_obj){
  tcr <- read.csv(paste(tcr_folder,"filtered_contig_annotations.csv", sep=""))
  
  # Subsets so only the first line of each barcode is kept,
  # as each entry for given barcode will have same clonotype.
  tcr <- tcr[!duplicated(tcr$barcode), ]
  
  # Only keep the barcode and clonotype columns. 
  # We'll get additional clonotype info from the clonotype table.
  tcr <- tcr[,c("barcode", "raw_clonotype_id")]
  names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"
  
  # Clonotype-centric info.
  clono <- read.csv(paste(tcr_folder,"clonotypes.csv", sep=""))
  
  # Slap the AA sequences onto our original table by clonotype_id.
  tcr <- merge(tcr, clono[, c("clonotype_id", "cdr3s_aa", "inkt_evidence", "mait_evidence")])
  
  # Renames rownames as barcode
  rownames(tcr) <-  tcr[,"barcode"]
  
  # Add to the Seurat object's metadata.
  clono_seurat <- AddMetaData(object=seurat_obj, metadata=tcr)
  return(clono_seurat)
}

## -----------------------------------------------------------------------------------------------------------
# load gene expression data from cellranger output


meta_data <- read.csv("~/BaseTIL_code/data/BaseTIL_metadata_table.csv")

sample_chr_vector <- c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", 
                       "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", "S19", 
                       "S20", "S21", "S22", "S23", "S24", "S25", "S26", "S27", "S28", "S29", 
                       "S30", "S32", "S33", "S34", "S35", "S36", "S37", "S39", 
                       "S40", "S41", "S42") 


for (sample in sample_chr_vector) {
  path_to_file <- paste0("~/BaseTIL/Cellranger_multi/run_Cellranger/", sample, "/outs")
  toc <-  Seurat::Read10X(paste0(path_to_file, "/per_sample_outs/", sample, "/count/sample_filtered_feature_bc_matrix"))
  sample.gex <- Seurat::CreateSeuratObject(counts = toc, project = paste(sample, ".gex", sep = ""))
  sample.gex <- add_clonotype(tcr_folder= paste("~/BaseTIL/Cellranger_multi/run_Cellranger/", sample, "/outs/per_sample_outs/", sample, "/vdj_t/", sep=""), seurat_obj=sample.gex)
  sample.gex$Patient <- meta_data[meta_data$Sample==sample, "Patient"]
  sample.gex$Type <- meta_data[meta_data$Sample==sample, "Type"]
  sample.gex$Response <- meta_data[meta_data$Sample==sample, "Response"]
  sample.gex$Source <- meta_data[meta_data$Sample==sample, "Source"]
  sample.gex$batchV <- meta_data[meta_data$Sample==sample, "Batch"]
  name_csv <- paste(sample, ".gex", sep = "")
  assign(name_csv, sample.gex)
}


## -----------------------------------------------------------------------------------------------------------
#merge single files to one

for (SX.gex in c(S1.gex, S2.gex, S3.gex, S4.gex, S5.gex, S6.gex, S7.gex, S8.gex, S9.gex, S10.gex, S11.gex, S12.gex, S13.gex, S14.gex, S15.gex, S16.gex, S17.gex, S18.gex, S19.gex, S20.gex, S21.gex, S22.gex, S23.gex, S24.gex, S25.gex, S26.gex, S27.gex, S28.gex, S29.gex, 
                 S30.gex, S32.gex, S33.gex, S34.gex, S35.gex, S36.gex, S37.gex, S39.gex, S40.gex, S41.gex, S42.gex)) {
  SX.gex <- NormalizeData(SX.gex, verbose = FALSE)
}

merged.gex <- merge(x = S1.gex, y= c(S2.gex, S3.gex, S4.gex, S5.gex, S6.gex, S7.gex,S8.gex, S9.gex, S10.gex, S11.gex, S12.gex, S13.gex, S14.gex, S15.gex, S16.gex, S17.gex, S18.gex, S19.gex, S20.gex, S21.gex, S22.gex, S23.gex, S24.gex, S25.gex, S26.gex, S27.gex, S28.gex, S29.gex, 
                                     S30.gex, S32.gex, S33.gex, S34.gex, S35.gex, S36.gex, S37.gex, S39.gex, S40.gex, S41.gex, S42.gex), 
                    add.cell.ids = sample_chr_vector)

merged.gex$Sample_Name <- paste(merged.gex$Patient, merged.gex$Type)




# Type column uses internal nomenclature without spaces, Tumor corresponds to baseline Tumor
# Type_new uses the nomenclature used in the manuscript
# Sample_Name is a combination between Patient and Type


merged.gex$Type <- factor(merged.gex$Type, levels = c("Tumor", "PreREP", "Expanded_TILs", "PBMC_7dpt", "Rebiopsy1", "Rebiopsy2", "PBMC_Ctrl"))

merged.gex$Type_new <- factor(merged.gex$Type, 
                           levels = c("Tumor","PreREP","Expanded_TILs", "PBMC_7dpt", "Rebiopsy1", "Rebiopsy2"),
                           labels= c("pre-ACT","Intermediate product","TIL product", "PBMC 7dpt", "post-ACT.1", "post-ACT.2"))


merged.gex$Sample_Name <- factor(merged.gex$Sample_Name, levels = c("UPN001 Tumor", "UPN001 PreREP", "UPN001 Expanded_TILs", "UPN001 PBMC_7dpt",  "UPN001 Rebiopsy1", "UPN001 Rebiopsy2",
                                                              "UPN002 Tumor", "UPN002 PreREP", "UPN002 Expanded_TILs", "UPN002 PBMC_7dpt",
                                                              "UPN003 Tumor", "UPN003 PreREP", "UPN003 Expanded_TILs", "UPN003 PBMC_7dpt",
                                                              "UPN006 Tumor", "UPN006 PreREP", "UPN006 Expanded_TILs", "UPN006 PBMC_7dpt",  "UPN006 Rebiopsy1", "UPN006 Rebiopsy2",
                                                              "UPN008 Tumor", "UPN008 PreREP", "UPN008 Expanded_TILs", "UPN008 PBMC_7dpt",  "UPN008 Rebiopsy1","UPN008 Rebiopsy2",
                                                              "UPN009 Tumor", "UPN009 PreREP", "UPN009 Expanded_TILs", "UPN009 PBMC_7dpt",  "UPN009 Rebiopsy1", "UPN009 Rebiopsy2",
                                                              "UPN011 Tumor", "UPN011 PreREP", "UPN011 Expanded_TILs", "UPN011 PBMC_7dpt",  "UPN011 Rebiopsy1"))


merged.gex$Sample_Name_new <- paste(merged.gex$Patient, merged.gex$Type_new)

merged.gex$Sample_Name_new <- factor(merged.gex$Sample_Name_new, levels = c("UPN001 pre-ACT", "UPN001 Intermediate product", "UPN001 TIL product", "UPN001 PBMC 7dpt",  "UPN001 post-ACT.1", "UPN001 post-ACT.2",
                                                                      "UPN002 pre-ACT", "UPN002 Intermediate product", "UPN002 TIL product", "UPN002 PBMC 7dpt",
                                                                      "UPN003 pre-ACT", "UPN003 Intermediate product", "UPN003 TIL product", "UPN003 PBMC 7dpt",
                                                                      "UPN006 pre-ACT", "UPN006 Intermediate product", "UPN006 TIL product", "UPN006 PBMC 7dpt",  "UPN006 post-ACT.1", "UPN006 post-ACT.2",
                                                                      "UPN008 pre-ACT", "UPN008 Intermediate product", "UPN008 TIL product", "UPN008 PBMC 7dpt",  "UPN008 post-ACT.1","UPN008 post-ACT.2",
                                                                      "UPN009 pre-ACT", "UPN009 Intermediate product", "UPN009 TIL product", "UPN009 PBMC 7dpt",  "UPN009 post-ACT.1", "UPN009 post-ACT.2",
                                                                      "UPN011 pre-ACT", "UPN011 Intermediate product", "UPN011 TIL product", "UPN011 PBMC 7dpt",  "UPN011 post-ACT.1"))




# prepost column combines post-ACT.1 and post-ACT.2 into post-ACT for visualizations where combination is desired
merged.gex$prepost <-  ifelse(as.character(merged.gex$Type_new) %in% c("post-ACT.1", "post-ACT.2"), "post-ACT", as.character(merged.gex$Type_new))
merged.gex$prepost <- factor(merged.gex$prepost, level = c("pre-ACT", "Intermediate product", "TIL product", "PBMC 7dpt", "post-ACT"))


rm(list = ls()[grep("^S", ls())])

# run garbage collect to free up memory
gc()



## -----------------------------------------------------------------------------------------------------------
# load tcr data from cellranger output



for (sample in sample_chr_vector) {
  path_to_file <- paste("~/BaseTIL/Cellranger_multi/run_Cellranger/", sample, "/outs/per_sample_outs/", sample, "/vdj_t/filtered_contig_annotations.csv", sep="")
  sample.data <- read.csv(path_to_file)
  name_csv <- paste(sample)
  assign(name_csv, sample.data)
}

contig_list <- list(S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11, S12, S13, S14, S15, S16, S17, S18, S19, S20, S21, S22, S23, S24, S25, S26, S27, S28, S29, 
                    S30, S32, S33, S34, S35, S36, S37, S39, S40, S41, S42)

combined <- combineTCR(
  contig_list,
  samples = sample_chr_vector, 
  ID = sample_chr_vector,
  removeMulti = F,
  filterMulti = TRUE,
  removeNA = F
)

for (SX in sample_chr_vector) {
  combined[[paste(SX, "_", SX, sep = "")]][["barcode"]] <- gsub(paste(SX, "_", SX, sep = ""), SX, combined[[paste(SX, "_", SX, sep = "")]][["barcode"]])
  
}

# merge tcr and gex objects

merged.gex <- combineExpression(combined, merged.gex, 
                                cloneCall="nt", 
                                chain = "TRB",
                                proportion = TRUE,
                                #cloneSize=c(Single=1, Small=10, Medium=50, Large=200, Hyperexpanded=5000)
)

rm(list = ls()[grep("^S", ls())])
rm(contig_list, combined)

#Annotates doublets based on TCRb gene
TCR <- as.data.frame(read.csv("~/BaseTIL_code/saveFiles/TCR_only_annotation.csv", row.names = 1))
merged.gex<- Seurat::AddMetaData(merged.gex, TCR$chain_pairing, col.name = "chain_pairing")

merged.gex@meta.data <- merged.gex@meta.data %>% 
  separate(CTaa, into = c("CDR3a", "CDR3b"), sep="_", remove = F)

merged.gex@meta.data <- merged.gex@meta.data %>% 
  separate(CTnt, into = c("nta", "ntb"), sep="_", remove = F)



## -----------------------------------------------------------------------------------------------------------
# load bulk TCRb seq data from co-culturing assay to identify tumor-reactive T cells


file_path = paste0("~/BaseTIL/TCR_activation_assay/", c("UPN001", "UPN001_2","UPN001_3" , "UPN002", "UPN002_2", "UPN003", "UPN006", "UPN008", "UPN009", "UPN011", "UPN011_2"), "/results")

immdata <- repLoad(file_path, .mode = "single")


# function calculates logFC enrichment of sorted vs. bulk TIL product

enrichmentPlot <-  function (df = immdata,
                             base = NULL,
                             enrich = NULL,
                             exportTable = FALSE,
                             threshold) {
  x.df <-
    as.data.frame(repFilter(
      immdata,
      .method = "by.meta",
      .query = list(Patient_Type = include(base))
    )$data[[1]])
  x.df <- x.df[, c("CDR3.nt", "Clones", "Proportion")]
  
  x.df <- aggregate(Clones ~ CDR3.nt, data = x.df, sum)
  x.df$Proportion <- x.df$Clones / sum(x.df$Clones)
  colnames(x.df)[2] <- base
  
  y.df <-
    as.data.frame(repFilter(
      immdata,
      .method = "by.meta",
      .query = list(Patient_Type = include(enrich))
    )$data[[1]])
  y.df <- y.df[, c("CDR3.nt", "Clones", "Proportion")]
  y.df <- aggregate(Clones ~ CDR3.nt, data = y.df, sum)
  y.df$Proportion <- y.df$Clones / sum(y.df$Clones)
  colnames(y.df)[2] <- enrich
  
  combined.df <- merge(x.df, y.df, by = "CDR3.nt",  all = F ) # all.y  = T,
  
  combined.df[base][is.na(combined.df[base])] <- 0
  combined.df[enrich][is.na(combined.df[enrich])] <- 0
  
  
  rownames(combined.df) <- combined.df$CDR3.nt
  combined.df$CDR3.nt <- NULL
  
  combined.df$logFC <- log2(combined.df$Proportion.y/combined.df$Proportion.x)
  #combined.df$logCPM <- log2((combined.df$UPN003_TIL + combined.df$UPN003_DIG))
  combined.df$logCPM <- edgeR::cpm(y=(combined.df[,base] + combined.df[,enrich]), log=T)
  
  combined.df$Enrichment <-  ifelse(combined.df$logFC >= threshold , "Reactive", NA)
  
  if (exportTable == TRUE) {
    return(combined.df)
  } else {
    plot<- ggplot(combined.df, aes(x = logCPM, y = logFC, color = Enrichment)) +
      geom_point(alpha = 1, size = 1) +
      scale_color_manual(values=c("red", "gray70")) +
      geom_hline(
        yintercept = c(threshold),
        linetype = "dashed",
        color = "gray"
      ) +
      labs(title = "", x = "Average Expression (log2 CPM)", y = "Log Fold Change") +
      # Add labels for logFC > 4
      theme_classic()    
    
    return(plot)
  }
  
}



# Threshold was increased for cases where few positive cells where sorted to reduce impurities

# UPN001_3
set_threshold <- 1
enrichmentPlot(immdata, base = "UPN001_3_TIL", enrich = "UPN001_3_CL", threshold = set_threshold)
enrich_UPN001_tumor<- enrichmentPlot(immdata, base = "UPN001_3_TIL", enrich = "UPN001_3_CL", threshold = set_threshold, exportTable = T)
merged.gex$logFC_UPN001_Tumor <- FALSE
merged.gex$Reactive_UPN001_Tumor <- FALSE
merged.gex@meta.data[merged.gex$Patient == "UPN001",]$logFC_UPN001_Tumor <-  enrich_UPN001_tumor[, "logFC"][match(merged.gex@meta.data[merged.gex$Patient == "UPN001",]$ntb, rownames(enrich_UPN001_tumor))]
merged.gex@meta.data[merged.gex$Patient == "UPN001",]$Reactive_UPN001_Tumor <-  as.character(enrich_UPN001_tumor[, "Enrichment"][match(merged.gex@meta.data[merged.gex$Patient == "UPN001",]$ntb, rownames(enrich_UPN001_tumor))])



# UPN002_2
set_threshold <- 1
enrichmentPlot(immdata, base = "UPN002_2_TIL", enrich = "UPN002_2_FTD", threshold = set_threshold)
enrich_UPN002_tumor<- enrichmentPlot(immdata, base = "UPN002_2_TIL", enrich = "UPN002_2_FTD", threshold = set_threshold, exportTable = T)
merged.gex$logFC_UPN002_Tumor <- FALSE
merged.gex$Reactive_UPN002_Tumor <- FALSE
merged.gex@meta.data[merged.gex$Patient == "UPN002",]$logFC_UPN002_Tumor <-  enrich_UPN002_tumor[, "logFC"][match(merged.gex@meta.data[merged.gex$Patient == "UPN002",]$ntb, rownames(enrich_UPN002_tumor))]
merged.gex@meta.data[merged.gex$Patient == "UPN002",]$Reactive_UPN002_Tumor <-  as.character(enrich_UPN002_tumor[, "Enrichment"][match(merged.gex@meta.data[merged.gex$Patient == "UPN002",]$ntb, rownames(enrich_UPN002_tumor))])


# UPN003
set_threshold <- 1 #-4.5
enrichmentPlot(immdata, base = "UPN003_TIL", enrich = "UPN003_DIG", threshold = set_threshold)
enrich_UPN003_tumor<- enrichmentPlot(immdata, base = "UPN003_TIL", enrich = "UPN003_DIG", threshold = set_threshold, exportTable = T)
merged.gex$logFC_UPN003_Tumor <- FALSE
merged.gex$Reactive_UPN003_Tumor <- FALSE
merged.gex@meta.data[merged.gex$Patient == "UPN003",]$logFC_UPN003_Tumor <-  enrich_UPN003_tumor[, "logFC"][match(merged.gex@meta.data[merged.gex$Patient == "UPN003",]$ntb, rownames(enrich_UPN003_tumor))]
merged.gex@meta.data[merged.gex$Patient == "UPN003",]$Reactive_UPN003_Tumor <-  as.character(enrich_UPN003_tumor[, "Enrichment"][match(merged.gex@meta.data[merged.gex$Patient == "UPN003",]$ntb, rownames(enrich_UPN003_tumor))])


# UPN006
set_threshold <- 2
enrichmentPlot(immdata, base = "UPN006_TIL", enrich = "UPN006_DIG", threshold = set_threshold)
enrich_UPN006_tumor<- enrichmentPlot(immdata, base = "UPN006_TIL", enrich = "UPN006_DIG", threshold = set_threshold, exportTable = T)
merged.gex$logFC_UPN006_Tumor <- FALSE
merged.gex$Reactive_UPN006_Tumor <- FALSE
merged.gex@meta.data[merged.gex$Patient == "UPN006",]$logFC_UPN006_Tumor <-  enrich_UPN006_tumor[, "logFC"][match(merged.gex@meta.data[merged.gex$Patient == "UPN006",]$ntb, rownames(enrich_UPN006_tumor))]
merged.gex@meta.data[merged.gex$Patient == "UPN006",]$Reactive_UPN006_Tumor <-  as.character(enrich_UPN006_tumor[, "Enrichment"][match(merged.gex@meta.data[merged.gex$Patient == "UPN006",]$ntb, rownames(enrich_UPN006_tumor))])


# UPN008
set_threshold <- 1#-4.5
enrichmentPlot(immdata, base = "UPN008_TIL", enrich = "UPN008_R2DIG", threshold = set_threshold)
enrich_UPN008_R2<- enrichmentPlot(immdata, base = "UPN008_TIL", enrich = "UPN008_R2DIG", threshold = set_threshold, exportTable = T)
merged.gex$logFC_UPN008_R2 <- FALSE
merged.gex$Reactive_UPN008_R2 <- FALSE
merged.gex@meta.data[merged.gex$Patient == "UPN008",]$logFC_UPN008_R2 <-  enrich_UPN008_R2[, "logFC"][match(merged.gex@meta.data[merged.gex$Patient == "UPN008",]$ntb, rownames(enrich_UPN008_R2))]
merged.gex@meta.data[merged.gex$Patient == "UPN008",]$Reactive_UPN008_R2 <-  as.character(enrich_UPN008_R2[, "Enrichment"][match(merged.gex@meta.data[merged.gex$Patient == "UPN008",]$ntb, rownames(enrich_UPN008_R2))])


# UPN009
set_threshold <- 0 #-3 #changed form -4.5
enrichmentPlot(immdata, base = "UPN009_TIL", enrich = "UPN009_DIG", threshold = set_threshold)
enrich_UPN009_tumor<- enrichmentPlot(immdata, base = "UPN009_TIL", enrich = "UPN009_DIG", threshold = set_threshold, exportTable = T)
merged.gex$logFC_UPN009_Tumor <- FALSE
merged.gex$Reactive_UPN009_Tumor <- FALSE
merged.gex@meta.data[merged.gex$Patient == "UPN009",]$logFC_UPN009_Tumor <-  enrich_UPN009_tumor[, "logFC"][match(merged.gex@meta.data[merged.gex$Patient == "UPN009",]$ntb, rownames(enrich_UPN009_tumor))]
merged.gex@meta.data[merged.gex$Patient == "UPN009",]$Reactive_UPN009_Tumor <-  as.character(enrich_UPN009_tumor[, "Enrichment"][match(merged.gex@meta.data[merged.gex$Patient == "UPN009",]$ntb, rownames(enrich_UPN009_tumor))])


set_threshold <- 0 #-2
enrichmentPlot(immdata, base = "UPN009_TIL", enrich = "UPN009_R2L", threshold = set_threshold)
enrich_UPN009_R1<- enrichmentPlot(immdata, base = "UPN009_TIL", enrich = "UPN009_R2L", threshold = set_threshold, exportTable = T)
merged.gex$logFC_UPN009_R1 <- FALSE
merged.gex$Reactive_UPN009_R1 <- FALSE
merged.gex@meta.data[merged.gex$Patient == "UPN009",]$logFC_UPN009_R1 <-  enrich_UPN009_R1[, "logFC"][match(merged.gex@meta.data[merged.gex$Patient == "UPN009",]$ntb, rownames(enrich_UPN009_R1))]
merged.gex@meta.data[merged.gex$Patient == "UPN009",]$Reactive_UPN009_R1 <-  as.character(enrich_UPN009_R1[, "Enrichment"][match(merged.gex@meta.data[merged.gex$Patient == "UPN009",]$ntb, rownames(enrich_UPN009_R1))])


set_threshold <- 1 #-0.5
enrichmentPlot(immdata, base = "UPN009_TIL", enrich = "UPN009_R2K", threshold = set_threshold)
enrich_UPN009_R2<- enrichmentPlot(immdata, base = "UPN009_TIL", enrich = "UPN009_R2K", threshold = set_threshold, exportTable = T)
merged.gex$logFC_UPN009_R2 <- FALSE
merged.gex$Reactive_UPN009_R2 <- FALSE
merged.gex@meta.data[merged.gex$Patient == "UPN009",]$logFC_UPN009_R2 <-  enrich_UPN009_R2[, "logFC"][match(merged.gex@meta.data[merged.gex$Patient == "UPN009",]$ntb, rownames(enrich_UPN009_tumor))]
merged.gex@meta.data[merged.gex$Patient == "UPN009",]$Reactive_UPN009_R2 <-  as.character(enrich_UPN009_R2[, "Enrichment"][match(merged.gex@meta.data[merged.gex$Patient == "UPN009",]$ntb, rownames(enrich_UPN009_tumor))])

set_threshold <- 1
enrichmentPlot(immdata, base = "UPN011_TIL", enrich = "UPN011_B", threshold = set_threshold)
enrich_UPN011_Tumor_B<- enrichmentPlot(immdata, base = "UPN011_TIL", enrich = "UPN011_B", threshold = set_threshold, exportTable = T)
merged.gex$logFC_UPN011_Tumor_B <- FALSE
merged.gex$Reactive_UPN011_Tumor_B <- FALSE
merged.gex@meta.data[merged.gex$Patient == "UPN011",]$logFC_UPN011_Tumor_B <-  enrich_UPN011_Tumor_B[, "logFC"][match(merged.gex@meta.data[merged.gex$Patient == "UPN011",]$ntb, rownames(enrich_UPN011_Tumor_B))]
merged.gex@meta.data[merged.gex$Patient == "UPN011",]$Reactive_UPN011_Tumor_B <-  as.character(enrich_UPN011_Tumor_B[, "Enrichment"][match(merged.gex@meta.data[merged.gex$Patient == "UPN011",]$ntb, rownames(enrich_UPN011_Tumor_B))])


set_threshold <- 1
enrichmentPlot(immdata, base = "UPN011_2_TIL", enrich = "UPN011_2_SB", threshold = set_threshold)
enrich_UPN011_Tumor_SB<- enrichmentPlot(immdata, base = "UPN011_2_TIL", enrich = "UPN011_2_SB", threshold = set_threshold, exportTable = T)
merged.gex$logFC_UPN011_Tumor_SB<- FALSE
merged.gex$Reactive_UPN011_Tumor_SB<- FALSE
merged.gex@meta.data[merged.gex$Patient == "UPN011",]$logFC_UPN011_Tumor_SB<-  enrich_UPN011_Tumor_SB[, "logFC"][match(merged.gex@meta.data[merged.gex$Patient == "UPN011",]$ntb, rownames(enrich_UPN011_Tumor_SB))]
merged.gex@meta.data[merged.gex$Patient == "UPN011",]$Reactive_UPN011_Tumor_SB<-  as.character(enrich_UPN011_Tumor_SB[, "Enrichment"][match(merged.gex@meta.data[merged.gex$Patient == "UPN011",]$ntb, rownames(enrich_UPN011_Tumor_SB))])


set_threshold <- 1 #changed from -1
enrichmentPlot(immdata, base = "UPN011_TIL", enrich = "UPN011_R1CL", threshold = set_threshold)
enrich_UPN011_R1<- enrichmentPlot(immdata, base = "UPN011_TIL", enrich = "UPN011_R1CL", threshold = set_threshold, exportTable = T)
merged.gex$logFC_UPN011_R1 <- FALSE
merged.gex$Reactive_UPN011_R1 <- FALSE
merged.gex@meta.data[merged.gex$Patient == "UPN011",]$logFC_UPN011_R1 <-  enrich_UPN011_R1[, "logFC"][match(merged.gex@meta.data[merged.gex$Patient == "UPN011",]$ntb, rownames(enrich_UPN011_R1))]
merged.gex@meta.data[merged.gex$Patient == "UPN011",]$Reactive_UPN011_R1 <-  as.character(enrich_UPN011_R1[, "Enrichment"][match(merged.gex@meta.data[merged.gex$Patient == "UPN011",]$ntb, rownames(enrich_UPN011_R1))])



# annotation of tumor-reactive Tcells:
# selfreactive = T cells were detected in the same lesion they showed reactivity against
# Tumor_reactive = clones showed reactivity against baseline tumor (pre-ACT samples)
# Post_reactive = clones showed reactivity against post treatment lesions (post-ACT samples)
# overall_reactive = clone showed reactivity


merged.gex$selfreactive <- FALSE
merged.gex$selfreactive <- ifelse(
  merged.gex$Sample_Name == "UPN001 Tumor" &
    merged.gex$Reactive_UPN001_Tumor == "Reactive" |
    merged.gex$Sample_Name == "UPN001 Rebiopsy1" &
    merged.gex$Reactive_UPN001_Tumor == "Reactive" |
    merged.gex$Sample_Name == "UPN002 Tumor" &
    merged.gex$Reactive_UPN002_Tumor == "Reactive" |
    merged.gex$Sample_Name == "UPN003 Tumor" &
    merged.gex$Reactive_UPN003_Tumor == "Reactive" |
    merged.gex$Sample_Name == "UPN006 Tumor" &
    merged.gex$Reactive_UPN006_Tumor == "Reactive" |
    merged.gex$Sample_Name == "UPN008 Rebiopsy2" &
    merged.gex$Reactive_UPN008_R2 == "Reactive" |
    merged.gex$Sample_Name == "UPN009 Tumor" &
    merged.gex$Reactive_UPN009_Tumor == "Reactive" |
    merged.gex$Sample_Name == "UPN009 Rebiopsy1" &
    merged.gex$Reactive_UPN009_R1 == "Reactive" |
    merged.gex$Sample_Name == "UPN009 Rebiopsy2" &
    merged.gex$Reactive_UPN009_R2 == "Reactive" |
    merged.gex$Sample_Name == "UPN011 Tumor" &
    merged.gex$Reactive_UPN011_Tumor_SB == "Reactive" |
    merged.gex$Sample_Name == "UPN011 Tumor" &
    merged.gex$Reactive_UPN011_Tumor_B == "Reactive" |
    merged.gex$Sample_Name == "UPN011 Rebiopsy1" &
    merged.gex$Reactive_UPN011_R1 == "Reactive"
  ,
  TRUE,
  FALSE
)

merged.gex$Tumor_reactive <- FALSE
merged.gex$Tumor_reactive <- ifelse(
  merged.gex$Reactive_UPN001_Tumor == "Reactive" |
    merged.gex$Reactive_UPN002_Tumor == "Reactive" |
    merged.gex$Reactive_UPN003_Tumor == "Reactive" |
    merged.gex$Reactive_UPN006_Tumor == "Reactive" |
    merged.gex$Reactive_UPN011_Tumor_SB == "Reactive" |
    merged.gex$Reactive_UPN011_Tumor_B == "Reactive"
  ,
  TRUE,
  FALSE
)




merged.gex$Post_reactive <- FALSE
merged.gex$Post_reactive <- ifelse(
  merged.gex$Reactive_UPN008_R2 == "Reactive" |
    merged.gex$Reactive_UPN009_R1 == "Reactive" |
    merged.gex$Reactive_UPN009_R2 == "Reactive" |
    merged.gex$Reactive_UPN011_R1 == "Reactive"
  ,
  TRUE,
  FALSE
)


merged.gex$overall_reactive <- FALSE
merged.gex$overall_reactive <- ifelse(
  merged.gex$Reactive_UPN001_Tumor == "Reactive" |
    merged.gex$Reactive_UPN002_Tumor == "Reactive" |
    merged.gex$Reactive_UPN003_Tumor == "Reactive" |
    merged.gex$Reactive_UPN006_Tumor == "Reactive" |
    merged.gex$Reactive_UPN008_R2 == "Reactive" |
    merged.gex$Reactive_UPN009_Tumor == "Reactive" |
    merged.gex$Reactive_UPN009_R1 == "Reactive" |
    merged.gex$Reactive_UPN009_R2 == "Reactive" |
    merged.gex$Reactive_UPN011_Tumor_SB == "Reactive" |
    merged.gex$Reactive_UPN011_Tumor_B == "Reactive" |
    merged.gex$Reactive_UPN011_R1 == "Reactive"
  ,
  TRUE,
  FALSE
)


merged.gex@meta.data$Tumor_reactive[is.na(merged.gex@meta.data$Tumor_reactive)] <- FALSE
merged.gex@meta.data$Tumor_reactive <- factor(merged.gex@meta.data$Tumor_reactive, levels = c(TRUE, FALSE))
table(merged.gex$Tumor_reactive)

merged.gex@meta.data$overall_reactive[is.na(merged.gex@meta.data$overall_reactive)] <- FALSE
merged.gex@meta.data$overall_reactive <- factor(merged.gex@meta.data$overall_reactive, levels = c(TRUE, FALSE))
table(merged.gex$overall_reactive)

merged.gex@meta.data$Post_reactive[is.na(merged.gex@meta.data$Post_reactive)] <- FALSE
merged.gex@meta.data$Post_reactive <- factor(merged.gex@meta.data$Post_reactive, levels = c(TRUE, FALSE))
table(merged.gex$Post_reactive)



saveRDS(merged.gex, "saveFiles/BaseTIL_raw.rds", compress = "gzip")

