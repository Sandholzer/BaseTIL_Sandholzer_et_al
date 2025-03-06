## -----------------------------------------------------------------------------------------------------------
#loading dependencies

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)

})

set.seed(1234567)
setwd("~/BaseTIL_code")
merged.gex <- readRDS(file = "saveFiles/BaseTIL_raw.rds")


## -----------------------------------------------------------------------------------------------------------
# calculate some QC metrics

merged.gex <-
  PercentageFeatureSet(merged.gex, "MT-", col.name = "percent_mito")
merged.gex <-
  PercentageFeatureSet(merged.gex, "^RP[SL]", col.name = "percent_ribo")

merged.gex <-
  PercentageFeatureSet(merged.gex, "^HB[^(P)]", col.name = "percent_hb")

merged.gex <-
  PercentageFeatureSet(merged.gex, "PECAM1|PF4", col.name = "percent_plat")




## -----------------------------------------------------------------------------------------------------------
#Plotting of Parameters for Filter selection
feats <-
  c("nFeature_RNA",
    "nCount_RNA",
    "percent_mito",
    "percent_ribo",
    "percent_hb", "percent_plat")
png(file="Graphs/QC_features_bevore_removal.png", width=1000, height=800 , type = "cairo")
VlnPlot(
  merged.gex,
  group.by = "Sample_Name",
  features = feats,
  pt.size = 0.1,
  ncol = 2
)+
  NoLegend()
dev.off()


# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
png(file="Graphs/QC_featuresVSumi_bevore_removal.png", width=500, height=1000, type = "cairo")
merged.gex@meta.data %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent_mito)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 1250) +
  geom_hline(yintercept = 650) +
  facet_wrap(~Sample_Name)
dev.off()

merged.gex$log10GenesPerUMI <- log10(merged.gex$nFeature_RNA) / log10(merged.gex$nCount_RNA)

png(file="Graphs/QC_log10GenesPerUMI_bevore_removal.png", width=600, height=500, type = "cairo")
merged.gex@meta.data %>%
  ggplot(aes(x=log10GenesPerUMI, color = Sample_Name, fill=Sample_Name)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)+
  facet_wrap(~Sample_Name)
dev.off()


## -----------------------------------------------------------------------------------------------------------
# filter the object to remove low quality files

table(merged.gex$Sample_Name)
merged.gex <- subset(merged.gex, subset = nFeature_RNA > 650 & nCount_RNA > 1250 & percent_mito < 15  & percent_hb < 10 & percent_ribo > 5 )
table(merged.gex$Sample_Name)

dim(merged.gex)

## -----------------------------------------------------------------------------------------------------------

# Only keeping those genes expressed in more than 5 cells
selected_f <- rownames(merged.gex)[Matrix::rowSums(merged.gex) > 5]
merged.gex <- subset(merged.gex, features = selected_f)


# # Filter non-coding genes
merged.gex <- merged.gex[!grepl("MALAT1", rownames(merged.gex)), ]
merged.gex <- merged.gex[!grepl("NEAT1", rownames(merged.gex)), ]
merged.gex <- merged.gex[!grepl("XIST", rownames(merged.gex)), ]

dim(merged.gex)


## -----------------------------------------------------------------------------------------------------------
print("QC Done!")
saveRDS(merged.gex, "saveFiles/BaseTIL_QC.rds", compress = "gzip")

