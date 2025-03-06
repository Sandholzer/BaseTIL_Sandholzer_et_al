args <- commandArgs(trailingOnly = TRUE)
print(paste("Subset to analyze:", args))

## -----------------------------------------------------------------------------------------------------------
#load dependencies

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(UCell)
  library(rio)
})

set.seed(1234567)
setwd("~/BaseTIL_code")


Subset <- args
#Subset <- "CD4"

gex.list <- readRDS(file = paste0("saveFiles/", Subset, "_integrated_3.rds"))


## -----------------------------------------------------------------------------------------------------------
# plotting some signatures 

signaturesHumanTILs <- list()

signaturesHumanTILs$Tfh <- c("CD4"  ,  "CD40LG" ,"TOX2"   ,"MAF"   , "CD200",  "BATF"  )
signaturesHumanTILs$CD8_NaiveLike <- c(  "CCR7",  "IL7R",  "SELL",  "TCF7",  "S1PR1", "LEF1") #"CD8A",  "CD8B",
signaturesHumanTILs$CD8_EffectorMemory <- c(  "GZMA",  "GZMK" , "CCL5" , "CXCR3") #CD8A",  "CD8B",
signaturesHumanTILs$Thelper <- c("CD40LG", "CD4",    "IL7R",   "RORA",   "ANXA1")
signaturesHumanTILs$CD4_NaiveLike <- c("CD40LG", "CD4",    "CCR7",   "SELL",   "IL7R",   "TCF7",   "LEF1")
signaturesHumanTILs$CD8_Tpex <- c("LAG3",  "XCL1" , "CRTAM", "TOX",   "ZEB2",  "PDCD1", "TCF7" ,"CCR7") # "CD8A",  "CD8B",  
signaturesHumanTILs$CD8_Tex <- c("LAG3",   "HAVCR2", "GZMB",   "PRF1",   "PDCD1",  "TIGIT", "CTLA4" ) #"CD8A" ,  "CD8B",   
signaturesHumanTILs$Treg <- c("CD4",   "IL2RA", "FOXP3")
signaturesHumanTILs$CD8_FOXP3 <- c( "FOXP3", "CD4-") #"CD8A",   "CD8B",
signaturesHumanTILs$cycling <- c("TOP2A", "MKI67", "STMN1")
signaturesHumanTILs$Tcell_MAIT <- c("KLRB1", "SLC4A10", "NCR3")
signaturesHumanTILs$Tcell_gd <- c("TRDC", "TRGC1", "TRGC2", "TRDV1", "TRAC-", "TRBC1-", "TRBC2-")
signaturesHumanTILs$Tcell_NK <- c("FGFBP2", "SPON2", "KLRF1", "FCGR3A", "KLRD1", "TRDC", "CD3E-", "CD3G-")

signaturesHumanTILs$T_Naiv_2 <- c("TCF7", "LEF1", "CCR7", "SELL", "MAL")
signaturesHumanTILs$T_Memory <- c("IL7R", "GPR183", "ZFP36L2", "CXCR4", "ZNF683", "CD52", "HOPX")
signaturesHumanTILs$T_ResMemory <- c("ZNF683", "CD52", "HOPX", "ID2", "CXCR6", "XCL1", "XCL2")
signaturesHumanTILs$T_EffectorMemory <- c("GZMK", "CXCR5", "CCR4", "CD28", "CXCR3", "GZMH", "CD27", "HLA-DRB1","IFNG", "TNF")
signaturesHumanTILs$TEMRA <- c("TBX21", "ASCL2", "CX3CR1", "KLRG1")
signaturesHumanTILs$T_EX <- c("TOX", "LAG3", "CASP8", "TIGIT", "CTLA4", "ITGAE", "PDCD1", "CXCL13","CD27")


gex.list <- AddModuleScore_UCell(gex.list, features = signaturesHumanTILs)

featnames <- paste0(names(signaturesHumanTILs), "_UCell")


png(paste0("Graphs/", Subset,"_basic_Annotation.png"), width = 2000, height = 2000 , type = "cairo")
FeaturePlot(gex.list, features = featnames, pt.size = 0.1, order = T, min.cutoff = 0.4, raster = F)
dev.off()




Jerby<- import_list("data/Jerby-Arnon_markers.xlsx")
GS.Jerby <- as.list(Jerby$TableS3A_cell.type.markers)
GS.Jerby <-
  GS.Jerby[c(
    "EXHAUSTED T CELL (SPECIFIC MARKERS)",
    "CYTOTOXIC T CELL (SPECIFIC MARKERS)",
    "CYTOTOXIC CD8 T CELL",
    "NAIVE T CELL (SPECIFIC MARKERS)",
    "T FOLLICULAR HELPER",
    "TH1",
    "TH2",
    "TH9",
    "TH17",
    "TH22",
    "TREG"
  )] 

names(GS.Jerby) <-
  c(
    "Tex_specificMarkers",
    "Cytotoxic_specificMarkers",
    "Cytotoxic_CD8",
    "Naive_T_specificMarkers",
    "Tfh",
    "TH1",
    "TH2",
    "TH9",
    "TH17",
    "TH22",
    "TREG"
  )

gex.list <- AddModuleScore_UCell(gex.list, features = GS.Jerby)

featnames <- paste0(names(GS.Jerby), "_UCell")

png(paste0("Graphs/", Subset,"_Jerby_Annotation.png"), width = 2000, height = 1500 , type = "cairo")
FeaturePlot(gex.list, features = featnames, pt.size = 0.1, order = T, min.cutoff = 0.2, raster = F) #, cols = 4
dev.off()




#For Annotation with Pan-cancer Zheng et al. depending if CD4 or CD8 subset is used


if (Subset == "CD8") {
 
  Zheng_CD8 <- read_csv("data/Pan_cancer_Zheng_CD8.csv")
  load("~/BaseTIL_code/data/Zheng_Japrin/exclude.gene.misc.human.v4.RData")
  

  Zheng_high_ES <- subset(Zheng_CD8, subset =  Zheng_CD8$comb.positive.freq > 20) #Filter based on Enrichment Score
    f.feat <- !(Zheng_high_ES$geneSymbol %in% all.gene.ignore.df[["seu.id"]]) &
    !(grepl("^(RP[LS]|Rp[ls])",Zheng_high_ES$geneSymbol,perl=T)) &
    !(grepl("^MT-|^MTRN|^MTAT|^MTND|^MRP",Zheng_high_ES$geneSymbol,perl=T)) &
    !(grepl("^RP[0-9]+-",Zheng_high_ES$geneSymbol,perl=T)) &
    !(Zheng_high_ES$geneSymbol %in% c("MALAT1", 'XIST', 'NEAT1', "LINC02694", "LINC01619", "LINC02506", "AL136456.1"))
  Zheng_high_ES <- Zheng_high_ES[f.feat,]
    Zheng_high_ES %>% group_by(cluster.name) %>% top_n(n = 20, wt = comb.ES) -> Zheng_high_ES
  gene_list <- split(Zheng_high_ES, f=Zheng_high_ES$cluster.name) #Create per Cluster list object
  gene_list <- lapply(gene_list, FUN = function(X){X <- X$geneSymbol}) #Take only geneSymbols for ModuleScore
  
 
  gex.list <- AddModuleScore_UCell(gex.list, features = gene_list, name = "_Zheng")

  featnames <- colnames(gex.list@meta.data)[grepl("_Zheng$", colnames(gex.list@meta.data))]
  
  
  png(paste0("Graphs/", Subset,"_Annotation_Zheng_CD8.png"), width = 2000, height = 2000 , type = "cairo")
  print(FeaturePlot(gex.list, features = featnames, pt.size = 0.1, order = T, min.cutoff = 0.4, raster = F))
  dev.off()
  
} else {
  
  Zheng_CD4 <- read_csv("data/Pan_cancer_Zheng_CD4.csv")
  load("~/BaseTIL_code/data/Zheng_Japrin/exclude.gene.misc.human.v4.RData")
  

  Zheng_high_ES <- subset(Zheng_CD4, subset =  Zheng_CD4$comb.positive.freq > 20) #Filter based on Enrichment Score
  f.feat <- !(Zheng_high_ES$geneSymbol %in% all.gene.ignore.df[["seu.id"]]) &
    !(grepl("^(RP[LS]|Rp[ls])",Zheng_high_ES$geneSymbol,perl=T)) &
    !(grepl("^MT-|^MTRN|^MTAT|^MTND|^MRP",Zheng_high_ES$geneSymbol,perl=T)) &
    !(grepl("^RP[0-9]+-",Zheng_high_ES$geneSymbol,perl=T)) &
    !(Zheng_high_ES$geneSymbol %in% c("MALAT1", 'XIST', 'NEAT1', "LINC02694", "LINC01619", "LINC02506", "AL136456.1"))
  Zheng_high_ES <- Zheng_high_ES[f.feat,]
  Zheng_high_ES %>% group_by(cluster.name) %>% top_n(n = 20, wt = comb.ES) -> Zheng_high_ES
  gene_list <- split(Zheng_high_ES, f=Zheng_high_ES$cluster.name) #Create per Cluster list object
  gene_list <- lapply(gene_list, FUN = function(X){X <- X$geneSymbol}) #Take only geneSymbols for ModuleScore
  

  gex.list <- AddModuleScore_UCell(gex.list, features = gene_list, name = "_Zheng")

  featnames <- colnames(gex.list@meta.data)[grepl("_Zheng$", colnames(gex.list@meta.data))]
  
  png(paste0("Graphs/",Subset,"_Annotation_Zheng_CD4.png"), width = 2000, height = 2000 , type = "cairo")
  print(FeaturePlot(gex.list, features = featnames, pt.size = 0.1, order = T, min.cutoff = 0.4, raster = F))
  dev.off()
  
}



## -----------------------------------------------------------------------------------------------------------
#Find ALL Markers for Cluster determination


gex.list <- Seurat::SetIdent(gex.list, value = gex.list$seurat_clusters)

plt <- list()
allMarkers <- FindAllMarkers(gex.list,  min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE) 

allMarkers <- subset(allMarkers, allMarkers$p_val_adj<0.05)
write.csv(allMarkers, file = paste0("Graphs/",Subset,"FindAllMarkers.csv"), row.names = F, quote = F)

top_specific_markers <- allMarkers %>%
  group_by(cluster) %>%
  top_n(3, avg_log2FC)



pdf(file = paste0("Graphs/",Subset,"_dottplot_top3_Markers.pdf"), width = 15, height = 8)
DotPlot(gex.list, features = unique(top_specific_markers$gene), group.by = "seurat_clusters", cols= "RdBu")+
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 1.0,
    hjust = 1
  )
  ) + 
  xlab(NULL) + 
  ylab(NULL) 
dev.off()


allMarkers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10


png(paste0("Graphs/",Subset,"_Heatmap.png"), type = "cairo", height = 1600, width = 1000)
DoHeatmap(subset(gex.list, downsample = 300), features = top10$gene, assay = "RNA", group.by= "seurat_clusters") #+ NoLegend()
dev.off()

