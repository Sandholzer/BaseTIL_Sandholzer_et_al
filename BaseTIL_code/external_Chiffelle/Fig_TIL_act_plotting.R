## -----------------------------------------------------------------------------------------------------------
#load dependencies

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(ggeasy)
  library(Seurat) 
  library(UCell)
  library(harmony)
  library(readr)
  library(ggeasy)
  library(ggpubr)
  library(viridis)
  library(purrr)
  library(tidyr)
  library(dplyr)
  library(stringr)
  library(ComplexHeatmap)
  library(circlize)
  
})

set.seed(12345678)

setwd("~/BaseTIL_code/external_Chiffelle")


acttil <- readRDS(file = "saveFiles/Chiffelle_TIL_act.rds")

acttil$RNR <- factor(acttil$RNR, levels = c("NR", "R"), labels = c("NRs", "Rs"))
act_aggro <- AggregateExpression(acttil, return.seurat = T, group.by = c("Sample_Name", "Response", "RNR"))
act_aggro$RNR <- factor(act_aggro$RNR, levels = c("NRs", "Rs"))

## -----------------------------------------------------------------------------------------------------------
#Plot Heatmap with interesting genes from clustering


selection <- c("RORC", "IL26", "IL17RE", "IL4I1","IL23R",  "AQP3","CCL20" ,"IL22", "IL21","LTK", "KLRB1", "ADAM12", "CCR6") #CCL4

#check if all markers are present
length(selection)
selection <- selection[selection %in% rownames(act_aggro@assays[["RNA"]]$data)]
length(selection)


mat<- GetAssayData(act_aggro, layer = "scale.data")
mat <- mat[selection,]

quantile(mat, c(0.05, 0.95))

#define order of heatmap columns
col.order <- c("patient10-ACTP_SD_NRs","patient4-ACTP_SD_NRs","patient12-ACTP_PD_NRs","patient11-ACTP_PD_NRs","patient6-ACTP_PD_NRs", "patient5-ACTP_PD_NRs", "patient1-ACTP_SD_NRs",
               "patient9-ACTP_R_Rs", "patient2-ACTP_R_Rs", "patient7-ACTP_R_Rs", "patient8-ACTP_R_Rs", "patient3-ACTP_R_Rs", "patient13-ACTP_R_Rs")
mat<- mat[ , col.order]
meta <- act_aggro@meta.data

#plot the heatmap
pdf("Fig/Fig_Heatmap_genes.pdf", height = 4, width = 6)

Heatmap(mat, name = "Z score",  
        column_split = factor(meta[col.order,]$RNR, levels = c("NRs", "Rs")),
        cluster_columns = F,
        show_column_dend = FALSE,
        cluster_column_slices = TRUE,
        column_title_gp = gpar(fontsize = 12),
        column_gap = unit(0.75, "mm"),
        cluster_rows = T,
        show_row_dend = FALSE,
        col=colorRamp2(c(-0.9, 0, 1.6), c("#3C4DA0", "white", "#CD1F2D")),
        #col=viridis(256),
        row_names_gp = gpar(fontsize = 12),
        column_title_rot = 0,
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("#CD2626", "#009ACD")))),
        show_column_names = F,
        use_raster = TRUE,
        raster_quality = 4,
        row_names_side = "left",
)



dev.off()



## -----------------------------------------------------------------------------------------------------------
# Estimate Type 17 percent in TIL products with gene selection

feat <- list()
feat$selection <- c("RORC", "IL26", "IL17RE" ,"IL23R",  "AQP3","CCL20" ,"IL22", "IL21","KLRB1", "ADAM12", "CCR6")

acttil <- AddModuleScore_UCell(acttil, features = feat, name = "")

#to define threshold for signatures score
RidgePlot(acttil, features = "selection", group.by = "Response") & geom_vline(xintercept = 0.09)

acttil@meta.data$isTh17 <- FALSE
acttil@meta.data[acttil@meta.data$selection > 0.09,"isTh17"] <- TRUE

#calculate frequency per sample
data_Th17<- dittoSeq::dittoBarPlot(acttil, var = "isTh17", group.by = "Patient", split.by = "RNR", data.out = T)
result_TIL_Th17<- data_Th17$data %>% filter(label==TRUE) %>% mutate(percent = percent*100)

# Calculate mean and standard error for each group
summary_data <- result_TIL_Th17 %>%
  group_by(RNR) %>%
  summarise(
    mean_value = mean(percent),
    sd_value = sd(percent),
    se_value = sd_value / sqrt(length(percent))
  )

bar_colors <- c("NRs" = "#CD2626", "Rs" = "#009ACD")


#plot the frequencies
pdf("Fig/TIL_Type17_total_percent_response.pdf", height = 3, width = 2)
ggplot(result_TIL_Th17, aes(x = RNR, y = percent, fill = RNR)) +
  geom_bar(stat = "summary", fun = "mean", alpha = 0.4, width = 0.7, position = position_dodge(width = 0.8), color = NA, show.legend = FALSE) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, aes(color = RNR), position = position_dodge(width = 0.8), show.legend = FALSE) +
  geom_jitter(aes(color = RNR), position = position_jitterdodge(jitter.width = 1, dodge.width = 0.8, seed =2), size = 1.5, show.legend = FALSE) +
  scale_fill_manual(name = "Response", values = bar_colors) +
  scale_color_manual(values = bar_colors) +
  theme_minimal() +
  labs(x = "", y = "Percent of TIL product") +
  ggtitle("Type 17")+
  theme_classic()+
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = 16),
    axis.line = element_line(color = "black"),
    plot.title = element_text(size = 15,hjust = 0.5, face = "bold")
  ) + geom_signif(comparisons = list(c("Rs", "NRs")),
                  map_signif_level = TRUE, textsize = 12, vjust = 0.5, test = "wilcox.test") +
  ylim(c(0,14))

dev.off()



## -----------------------------------------------------------------------------------------------------------
# Estimate Type 17 percent in pre-ACT with gene selection

prepost <- readRDS("~/BaseTIL/Chiffelle/saveFiles/prePost_Tcells_clustered.rds")
prepost$RNR <- factor(prepost$RNR, levels = c("NR", "R"), labels = c("NRs", "Rs"))

pre_act <- subset(prepost, prepost == "pre")
pre_act <- subset(pre_act, seurat_clusters != 9) #remove NK cells



feat <- list()
feat$selection <- c("RORC", "IL26", "IL17RE" ,"IL23R",  "AQP3","CCL20" ,"IL22", "IL21"
                    ,"KLRB1", "ADAM12", "CCR6") #


pre_act <- AddModuleScore_UCell(pre_act, features = feat, name = "")

RidgePlot(pre_act, features = "selection", group.by = "Response") & geom_vline(xintercept = 0.15)

pre_act@meta.data$isTh17 <- FALSE
pre_act@meta.data[pre_act@meta.data$selection >= 0.15,"isTh17"] <- TRUE

dittoSeq::dittoBarPlot(pre_act, var = "isTh17", group.by = "Patient", split.by = "Response")

data_Th17<- dittoSeq::dittoBarPlot(pre_act, var = "isTh17", group.by = "Patient", split.by = "RNR", data.out = T)

result_tumor_Th17<- data_Th17$data %>% filter(label==TRUE) %>% mutate(percent = percent*100)

# Calculate mean and standard error for each group
summary_data <- result_tumor_Th17 %>%
  group_by(RNR) %>%
  summarise(
    mean_value = mean(percent),
    sd_value = sd(percent),
    se_value = sd_value / sqrt(length(percent))
  )

bar_colors <- c("NRs" = "#CD2626", "Rs" = "#009ACD")



pdf("Fig/Tumor_Type17_total_percent_response.pdf", height = 3, width = 2)

ggplot(result_tumor_Th17, aes(x = RNR, y = percent, fill = RNR)) +
  geom_bar(stat = "summary", fun = "mean", alpha = 0.4, width = 0.7, position = position_dodge(width = 0.8), color = NA, show.legend = FALSE) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, aes(color = RNR), position = position_dodge(width = 0.8), show.legend = FALSE) +
  geom_jitter(aes(color = RNR), position = position_jitterdodge(jitter.width = 1, dodge.width = 0.8), size = 1.5, show.legend = FALSE) +
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
  ) + geom_signif(comparisons = list(c("Rs", "NRs")),
                  map_signif_level = TRUE, textsize = 12, vjust = 0.5, test = "wilcox.test") +
  ylim(c(0,7))

dev.off()






## -----------------------------------------------------------------------------------------------------------
# Estimate Type 17 changes during TIL expansion


result_TIL_Th17$Type <- "TIL"
result_tumor_Th17$Type <- "pre-ACT"

results_merged <- merge(result_TIL_Th17, result_tumor_Th17, by=c("grouping", "label", "RNR"), suffix= c(".TIL", ".Tum"))
results_merged <- results_merged %>% select(grouping, percent.TIL, percent.Tum, RNR )

longer_results<- pivot_longer(results_merged, cols = c("percent.TIL", "percent.Tum"))
longer_results$Type <- ifelse(longer_results$name == "percent.TIL", "TIL product", "pre-ACT") 


pdf("Fig/Type17_freq_expansion.pdf", height = 4, width = 4)
longer_results %>% ggplot(aes(
  factor(
    Type,
    levels = c(
      "pre-ACT",
      "TIL product"
    )
  ),
  value,
  color = RNR,
  group = grouping
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

