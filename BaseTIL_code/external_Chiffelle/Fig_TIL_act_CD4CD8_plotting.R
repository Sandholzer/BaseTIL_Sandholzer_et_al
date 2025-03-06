## -----------------------------------------------------------------------------------------------------------
#load dependencies

suppressPackageStartupMessages({
  library(Seurat)
  library(readr)
  library(dplyr)
  library(UCell)
  library(ggplot2)
  library(ggpubr)
  library(ComplexHeatmap)
  library(circlize)
  
})

set.seed(12345678)

setwd("~/BaseTIL_code/external_Chiffelle")


## -----------------------------------------------------------------------------------------------------------
#load datasets

acttil.CD8 <- readRDS("~/BaseTIL/Chiffelle/saveFiles/acttil_CD8_clustered.rds")
acttil.CD8$RNR <- factor(acttil.CD8$RNR, levels = c("NR", "R"), labels = c("NRs", "Rs"))

acttil.CD4 <- readRDS("~/BaseTIL/Chiffelle/saveFiles/acttil_CD4_clustered.rds")
acttil.CD4$RNR <- factor(acttil.CD4$RNR, levels = c("NR", "R"), labels = c("NRs", "Rs"))


## -----------------------------------------------------------------------------------------------------------
#load Type 17 specific gene lists


feat <- list()
markers_CD8<- read_csv("~/BaseTIL_code/CD8_Plotting/Fig/DGE_results/CD8_TIL_prod_cluster.csv")
feat$Tc17_1p <- markers_CD8%>% filter(cluster == 8,p_val_adj < 1e-300, avg_log2FC >1 ) %>% pull(gene)

markers_CD4<- read_csv("~/BaseTIL_code/CD4_Plotting/Fig/DGE_results/CD4_TIL_prod_cluster.csv")
feat$Th17_1p <- markers_CD4%>% filter(cluster == 2,p_val_adj < 1e-300, avg_log2FC >1 ) %>% pull(gene)


acttil.CD8 <- AddModuleScore_UCell(acttil.CD8, features = feat, name = "")


acttil.CD4 <- AddModuleScore_UCell(acttil.CD4, features = feat, name = "")

## -----------------------------------------------------------------------------------------------------------
#define cells with high Th17 signature score in CD4


acttil.CD4@meta.data$isTh17 <- FALSE
acttil.CD4@meta.data[acttil.CD4@meta.data$Th17_1p > 0.2,"isTh17"] <- TRUE

data_Th17<- dittoSeq::dittoBarPlot(acttil.CD4, var = "isTh17", group.by = "Patient", split.by = "RNR", data.out = T)

result_Th17<- data_Th17$data %>% filter(label==TRUE) %>% mutate(percent = percent*100)

# Calculate mean and standard error for each group
summary_data <- result_Th17 %>%
  group_by(RNR) %>%
  summarise(
    mean_value = mean(percent),
    sd_value = sd(percent),
    se_value = sd_value / sqrt(length(percent))
  )

bar_colors <- c("NRs" = "#CD2626", "Rs" = "#009ACD")

# Create the plot
pdf("Fig/CD4_ACT_Th17_percent_response.pdf", height = 3, width = 2)

ggplot(result_Th17, aes(x = RNR, y = percent, fill = RNR)) +
  geom_bar(stat = "summary", fun = "mean", alpha = 0.4, width = 0.7, position = position_dodge(width = 0.8), color = NA, show.legend = FALSE) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, aes(color = RNR), position = position_dodge(width = 0.8), show.legend = FALSE) +
  geom_jitter(aes(color = RNR), position = position_jitterdodge(jitter.width = 1, dodge.width = 0.8), size = 1.5, show.legend = FALSE) +
  scale_fill_manual(name = "Response", values = bar_colors) +
  scale_color_manual(values = bar_colors) +
  theme_minimal() +
  labs(x = "", y = "Percent of CD4 T cells") +
  ggtitle("Th17")+
  theme_classic()+
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = 16),
    axis.line = element_line(color = "black"),
    plot.title = element_text(size = 15,hjust = 0.5, face = "bold")
  ) + geom_signif(comparisons = list(c("Rs", "NRs")),
                  map_signif_level = TRUE, textsize = 12, vjust = 0.5, test = "wilcox.test") +
  scale_y_continuous(limits = c(0, 35))

dev.off()






## -----------------------------------------------------------------------------------------------------------
#define cells with high Tc17 signature score in CD8

acttil.CD8@meta.data$isTc17 <- FALSE
acttil.CD8@meta.data[acttil.CD8@meta.data$Tc17_1p > 0.085,"isTc17"] <- TRUE

data_Tc17<- dittoSeq::dittoBarPlot(acttil.CD8, var = "isTc17", group.by = "Patient", split.by = "RNR", data.out = T)

result_Tc17<- data_Tc17$data %>% filter(label==TRUE) %>% mutate(percent = percent*100)

# Calculate mean and standard error for each group
summary_data <- result_Tc17 %>%
  group_by(RNR) %>%
  summarise(
    mean_value = mean(percent),
    sd_value = sd(percent),
    se_value = sd_value / sqrt(length(percent))
  )

bar_colors <- c("NRs" = "#CD2626", "Rs" = "#009ACD")

# Create the plot
pdf("Fig/CD8_ACT_Tc17_percent_response.pdf", height = 3, width = 2)

ggplot(result_Tc17, aes(x = RNR, y = percent, fill = RNR)) +
  geom_bar(stat = "summary", fun = "mean", alpha = 0.4, width = 0.7, position = position_dodge(width = 0.8), color = NA, show.legend = FALSE) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, aes(color = RNR), position = position_dodge(width = 0.8), show.legend = FALSE) +
  geom_jitter(aes(color = RNR), position = position_jitterdodge(jitter.width = 1, dodge.width = 0.8, seed=4), size = 1.5, show.legend = FALSE) +
  scale_fill_manual(name = "Response", values = bar_colors) +
  scale_color_manual(values = bar_colors) +
  theme_minimal() +
  labs(x = "", y = "Percent of CD8 T cells") +
  ggtitle("Tc17")+
  theme_classic()+
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = 16),
    axis.line = element_line(color = "black"),
    plot.title = element_text(size = 16,hjust = 0.5, face = "bold")
  ) + geom_signif(comparisons = list(c("Rs", "NRs")),
                  map_signif_level = TRUE, textsize = 12, vjust = 0.5, test = "wilcox.test") +
  ylim(c(0,9))


dev.off()



