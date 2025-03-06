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
  library(RColorBrewer)
  library(dittoSeq)
})

set.seed(1234567)
setwd("~/BaseTIL_code/CD4_Plotting")

gex.CD4 <- readRDS(file = "../saveFiles/CD4_annotated.rds")

#remove control PBMC sample from analysis
gex.CD4 <- subset(gex.CD4, Type != "PBMC_Ctrl")

gex.CD4$Type_temp <- factor(gex.CD4$prepost,
                            levels = rev(c("pre-ACT","Intermediate product","TIL product", "PBMC 7dpt", "post-ACT")))

## -----------------------------------------------------------------------------------------------------------
# plotting distribution per sample Type

pdf("Fig/CD4_Type_distribution.pdf", height = 5, width = 5)
dittoBarPlot(
  gex.CD4,
  group.by = "Celltype",
  var = "Type_temp",
  retain.factor.levels = T) +
  scale_fill_manual(
    values = c(
      "pre-ACT" = "#E64B35B2",
      "Intermediate product" = "#4DBBD5B2",
      "TIL product" =  "#00A087B2",
      "PBMC 7dpt" = "#3C5488B2",
      "post-ACT" =  "#F39B7FB2"
    ), 
    breaks = c("pre-ACT", "Intermediate product", "TIL product", "PBMC 7dpt", "post-ACT")) + 
  coord_flip() +
  xlab(NULL) + theme(
    legend.position = "right", legend.text = element_text(size=15),
    legend.title = element_blank(),
    axis.text.x = element_text(
      size=15,
      angle = 0,
      vjust = 0.5,
      hjust = 0.5
    ),
    axis.title=element_text(size=15),
    axis.text.y = element_text(size=15)
  ) + scale_x_discrete(limits= rev(c("Tn", "Tm", "Teff", "Th1/Tex", "Tfh", "Th2", "Treg", "Th17", "Tprol")))+
  ggtitle(NULL)
dev.off()









## -----------------------------------------------------------------------------------------------------------
# plotting Frequency per cluster


table <- dittoBarPlot(gex.CD4, var = "Celltype", group.by = "Sample_Name_new", retain.factor.levels = F,data.out=T)
cell_freq <- table$data


cell_freq$percent <- cell_freq$percent*100
cell_freq[c('Patient', 'Type')] <- str_split_fixed(cell_freq$grouping, ' ', 2)


cell_freq[cell_freq$Type %in% c("post-ACT.1", "post-ACT.2"),"Type"] <- "post-ACT"


plt <- list()




Colors <- RColorBrewer::brewer.pal(n=12, "Paired")


names(Colors) <- c("Tn", "Tm", "Teff", "Th1/Tex", "Tfh", "Th2", "Treg", "Th17", "Tprol")

for (Type_iterater in unique(cell_freq$Type)) {
  
  cell_freq.sub <- cell_freq[cell_freq$Type == Type_iterater,]
  
  summary_data <- cell_freq.sub %>%
    group_by(label) %>%
    summarise(mean_percentage = mean(percent),
              sd_percentage = sd(percent),
              n = n())
  
  
  plt[[Type_iterater]] <- ggplot(summary_data, aes(x = label, y = mean_percentage, fill = label)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = mean_percentage - sd_percentage/sqrt(n), 
                      ymax = mean_percentage + sd_percentage/sqrt(n)),
                  width = 0.2, position = position_dodge(width = 0.9)) +
    # geom_point(data = cell_freq.sub, aes(x = label, y = percent, color = factor(Patient)), 
    #            position = position_dodge(width = 0.5), 
    #            size = 1.5, alpha = 0.7) + 
    ylim(c(0, 60)) +
    ggtitle(Type_iterater)+
    labs( #paste0(Type_iterater)
      x = "",
      y = "Percent") +
    theme_classic() +
    theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size = 12),
          plot.title = element_text(size=14, face = "bold", hjust = 0.5)) + 
    scale_x_discrete(limits= c("Tn", "Tm", "Teff", "Th1/Tex", "Tfh", "Th2", "Treg", "Th17", "Tprol"))+
    scale_fill_manual(values=Colors) + 
    NoLegend()
  
}






pdf(paste0("Fig/CD4_Cluster_Frequency_perType.pdf"), height = 3, width = 15)
ggarrange(plotlist =  c(plt["pre-ACT"], plt["Intermediate product"],plt["TIL product"] ,plt["PBMC 7dpt"], plt["post-ACT"]), 
          nrow = 1, common.legend = T, legend = "none")

dev.off()




## -----------------------------------------------------------------------------------------------------------
# plotting Treg frequencies before and after transfer


data<- dittoSeq::dittoBarPlot(subset(gex.CD4, Type_new %in% c("TIL product", "PBMC 7dpt")), var = "Celltype", group.by = "Sample_Name_new", split.by = "Response", data.out = T)
data <- data$data
data[c('Patient', 'Type_new')] <- str_split_fixed(data$grouping, ' ', 2)

result<- data %>% filter(Type_new %in% c("TIL product", "PBMC 7dpt"), label == "Treg") %>% mutate(percent = percent*100)
result$RNR <- ifelse(result$Response %in% c("PD", "SD"), "NRs", "Rs")

result$RNR <- factor(result$RNR, levels = c("NRs", "Rs"))
result$Type_new <- factor(result$Type_new, levels = c("TIL product", "PBMC 7dpt"))


# Define colors
bar_colors <- c("NRs" = "#CD2626", "Rs" = "#009ACD")


# Create the plot
pdf("Fig/CD4_Transfer_Treg_percent_response.pdf", height = 3, width = 4)
# Convert group to factor for better plotting

# Create the bar plot
  p <- ggplot(result, aes(x = Type_new, y = percent, fill = RNR)) +
    geom_bar(stat = "summary", fun = "mean", alpha = 0.4, width = 0.7, position = position_dodge(width = 0.8), color = NA) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, aes(color = RNR), position = position_dodge(width = 0.8), show.legend = FALSE) +
    geom_jitter(aes(color = RNR), position = position_jitterdodge(jitter.width = 0.35, dodge.width = 0.8, seed = 1), size = 1.5, show.legend = FALSE) +
    scale_fill_manual(name = "Response", values = bar_colors) +
    scale_color_manual(values = bar_colors) +
    theme_minimal() +
    labs(x = "", y = "Percent of CD4 T cells") +
    ggtitle("Tregs during transfer")+
    theme_classic()+
    theme(
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.text = element_text(size = 10),
      # panel.grid.major = element_blank(),
      # panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    ) 
  
  p
  p+stat_compare_means(aes(group = RNR), label = "p.signif", method = "t.test")
dev.off()






## -----------------------------------------------------------------------------------------------------------
# plotting Treg frequencies before and after therapy


data<- dittoSeq::dittoBarPlot(subset(gex.CD4, Type_new %in% c("post-ACT.1", "post-ACT.2", "pre-ACT")), var = "Celltype", group.by = "Sample_Name_new", split.by = "Response", data.out = T)
data <- data$data
data[c('Patient', 'Type_new')] <- str_split_fixed(data$grouping, ' ', 2)

result<- data %>% filter(Type_new %in% c("post-ACT.1", "post-ACT.2", "pre-ACT"), label == "Treg") %>% mutate(percent = percent*100)
result$prepost <- ifelse(result$Type_new == "pre-ACT", "pre-ACT", "post-ACT")
result$RNR <- ifelse(result$Response %in% c("PD", "SD"), "NRs", "Rs")
result$prepost <- factor(  result$prepost, level= c("pre-ACT", "post-ACT"))
result$RNR <- factor(result$RNR, levels = c("NRs", "Rs"))


# Define colors
bar_colors <- c("NRs" = "#CD2626", "Rs" = "#009ACD")



# Create the plot
pdf("Fig/CD4_Tumor_Treg_percent_response.pdf", height = 3, width = 4)

p <- ggplot(result, aes(x = prepost, y = percent, fill = RNR)) +
  geom_bar(stat = "summary", fun = "mean", alpha = 0.4, width = 0.7, position = position_dodge(width = 0.8), color = NA) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, aes(color = RNR), position = position_dodge(width = 0.8), show.legend = FALSE) +
  geom_jitter(aes(color = RNR), position = position_jitterdodge(jitter.width = 0.35, dodge.width = 0.8)
              , size = 1.5, show.legend = FALSE) +
  scale_fill_manual(name = "Response", values = bar_colors) +
  scale_color_manual(values = bar_colors) +
  theme_minimal() +
  labs(x = "", y = "Percent of CD4 T cells") +
  ggtitle("Tregs in lesions")+
  theme_classic()+
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 10),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) 

p
p+stat_compare_means(aes(group = RNR), label = "p.signif", method = "t.test")

dev.off()















