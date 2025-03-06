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
setwd("~/BaseTIL_code/CD8_Plotting")

gex.CD8 <- readRDS(file = "../saveFiles/CD8_annotated.rds")

#remove control PBMC sample from analysis
gex.CD8 <- subset(gex.CD8, Type != "PBMC_Ctrl")

gex.CD8$Type_temp <- factor(gex.CD8$prepost,
                            levels = rev(c("pre-ACT","Intermediate product","TIL product", "PBMC 7dpt", "post-ACT")))


## -----------------------------------------------------------------------------------------------------------
# plotting distribution per sample Type

pdf("Fig/CD8_Type_distribution.pdf", height = 5, width = 5)
dittoBarPlot(
  gex.CD8,
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
  ) + scale_x_discrete(limits= rev(c("Tn", "Tm", "Tem", "Teff", "Temra", "Tex", "Tc17", "FOXP3", "Tprol1", "Tprol2", "Tprol3")))+
  ggtitle(NULL)
  dev.off()


  
  
  
  
  
  
  
## -----------------------------------------------------------------------------------------------------------
# plotting Frequency per cluster
  
  table <- dittoBarPlot(gex.CD8, var = "Celltype", group.by = "Sample_Name_new", retain.factor.levels = F,data.out=T)
  cell_freq <- table$data

  cell_freq$percent <- cell_freq$percent*100
  cell_freq[c('Patient', 'Type')] <- str_split_fixed(cell_freq$grouping, ' ', 2)
  
  
  cell_freq[cell_freq$Type %in% c("post-ACT.1", "post-ACT.2"),"Type"] <- "post-ACT"
  
  
  plt <- list()
  

  Colors <- RColorBrewer::brewer.pal(n=12, "Paired")
  names(Colors) <- c("Tn", "Tm", "Tem", "Teff", "Temra", "Tex", "Tc17", "FOXP3", "Tprol1", "Tprol2", "Tprol3")
  
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
      scale_x_discrete(limits= c("Tn", "Tm", "Tem", "Teff", "Temra", "Tex", "Tc17", "FOXP3", "Tprol1", "Tprol2", "Tprol3"))+
      scale_fill_manual(values=Colors) + 
      NoLegend()
    
  }
  

  
  pdf(paste0("Fig/CD8_Cluster_Frequency_perType.pdf"), height = 3, width = 15)
  ggarrange(plotlist =  c(plt["pre-ACT"], plt["Intermediate product"],plt["TIL product"] ,plt["PBMC 7dpt"], plt["post-ACT"]), 
            nrow = 1, common.legend = T, legend = "none")
   dev.off()
  
  
  
  
  
## -----------------------------------------------------------------------------------------------------------
# plotting Frequency per cluster of tumor reactive
  
  gex.sub <- subset(gex.CD8, overall_reactive == TRUE)
  
  table <- dittoBarPlot(gex.sub, var = "Celltype", group.by = "Sample_Name_new", retain.factor.levels = F,data.out=T)
  
  
  
  cell_freq <- table$data
  # test$Response <- NULL
  
  cell_freq$percent <- cell_freq$percent*100
  cell_freq[c('Patient', 'Type')] <- str_split_fixed(cell_freq$grouping, ' ', 2)
  
  #cell_freq <- cell_freq[cell_freq$Patient != "UPN006",]
  cell_freq <- cell_freq[cell_freq$grouping != "UPN001 post-ACT.2",]
  cell_freq <- cell_freq[cell_freq$grouping != "UPN008 post-ACT.1",]
  
  cell_freq[cell_freq$Type %in% c("post-ACT.1", "post-ACT.2"),"Type"] <- "post-ACT"
  
  
  plt <- list()
  
  
  
  library(RColorBrewer)
  Colors <- RColorBrewer::brewer.pal(n=12, "Paired")
  names(Colors) <- c("Tn", "Tm", "Tem", "Teff", "Temra", "Tex", "Tc17", "FOXP3", "Tprol1", "Tprol2", "Tprol3")
  
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
      ylim(c(0, 80)) +
      ggtitle(Type_iterater)+
      labs( #paste0(Type_iterater)
        x = "",
        y = "Percent") +
      theme_classic() +
      theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1),
            axis.text.y = element_text(size=12),
            axis.title.y = element_text(size = 12),
            plot.title = element_text(size=14, face = "bold", hjust = 0.5)) + 
      scale_x_discrete(limits= c("Tn", "Tm", "Tem", "Teff", "Temra", "Tex", "Tc17", "FOXP3", "Tprol1", "Tprol2", "Tprol3"))+
      scale_fill_manual(values=Colors) + 
      NoLegend()
    
  }
  

  
  pdf(paste0("Fig/CD8_Overall_TR_Cluster_Frequency_perType.pdf"), height = 3, width = 15)
  ggarrange(plotlist =  c(plt["pre-ACT"], plt["Intermediate product"],plt["TIL product"] ,plt["PBMC 7dpt"], plt["post-ACT"]), 
            nrow = 1, common.legend = T, legend = "none")

  dev.off()
  
  
  
  
