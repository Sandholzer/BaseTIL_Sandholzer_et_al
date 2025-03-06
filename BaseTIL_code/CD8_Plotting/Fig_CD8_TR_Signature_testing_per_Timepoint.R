## -----------------------------------------------------------------------------------------------------------
#load dependencies

suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(UCell)
  library(ggpubr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(tidyverse)
  library(ggpubr)
})

set.seed(1234567)
setwd("~/BaseTIL_code/CD8_Plotting/")

gex.CD8 <- readRDS(file = "~/BaseTIL_code/saveFiles/CD8_annotated.rds")

## -----------------------------------------------------------------------------------------------------------
#subset only tumor-reactive T cells and remove sampels with to low numbers
TR.CD8 <- subset(gex.CD8, overall_reactive == TRUE & Patient %in% c("UPN001","UPN002", "UPN003", "UPN006", "UPN011", "UPN009", "UPN008")) #
TR.CD8 <- subset(TR.CD8, Sample_Name %in% c("UPN008 Rebiopsy1", "UPN001 Rebiopsy2"), invert = T)

TR.CD8$Type_new <- as.character(TR.CD8$Type_new)
TR.CD8$Type_post <- ifelse(TR.CD8$Type_new %in% c("post-ACT.1", "post-ACT.2"), "post-ACT", TR.CD8$Type_new)
Idents(TR.CD8) <- TR.CD8$Type_post



#aggregate based on patient and timepoint
aggro_seur <- AggregateExpression(TR.CD8, return.seurat = T, group.by = c("Patient", "Type_post"))

## -----------------------------------------------------------------------------------------------------------
#load published signatures
published_signatures <- read_csv("~/BaseTIL_code/data/Published_signatures.csv")


gene.sets <- published_signatures %>%
  pivot_longer(cols = everything(), names_to = "Column") %>%
  group_by(Column) %>%
  summarize(Values = list(value)) %>%
  deframe()

gene.sets<- lapply(gene.sets, function(x) x <- x[!is.na(x)])

selection <- names(gene.sets)
selection <- c("Mel_Exhaust_Tirosh", "Li.CD8.DYS", "Krishna.ACT.Stem.Like", "Krishna.ACT.Term.Diff")
gene.sets <- gene.sets[selection]

# caluclate signatrue score
aggro_seur <- AddModuleScore_UCell(aggro_seur, features = gene.sets, name = "")


meta.data <- data.frame(aggro_seur[[]])
meta.data$Type_post<- factor(meta.data$Type_post, levels = c("pre-ACT","Intermediate product","TIL product", "PBMC 7dpt", "post-ACT"),
                       labels = c("pre-ACT","Intermediate\nproduct","TIL product", "PBMC 7dpt", "post-ACT"))

#define significance testing
significance <- list(c("pre-ACT", "Intermediate\nproduct"),
                     c("Intermediate\nproduct", "TIL product"),
                     c("TIL product", "PBMC 7dpt"),
                     c("PBMC 7dpt", "post-ACT"))


#define colors for barplots
group_colors <- c("#E64B35B2", "#4DBBD5B2", "#00A087B2", "#3C5488B2", "#F39B7FB2")
names(group_colors)<- c("pre-ACT","Intermediate\nproduct","TIL product", "PBMC 7dpt", "post-ACT")

## -----------------------------------------------------------------------------------------------------------
# plotting for every timepoint
for (iterater in selection) {
  plot<-  ggplot(meta.data, aes(x = Type_post, y = .data[[iterater]], fill = Type_post)) +
  geom_bar(stat = "summary", fun = "mean", alpha = 0.5, width = 0.7, position = position_dodge(width = 0.8), color = NA, show.legend = FALSE) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, aes(color = Type_post), position = position_dodge(width = 0.8), show.legend = FALSE) +
  geom_jitter(aes(color = Type_post), position = position_jitterdodge(jitter.width = 1, dodge.width = 1, seed = 2), size = 1.5, show.legend = FALSE) +
  scale_fill_manual(name = "Response", values = group_colors) +
  scale_color_manual(values = group_colors) +
  theme_minimal() +
  labs(x = "", y = "Signature Score") +
  theme_classic()+
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust =1)) + 
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.15)))+ 
    stat_compare_means(
      aes(group = Type_post),
      label = "p.signif",
      method = "wilcox.test",
      comparisons = significance,
      size=7  , hide.ns = F
    )


pdf(paste0("Fig/Signature_", iterater, "_reactive_Type_bars.pdf"), height = 4, width = 4.5)
print(plot)
dev.off()
}
