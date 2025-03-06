## -----------------------------------------------------------------------------------------------------------
#load dependencies

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
})

set.seed(1234567)
setwd("~/BaseTIL_code/Scenic")

gex.CD8 <- readRDS(file = "../saveFiles/CD8_annotated.rds")

## -----------------------------------------------------------------------------------------------------------
#save full loom
SaveLoom(gex.CD8, "~/BaseTIL_code/Scenic/pySenic_data/expression_matrix/CD8_full.loom", overwrite = T)


## -----------------------------------------------------------------------------------------------------------
#save loom per timepoint

for (type_vec in unique(gex.CD8$Type)) {
  seuratobject <- subset(gex.CD8, Type == type_vec)

  SaveLoom(seuratobject, paste0("~/BaseTIL_code/Scenic/pySenic_data/expression_matrix/CD8_",type_vec, ".loom"), overwrite = T)
}
