## -----------------------------------------------------------------------------------------------------------
#load dependencies
suppressPackageStartupMessages({
  library(readr)
})
set.seed(1234567)
setwd("~/BaseTIL_code/TCR_bulk_sequencing")



#Metadata for UPN001_3
metadata<- as_tibble(data.frame(
  Sample = c("TIL_001_3.clones_TRB", "CL_001_3.clones_TRB"),
  Type = c("TIL", "CL"),
  Patient = c("UPN001_3", "UPN001_3"),
  stringsAsFactors = F
))

metadata$Patient_Type <- paste(metadata$Patient, metadata$Type,  sep = "_")

write_tsv(metadata,"/UPN001_3/results/metadata.txt")




#Metadata for UPN002_2
metadata<- as_tibble(data.frame(
  Sample = c("TIL_002_2.clones_TRB", "CL_002_2.clones_TRB", "FTD_002_2.clones_TRB"),
  Type = c("TIL", "CL", "FTD"),
  Patient = c("UPN002_2","UPN002_2","UPN002_2"),
  stringsAsFactors = F
))

metadata$Patient_Type <- paste(metadata$Patient, metadata$Type,  sep = "_")

write_tsv(metadata,"/UPN002_2/results/metadata.txt")




#Metadata for UPN003
metadata<- as_tibble(data.frame(
  Sample = c( "TIL_003.clones_TRB", "CD3_003.clones_TRB", "DIG_003.clones_TRB"),
  Type = c("TIL", "CD3", "DIG"),
  Patient = c("UPN003", "UPN003", "UPN003"),
  stringsAsFactors = F
))

metadata$Patient_Type <- paste(metadata$Patient, metadata$Type,  sep = "_")

write_tsv(metadata,"/UPN003/results/metadata.txt")



#Metadata for UPN006
metadata<- as_tibble(data.frame(
  Sample = c( "TIL_006.clones_TRB", "DIG_006.clones_TRB"),
  Type = c("TIL", "DIG"),
  Patient = c("UPN006", "UPN006"),
  stringsAsFactors = F
))

metadata$Patient_Type <- paste(metadata$Patient, metadata$Type,  sep = "_")

write_tsv(metadata,"/UPN006/results/metadata.txt")




#Metadata for UPN008
metadata<- as_tibble(data.frame(
  Sample = c( "TIL_008.clones_TRB", "R2FTD_008.clones_TRB", "R2DIG_008.clones_TRB", "R2CL_008.clones_TRB"),
  Type = c("TIL", "R2FTD", "R2DIG", "R2CL"),
  Patient = c("UPN008", "UPN008", "UPN008", "UPN008"),
  stringsAsFactors = F
))

metadata$Patient_Type <- paste(metadata$Patient, metadata$Type,  sep = "_")

write_tsv(metadata,"/UPN008/results/metadata.txt")



#Metadata for UPN009
metadata<- as_tibble(data.frame(
  Sample = c( "TIL_009.clones_TRB", "DIG_009.clones_TRB", "R2L_009.clones_TRB", "R2K_009.clones_TRB"),
  Type = c("TIL", "DIG", "R2L", "R2K"),
  Patient = c("UPN009", "UPN009", "UPN009", "UPN009"),
  stringsAsFactors = F
))

metadata$Patient_Type <- paste(metadata$Patient, metadata$Type,  sep = "_")

write_tsv(metadata,"/UPN009/results/metadata.txt")



#Metadata for UPN011_2
metadata<- as_tibble(data.frame(
  Sample = c( "TIL_011_2.clones_TRB", "SB_011_2.clones_TRB", "B_011_2.clones_TRB"),
  Type = c("TIL", "SB", "B"),
  Patient = c("UPN011_2", "UPN011_2", "UPN011_2"),
  stringsAsFactors = F
))

metadata$Patient_Type <- paste(metadata$Patient, metadata$Type,  sep = "_")

write_tsv(metadata,"/UPN011_2/results/metadata.txt")



