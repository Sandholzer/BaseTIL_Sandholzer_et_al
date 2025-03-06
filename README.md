## Longitudinal profiling of tumor-reactive T cells during TIL therapy reveals resistance linked to co-transfer of Type 17 T cells

Authors: 
Michael T. Sandholzer<1*></sup>, Alessia G. Liner1†, Clara Serger1†, Sarp Uzun3, David König1,2,8, Helen Thut1, Reto Ritschard1, Andreas Zingg1, Natalia Rodrigues Mantuano1, Benjamin Kasenda2, Katharina Glatz3, Elisabeth Kappos4, Matthias Matter3, Andreas Holbro6,8, Jakob Passweg6,8, Nina Khanna7,8, Lukas Jeker1,8, Mascha Binder1,2,8, Alfred Zippelius1,2,8, and Heinz Läubli1,2,8*
Affiliations:
1Department of Biomedicine, University of Basel and University Hospital Basel, Switzerland.
2Division of Medical Oncology, University Hospital Basel, Switzerland.
3Department of Pathology, University Hospital Basel, Switzerland.
4Departement of Plastic, Reconstructive, Aesthetic and Hand surgery, University Hospital Basel and University of Basel 
5Division of Dermatology, University Hospital Basel, Switzerland.
6Division of Hematology, University Hospital Basel, Switzerland.
7Division of Infectious Diseases, University Hospital Basel, Switzerland.
8Innovation Focus Cell Therapies, University Hospital Basel, Switzerland

Abstract 
Adoptive cell therapy (ACT) using ex vivo-expanded autologous tumor-infiltrating lymphocytes (TILs) yields durable responses in metastatic melanoma, yet a significant proportion of patients fail to achieve therapeutic benefit, and the mechanisms driving treatment resistance remain poorly understood. Here, we provide a comprehensive analysis of tumor-reactive T cell dynamics during TIL therapy using single-cell RNA and T cell receptor (TCR) sequencing from seven patients in the phase I BaseTIL trial. Tumor-reactive T cells expanded early during ex vivo TIL production, transitioning from exhausted to effector or memory T cell states. Post-transfer, these cells formed a memory-like reservoir in peripheral blood but reacquired exhaustion markers upon infiltrating metastatic lesions, particularly in non-responders. Responders retained a pool of less differentiated tumor-reactive TILs with stem-like characteristics in tumor lesions post-ACT. In contrast, non-responders exhibited increased proportions of regulatory T cells (Tregs) in both circulation and metastatic lesions, accompanied by an enrichment of Type 17 bystander T cells (Th17 and Tc17) in their TIL infusion products. Using a melanoma mouse model, we demonstrated that the transfer of Th17 cells drives Treg accumulation and significantly impairs tumor control. Our findings offer novel insights into clonal dynamics of TIL therapy and identify Th17-driven Treg expansion as a critical resistance mechanism, highlighting actionable strategies, such as optimizing TIL expansion protocols and targeting Type 17-driven pathways, to enhance ACT efficacy.


## Goal of the github
This github contains the scripts to recapitulate the anaylsis performed for this article. 


## Run cellranger multi
sbatch ~/BaseTIL_code/Cellranger_multi/run_cellranger/run_cellranger.sh

## Run velocyto
sbatch ~/BaseTIL_code/Cellranger_multi/velocyto/rename_bam_bcmatrix.sh
sbatch ~/BaseTIL_code/Cellranger_multi/velocyto/run_scVelo_loop.sh

## Run mixcr for alligning bulk TCRseq data
sbatch ./TCR_bulk_sequencing/run_mixcr.sh
Rscript ./TCR_bulk_sequencing/make_metadata.R #creates meta data table to load bulk TCRseq data into immunarch


## T cell clustering and splitting into CD4 and CD8
Rscript ./BaseTIL_00_createRAW.R
Rscript ./BaseTIL_01_QC.R
Rscript ./BaseTIL_02_integration_Tcells.R
Rscript ./BaseTIL_03a_Subsetting_CD4_CD8.R
Rscript ./BaseTIL_03b_Subsetting_purification_on_cluster.R
Rscript ./BaseTIL_04_Subsetting_purification_on_genes.R


## Go on with CD8 subset analysis
Rscript ./BaseTIL_05_PlottingForAnnotation_CDX.R CD8
Rscript ./CD8_Plotting/BaseTIL_06_Annotation_CD8.R
Rscript ./CD8_Plotting/Fig_CD8_Zheng_annotation_heatmap.R
Rscript ./CD8_Plotting/Fig_CD8_frequency_plotting.R

Rscript ./CD8_Plotting/Fig_CD8_DGE_Tumor_reactive_signature.R #Creates Signature used in Expansion_Contraction
Rscript ./CD8_Plotting/Fig_CD8_Expansion_Contraction.R #saves tumor reactive object for further analysis
Rscript ./CD8_Plotting/Fig_CD8_TR_TCR_on_UMAP.R 
Rscript ./CD8_Plotting/Fig_CD8_DGE_reactive_in_Timepoint.R

Rscript ./Scenic/save_expression_asLoom.R # saves the filtered gex data for CD8 as loom files per timepoint

##Run scenic per timepoint (full run exeeds processing power)
sbatch ./Scenic/run_pyScenic_demand.sh Tumor
sbatch ./Scenic/run_pyScenic_demand.sh PreREP
sbatch ./Scenic/run_pyScenic_demand.sh Expanded_TILs
sbatch ./Scenic/run_pyScenic_demand.sh PBMC_7dpt
sbatch ./Scenic/run_pyScenic_demand.sh Rebiopsy1
sbatch ./Scenic/run_pyScenic_demand.sh Rebiopsy2

Rscript ./CD8_Plotting/Fig_CD8_scenic_DGE_reactive_in_Timepoint.R # make sure scenic pipeline was executed first

#Plotting on CD8 tumor reactive subset
Rscript ./CD8_Plotting/Fig_CD8_TR_Heatmap_genes_between_Timpoint.R
Rscript ./CD8_Plotting/Fig_CD8_TR_Signature_testing_per_Timepoint.R


#Reclustering and plotting of CD8 tumor-reactive subset
Rscript ./CD8_Plotting/Fig_CD8_TR_reclustering.R
Rscript ./CD8_Plotting/Fig_CD8_TR_reclustering_plotting.R
Rscript ./CD8_Plotting/Fig_CD8_TR_Pseudotime.R
sbatch ./CD8_Plotting/run_velocity.sh
Rscript ./CD8_Plotting/Fig_CD8_TCR_Overlap.R
Rscript ./CD8_Plotting/Fig_CD8_ExpTIL_Tc17.R  # TIL product subclustering CD8


## Go on with CD4 subset analysis
Rscript ./BaseTIL_05_PlottingForAnnotation_CDX.R CD4
Rscript ./CD4_Plotting/BaseTIL_06_Annotation_CD4.R
Rscript ./CD4_Plotting/Fig_CD4_Zheng_annotation_heatmap.R
Rscript ./CD4_Plotting/Fig_CD4_frequency_plotting.R
Rscript ./CD4_Plotting/Fig_CD4_Expansion_Contraction.R
Rscript ./CD4_Plotting/Fig_CD4_TIL_Bulk_overlap.R
Rscript ./CD4_Plotting/Fig_CD4_Treg_infiltration_postACT.R
Rscript ./CD4_Plotting/Fig_CD4_DGE_Treg_PrevsPost.R
Rscript ./CD4_Plotting/Fig_CD4_ExpTIL_Th17.R  # TIL product subclustering CD4
Rscript ./Fig_Type17_freq_exp.R # Calculate total Type17 frequencies 

## TCR plotting 
Rscript ./Fig_Frequency_Tumorreactive.R
Rscript ./Fig_TCR_richness.R
Rscript ./Fig_TCR_Overlap.R


## Chiffelle
Rscript ./external_Chiffelle/TIL_act_clustering.R
Rscript ./external_Chiffelle/TIL_act_CD4_CD8_clustering.R
Rscript ./external_Chiffelle/PrePost_clustering.R

Rscript ./external_Chiffelle/Fig_TIL_DGE_Response.R 
Rscript ./external_Chiffelle/Fig_TIL_act_plotting.R
Rscript ./external_Chiffelle/Fig_TIL_act_CD4CD8_plotting.R

## Thompson
Rscript ./external_ThompsonThompson_Th17.R
