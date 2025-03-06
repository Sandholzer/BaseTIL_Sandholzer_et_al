#!/bin/bash

#SBATCH --job-name=FullAnalysis_BaseTIL                #This is the name of your job
#SBATCH --cpus-per-task=15                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=16G              #This is the memory reserved per core.
#Total memory reserved: 240GB

# Are you sure that you need THAT much memory?

#SBATCH --time=24:00:00        #This is the time that your task will run
#SBATCH --qos=1day           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=myrun.oe%j     #This is the joined STDOUT and STDERR file
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=michael.sandholzer@unibas.ch        #You will be notified via email when your task ends or fails

#This job runs from the current working directory


#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $SLURM_JOBID stores the ID number of your job.


#load your required modules below
#################################
ml purge
module load R/4.4.1-foss-2023b

#export your required environment variables below
#################################################
export OPENBLAS_NUM_THREADS=1


#add your command lines below
#############################


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
