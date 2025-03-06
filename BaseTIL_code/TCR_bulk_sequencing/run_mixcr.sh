#!/bin/bash

#SBATCH --job-name=mixcr                   #This is the name of your job
#SBATCH --cpus-per-task=8                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=8G              #This is the memory reserved per core.
#Total memory reserved: 64GB

#SBATCH --time=6:00:00        #This is the time that your task will run
#SBATCH --qos=6hours           #You will run in this queue

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
conda activate mixcr

# mixcr needs to be connected to the internet onece for licence so start at terminal first
# mixcr analyze thermofisher-human-rna-trb-oncomine-sr --help

#export your required environment variables below
#################################################


#add your command lines below
#############################

#UPN001_3
mixcr analyze thermofisher-human-rna-trb-oncomine-sr \
      ~/BaseTIL_code/TCR_bulk_sequencing/UPN001_3/CL_001_03_TCR_TIL_39_IonDual_A11_0181.fastq \
      ~/BaseTIL_code/TCR_bulk_sequencing/UPN001_3/results/CL_001_3

mixcr analyze thermofisher-human-rna-trb-oncomine-sr \
     ~/BaseTIL_code/TCR_bulk_sequencing/UPN001_3/TIL_001_3_TCR_TIL_38_IonDual_H12_0196.fastq \
     ~/BaseTIL_code/TCR_bulk_sequencing/UPN001_3/results/TIL_001_3


#UPN002_2
mixcr analyze thermofisher-human-rna-trb-oncomine-sr \
      ~/BaseTIL_code/TCR_bulk_sequencing/UPN002_2/TIL_002_2_TCR_TIL_35_IonDual_E12_0193.fastq \
      ~/BaseTIL_code/TCR_bulk_sequencing/UPN002_2/results/TIL_002_2

mixcr analyze thermofisher-human-rna-trb-oncomine-sr \
     ~/BaseTIL_code/TCR_bulk_sequencing/UPN002_2/FTD_002_2_TCR_TIL_36_IonDual_F12_0194.fastq \
     ~/BaseTIL_code/TCR_bulk_sequencing/UPN002_2/results/FTD_002_2

mixcr analyze thermofisher-human-rna-trb-oncomine-sr \
     ~/BaseTIL_code/TCR_bulk_sequencing/UPN002_2/CL_002_2_TCR_TIL_37_IonDual_G12_0195.fastq \
     ~/BaseTIL_code/TCR_bulk_sequencing/UPN002_2/results/CL_002_2
  
     
#UPN003
mixcr analyze thermofisher-human-rna-trb-oncomine-sr \
     ~/BaseTIL_code/TCR_bulk_sequencing/UPN003/TIL_003_TCR_TIL_9_IonDual_B9_0166.fastq \
     ~/BaseTIL_code/TCR_bulk_sequencing/UPN003/results/TIL_003

mixcr analyze thermofisher-human-rna-trb-oncomine-sr \
     ~/BaseTIL_code/TCR_bulk_sequencing/UPN003/CD3_003_TCR_TIL_10_IonDual_C9_0167.fastq \
     ~/BaseTIL_code/TCR_bulk_sequencing/UPN003/results/CD3_003

mixcr analyze thermofisher-human-rna-trb-oncomine-sr \
     ~/BaseTIL_code/TCR_bulk_sequencing/UPN003/DIG_003_TCR_TIL_11_IonDual_D9_0168.fastq \
     ~/BaseTIL_code/TCR_bulk_sequencing/UPN003/results/DIG_003
     
   
    
#UPN006
mixcr analyze thermofisher-human-rna-trb-oncomine-sr \
      ~/BaseTIL_code/TCR_bulk_sequencing/UPN006/TIL_006_TCR_TIL_33_IonDual_C12_0191.fastq \
      ~/BaseTIL_code/TCR_bulk_sequencing/UPN006/results/TIL_006
      
mixcr analyze thermofisher-human-rna-trb-oncomine-sr \
      ~/BaseTIL_code/TCR_bulk_sequencing/UPN006/DIG_006_TCR_TIL_34_IonDual_D12_0192.fastq \
      ~/BaseTIL_code/TCR_bulk_sequencing/UPN006/results/DIG_006
 
 
     
#UPN008
mixcr analyze thermofisher-human-rna-trb-oncomine-sr \
     ~/BaseTIL_code/TCR_bulk_sequencing/UPN008/TIL_008_TCR_TIL_27_IonDual_G2_0115.fastq \
     ~/BaseTIL_code/TCR_bulk_sequencing/UPN008/results/TIL_008

mixcr analyze thermofisher-human-rna-trb-oncomine-sr \
     ~/BaseTIL_code/TCR_bulk_sequencing/UPN008/R2CL_008_TCR_TIL_30_IonDual_F5_0138.fastq \
     ~/BaseTIL_code/TCR_bulk_sequencing/UPN008/results/R2CL_008

mixcr analyze thermofisher-human-rna-trb-oncomine-sr \
     ~/BaseTIL_code/TCR_bulk_sequencing/UPN008/R2DIG_008_TCR_TIL_29_IonDual_E5_0137.fastq \
     ~/BaseTIL_code/TCR_bulk_sequencing/UPN008/results/R2DIG_008

mixcr analyze thermofisher-human-rna-trb-oncomine-sr \
     ~/BaseTIL_code/TCR_bulk_sequencing/UPN008/R2FTD_008_TCR_TIL_28_IonDual_H2_0116.fastq \
     ~/BaseTIL_code/TCR_bulk_sequencing/UPN008/results/R2FTD_008

     
#UPN009
mixcr analyze thermofisher-human-rna-trb-oncomine-sr \
     ~/BaseTIL_code/TCR_bulk_sequencing/UPN009/TIL_009_TCR_TIL_18_IonDual_B3_0118.fastq \
     ~/BaseTIL_code/TCR_bulk_sequencing/UPN009/results/TIL_009

mixcr analyze thermofisher-human-rna-trb-oncomine-sr \
     ~/BaseTIL_code/TCR_bulk_sequencing/UPN009/DIG_009_TCR_TIL_21_IonDual_E3_0121.fastq \
     ~/BaseTIL_code/TCR_bulk_sequencing/UPN009/results/DIG_009

mixcr analyze thermofisher-human-rna-trb-oncomine-sr \
     ~/BaseTIL_code/TCR_bulk_sequencing/UPN009/R2L_009_TCR_TIL_20_IonDual_D3_0120.fastq \
     ~/BaseTIL_code/TCR_bulk_sequencing/UPN009/results/R2L_009

mixcr analyze thermofisher-human-rna-trb-oncomine-sr \
     ~/BaseTIL_code/TCR_bulk_sequencing/UPN009/R2K_009_TCR_TIL_19_IonDual_C3_0119.fastq \
     ~/BaseTIL_code/TCR_bulk_sequencing/UPN009/results/R2K_009

     
     
#UPN011_2    
mixcr analyze thermofisher-human-rna-trb-oncomine-sr \
     ~/BaseTIL_code/TCR_bulk_sequencing/UPN011_2/TIL_011_2_TCR_TIL_40_IonDual_B11_0182.fastq \
     ~/BaseTIL_code/TCR_bulk_sequencing/UPN011_2/results/TIL_011_2

mixcr analyze thermofisher-human-rna-trb-oncomine-sr \
     ~/BaseTIL_code/TCR_bulk_sequencing/UPN011_2/SB_011_2TCR_TIL_41_IonDual_C11_0183.fastq \
     ~/BaseTIL_code/TCR_bulk_sequencing/UPN011_2/results/SB_011_2

mixcr analyze thermofisher-human-rna-trb-oncomine-sr \
     ~/BaseTIL_code/TCR_bulk_sequencing/UPN011_2/B_011_2_TCR_TIL_42_IonDual_D11_0184.fastq \
     ~/BaseTIL_code/TCR_bulk_sequencing/UPN011_2/results/B_011_2

     
     