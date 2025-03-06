#!/bin/bash

#SBATCH --job-name=velocyto                 #This is the name of your job
#SBATCH --cpus-per-task=10                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=32G              #This is the memory reserved per core.
#Total memory reserved: 320GB

# Are you sure that you need THAT much memory?

#SBATCH --time=168:00:00        #This is the time that your task will run
#SBATCH --qos=1week           #You will run in this queue

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
module load SAMtools/1.15-GCC-10.3.0


#export your required environment variables below
#################################################
conda activate tutorial

#add your command lines below
#############################


for i in {1..42};
do echo 'S'$i
transcriptome="~/BaseTIL_code/Cellranger_multi/reference/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
cellranger_output="~/BaseTIL_code/Cellranger_multi/run_Cellranger/S"$i"/outs/per_sample_outs/S"$i
 
velocyto run10x --samtools-memory 8704 \
                --samtools-threads 30 \
                $cellranger_output \
                $transcriptome

done



