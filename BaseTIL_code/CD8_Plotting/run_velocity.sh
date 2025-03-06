#!/bin/bash

#SBATCH --job-name=scv                 #This is the name of your job
#SBATCH --cpus-per-task=15                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=32G              #This is the memory reserved per core.
#Total memory reserved: 480GB

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



#export your required environment variables below
#################################################
conda activate velocity2

#add your command lines below
#############################
python .CD8_Plotting/Fig_RNA_velocity.py
