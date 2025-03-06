#!/bin/bash

#SBATCH --job-name=pyScenic                #This is the name of your job
#SBATCH --cpus-per-task=20                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=16G              #This is the memory reserved per core.
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
#module load Boost/1.69.0-foss-2018b-Python-3.6.6
#module load GDAL/2.4.1-foss-2018b-Python-3.6.6
#module load protobuf/3.13.0-GCCcore-7.3.0
#module load R/4.1.2-foss-2018b-bare
#module load R/4.2.2-foss-2021a-bare
#module load R/4.3.0-foss-2021a-bare

source ~/miniconda3/etc/profile.d/conda.sh
conda activate ~/mambaforge/envs/pyscenic

#export your required environment variables below
#################################################


#add your command lines below
#############################

#To correct loom files for row and colnames
python ~/BaseTIL_code/Scenic/pyScenic_loomcorrection.py



samp=$1
echo 'pyscenic run for' $samp

# setup directories
f_loom_path_scen='~/BaseTIL_code/Scenic/pySenic_data/expression_matrix/CD8_'$samp'_corrected.loom'
adj_out='~/BaseTIL_code/Scenic/pySenic_data/output/adj_'$samp'.tsv'
reg_out='~/BaseTIL_code/Scenic/pySenic_data/output/reg_'$samp'.csv'
f_db_names='~/BaseTIL_code/Scenic/pySenic_data/databases/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather ~/BaseTIL_code/Scenic/pySenic_data/databases/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather ~/BaseTIL_code/Scenic/pySenic_data/databases/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather ~/BaseTIL_code/Scenic/pySenic_data/databases/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'
f_pyscenic_output='~/BaseTIL_code/Scenic/pySenic_data/output/pyscenic_'$samp'_output.loom'


# run pysenic steps

pyscenic grn \
          $f_loom_path_scen \
          ~/BaseTIL_code/Scenic/pySenic_data/databases/allTFs_hg38.txt \
          --output $adj_out\
          --num_workers 20


pyscenic ctx \
    $adj_out \
    $f_db_names \
    --annotations_fname ~/BaseTIL_code/Scenic/pySenic_data/databases/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
    --expression_mtx_fname $f_loom_path_scen \
    --output $reg_out \
    --mask_dropouts \
    --num_workers 20


pyscenic aucell \
    $f_loom_path_scen \
    $reg_out \
    --output $f_pyscenic_output \
    --num_workers 20


