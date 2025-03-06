#! /bin/bash


for i in {1..42};
do 
echo 'S'$i
ln -s ~/BaseTIL_code/Cellranger_multi/run_Cellranger/'S'$i/outs/per_sample_outs/'S'$i/count/sample_alignments.bam ~/BaseTIL_code/Cellranger_multi/run_Cellranger/'S'$i/outs/per_sample_outs/'S'$i/count/possorted_genome_bam.bam
ln -s ~/BaseTIL_code/Cellranger_multi/run_Cellranger/'S'$i/outs/per_sample_outs/'S'$i/count/sample_filtered_feature_bc_matrix ~/BaseTIL_code/Cellranger_multi/run_Cellranger/'S'$i/outs/per_sample_outs/'S'$i/count/filtered_feature_bc_matrix
ln -s ~/BaseTIL_code/Cellranger_multi/run_Cellranger/'S'$i/outs/per_sample_outs/'S'$i/count ~/BaseTIL_code/Cellranger_multi/run_Cellranger/'S'$i/outs/per_sample_outs/'S'$i/outs

done
