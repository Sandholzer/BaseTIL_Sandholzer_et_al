import warnings

warnings.filterwarnings(
    "ignore",
    ".*IProgress not found*",
)
from palmotif import compute_motif, svg_logo
import scanpy as sc
import scirpy as ir
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
import anndata




adatas_tcr = []
sample_chr_vector = ["S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", 
                       "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", "S19", 
                       "S20", "S21", "S22", "S23", "S24", "S25", "S26", "S27", "S28", "S29", 
                       "S30", "S32", "S33", "S34", "S35", "S36", "S37", "S39", 
                       "S40", "S41", "S42"] 

metadata = pd.read_csv('~/BaseTIL_code/data/BaseTIL_metadata_table.csv')
metadata['Batch'] = metadata['Batch'].astype('Int64')

for sample in sample_chr_vector:
    path_tcr_csv = f"~/BaseTIL_code/Cellranger_multi/run_Cellranger/{sample}/outs/per_sample_outs/{sample}/vdj_t/filtered_contig_annotations.csv"
    adata_tcr = ir.io.read_10x_vdj(path_tcr_csv)
    print(f"Amount cells: {len(adata_tcr)}")
    prefix = f"{sample}_"
    adata_tcr.obs_names = [prefix + barcode for barcode in adata_tcr.obs_names]
    adata_tcr.obs['orig.ident'] = f"{sample}"
    adata_tcr.obs['Patient'] = metadata.loc[metadata['Sample'] == sample]['Patient'].item()
    adata_tcr.obs['Type'] = metadata.loc[metadata['Sample'] == sample]['Type'].item()
    adata_tcr.obs['Response'] = metadata.loc[metadata['Sample'] == sample]['Response'].item()
    adata_tcr.obs['Source'] = metadata.loc[metadata['Sample'] == sample]['Source'].item()  
    adata_tcr.obs['batchV'] = metadata.loc[metadata['Sample'] == sample]['Batch'].item()      
    adata_tcr.obs_names_make_unique()
    adatas_tcr.append(adata_tcr)
    
adata_tcr

adata_tcr = anndata.concat(adatas_tcr)

adata_tcr.obs_names_make_unique()


ir.tl.chain_qc(adata_tcr)

_ = ir.pl.group_abundance(adata_tcr, groupby="orig.ident", target_col="chain_pairing")

_ = ir.pl.group_abundance(
   adata_tcr, groupby="orig.ident", target_col="chain_pairing", normalize=True
)

print(adata_tcr.obs["chain_pairing"].value_counts())



observation_ids = ["chain_pairing" ,'is_cell', 'high_confidence', 'multi_chain', 'extra_chains', 'IR_VJ_1_c_call', 'IR_VJ_2_c_call', 'IR_VDJ_1_c_call', 'IR_VDJ_2_c_call', 'IR_VJ_1_consensus_count', 'IR_VJ_2_consensus_count', 'IR_VDJ_1_consensus_count', 'IR_VDJ_2_consensus_count', 'IR_VJ_1_d_call', 'IR_VJ_2_d_call', 'IR_VDJ_1_d_call', 'IR_VDJ_2_d_call', 'IR_VJ_1_duplicate_count', 'IR_VJ_2_duplicate_count', 'IR_VDJ_1_duplicate_count', 'IR_VDJ_2_duplicate_count', 'IR_VJ_1_j_call', 'IR_VJ_2_j_call', 'IR_VDJ_1_j_call', 'IR_VDJ_2_j_call', 'IR_VJ_1_junction', 'IR_VJ_2_junction', 'IR_VDJ_1_junction', 'IR_VDJ_2_junction', 'IR_VJ_1_junction_aa', 'IR_VJ_2_junction_aa', 'IR_VDJ_1_junction_aa', 'IR_VDJ_2_junction_aa', 'IR_VJ_1_locus', 'IR_VJ_2_locus', 'IR_VDJ_1_locus', 'IR_VDJ_2_locus', 'IR_VJ_1_productive', 'IR_VJ_2_productive', 'IR_VDJ_1_productive', 'IR_VDJ_2_productive', 'IR_VJ_1_v_call', 'IR_VJ_2_v_call', 'IR_VDJ_1_v_call', 'IR_VDJ_2_v_call', 'has_ir']
#multiple_observations = adata.obs[observation_ids]
multiple_observations = adata_tcr.obs[observation_ids]

multiple_observations.to_csv("~/BaseTIL_code/saveFiles/TCR_only_annotation.csv", index=True) 

