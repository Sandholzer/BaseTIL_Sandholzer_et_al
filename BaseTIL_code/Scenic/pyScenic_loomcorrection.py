import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp

types_group = ['full', 'Tumor', 'PreREP', 'Expanded_TILs', 'PBMC_7dpt', 'Rebiopsy1', 'Rebiopsy2']

for type_value in types_group:
    input=f'~/BaseTIL_code/Scenic/pySenic_data/expression_matrix/CD8_{type_value}.loom'
    output=f'~/BaseTIL_code/Scenic/pySenic_data/expression_matrix/CD8_{type_value}_corrected.loom'
    
    adata = sc.read_loom( input )
    row_attrs = {
    "Gene": np.array(adata.var.index) ,
    }
    col_attrs = {
    "CellID": np.array(adata.obs.index) ,
    "nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
    }
    
    for key,values in adata.obs.items():
        col_attrs[key]=np.array(values)
    
    lp.create( output, adata.X.transpose(), row_attrs, col_attrs )
