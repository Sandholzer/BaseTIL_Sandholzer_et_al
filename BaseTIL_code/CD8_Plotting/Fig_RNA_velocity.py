## -----------------------------------------------------------------------------------------------------------
#load dependencies

import re
import scanpy as sc
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd
import scvelo as scv
import matplotlib.pyplot as pl
import anndata as ad
import cellrank as cr

scv.set_figure_params(style="scvelo")
pl.rcParams["figure.figsize"] = (5,4)
pl.rcParams['axes.grid'] = False

np.random.seed(19680801)
sc.settings.figdir = "~/BaseTIL_code/CD8_Plotting/Graphs"
scv.settings.figdir = "~/BaseTIL_code//CD8_Plotting/Graphs"
scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2
pl.style.use('ggplot')

adata = sc.read_h5ad(filename = "~/BaseTIL_code/saveFiles/CD8_tumor_reactive.h5ad")
adata.obs['TR_anno'] = adata.obs['TR_anno'].astype('category',copy=False)

## -----------------------------------------------------------------------------------------------------------
#plot umap 
sc.pl.umap(adata, color=['TR_anno'], frameon=False, save='cluster.pdf') #, palette=Colorss


#load pre-ACT and post-ACT samples
sample_name = ["S1", "S22", "S26", "S27", "S32", "S34", "S37", "S39", "S40", "S41", "S42",  "S5",  "S6",  "S7"]   

ldata_list = []

for sample in sample_name:
  print(f"Start {sample}")
  ldata= scv.read(f'~/BaseTIL/Cellranger_multi/run_Cellranger/{sample}/outs/per_sample_outs/{sample}/velocyto/{sample}.loom', cache=True)
  barcodes = [bc.split(':')[1] for bc in ldata.obs.index.tolist()]
  barcodes = [f'{sample}_' + bc[0:len(bc)-1] + '-1' for bc in barcodes]
  ldata.obs.index = barcodes
  ldata.var_names_make_unique()
  ldata_list.append(ldata)

# concatenate the looms
ldata = ldata_list[0].concatenate(ldata_list[1:len(ldata_list)], index_unique=None)
del ldata_list

#merge spliced/unspliced matrix with existing object
adata = scv.utils.merge(adata, ldata)

# check splicing proportions
scv.pl.proportions(adata, groupby='TR_anno', save=f"glo_Spliced.pdf")


## clean some genes
ignoreGenes = pd.read_csv("~/BaseTIL/data/all.gene.ignore.df.txt.gz", compression="gzip", sep="\t", header=0, index_col=0)
flag = [not i in list(ignoreGenes['geneSymbol']) for i in list(adata.var.index)]
adata = adata[:,flag]
flag = [not bool(re.match('^RP[LS]', i)) for i in adata.var_names]
adata = adata[:,flag]
adata = adata[:,adata.var_names != "MALAT1"]

# normalize and filter
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000) 
sc.pp.neighbors(adata)
scv.pp.moments(adata, n_pcs=None, n_neighbors=None)


#calculate splicing dynamics
scv.tl.recover_dynamics(adata,var_names='all', n_jobs=12)
scv.tl.velocity(adata, mode="dynamical", n_jobs=12)
scv.tl.velocity_graph(adata, n_jobs=12)


#plot velocity embeddings on UMAP
Colorss=["#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C"]
scv.set_figure_params(style="scvelo")
pl.rcParams["figure.figsize"] = (5,4)
pl.rcParams['axes.grid'] = False


scv.pl.velocity_embedding_stream(adata, basis="umap",title="", color="TR_anno", save="TR_velo_stream.pdf", size=40, palette=Colorss, frameon=False, legend_loc='none')
scv.pl.velocity_embedding_grid(adata, basis='umap', color='TR_anno', save="glo_embedding_grid.pdf", title='', scale=0.5, density=0.5, arrow_size=2, linewidth=1, arrow_length=2)
