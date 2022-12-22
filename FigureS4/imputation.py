####  imputation use gimVI, spaGE, and Tangram method

import sys
import scanpy
import scanpy as sc
import anndata
from anndata import AnnData
import pathlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import tangram as tg
import scipy.stats as st
from scipy.stats import spearmanr
from scvi.external import GIMVI
import SpaGE
from SpaGE.main import SpaGE
import loompy
import os
import argparse

## functions
def Norm_rna(x):
	return np.log(((x/np.sum(x))*1000000)+1)

def Norm_spatial(x):
	return np.log(((x/np.sum(x))*np.median(cell_count)) + 1)

def get_gene_id(data,gene):
	for idx,i in enumerate(data):
		if i==gene:
			return idx
		
if __name__=="__main__":
	
	parse = argparse.ArgumentParser(description=__doc__)
	
	parse.add_argument("-seq", "--sc_rna", type=str, action="store", help="input file for scRNA-seq data")
	parse.add_argument("-sm", "--sc_rna_meta", type=str, action="store", help="input file for scRNA-seq meta data")
	parse.add_argument("-scs", "--segmentation_cell_stats", type=str, action="store", help="segmentation file of spatial data")
	parse.add_argument("-sc", "--segmentation_counts", type=str, action="store", help="segmentation count file of spatial data")
	parse.add_argument("-o", "--output_dir", type=str, help='results directory',default='./results')
	parse.add_argument("-m", "--method", type=str, help='imputation method',default='spaGE')
	args = parse.parse_args()
	
## read rna data
chunks = pd.read_table(args.sc_rna,index_col=0,chunksize=1000000)
rna_data = pd.concat(chunks)
gene_count = np.sum(rna_data>0, axis=1)
rna_data = rna_data.loc[gene_count >=10, :]
rna_data = rna_data.apply(Norm_rna,axis=0)

## read spatial data
spatial_meta = pd.read_csv(args.segmentation_cell_stats)
spatial_data = pd.read_csv(args.segmentation_counts,sep='\t')
spatial_data = spatial_data.set_index('gene') # a total of 193 genes, with 169 overlapped genes
cell_count = np.sum(spatial_data,axis=0)
spatial_data = spatial_data.apply(Norm_spatial,axis=0)

## adata format is only needed for tangram and gimVI method
obs = pd.DataFrame()
obs['x_coord'] = spatial_meta.x
obs['y_coord'] = spatial_meta.y
obs['labels'] = spatial_meta.cluster
spatial_adata =  anndata.AnnData(spatial_data.T.values, obs = obs)

spatial_adata.var_names = spatial_data.index.tolist()
spatial_adata.obs_names = spatial_data.columns.tolist()

seq_meta = pd.read_csv(args.sc_rna_meta,compression='gzip',sep='\t')
seq_meta = seq_meta[seq_meta['Group']=='Control']

rna_adata = anndata.AnnData(rna_data.T)
rna_adata.obs['cell_subclass'] = seq_meta.celltype
rna_adata.var_names = rna_data.index.tolist()
rna_adata.obs_names = rna_data.columns.tolist()

overlap_genes = list(set(rna_adata.var_names) & set(spatial_adata.var_names))
rna_adata = rna_adata[:,overlap_genes].copy()#only use genes in both datasets

## define the genes you want to impute 

gene_set = ['Slc12a1','Ncam1','Aqp4']
plt.style.use('dark_background')

if args.method == "spaGE":
	
	for i in gene_set:
		Imp_Genes = SpaGE(spatial_data.T.drop(i,axis=1),rna_data.T,n_pv=30,
								genes_to_predict = [i])
		
		fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,5))
		ax1.axis('off')
		cmap = spatial_data.T[i]
		cmap[cmap > np.percentile(cmap,99)] = np.percentile(cmap,99)
		ax1.scatter(spatial_meta['x'],spatial_meta['y'],s=0.1,c=cmap)
		ax1.set_title('Measured ' + i, fontsize = 12)
		ax1.set_ylabel(i)
		
		ax2.axis('off')
		cmap = Imp_Genes[i]
		cmap[cmap > np.percentile(cmap,99)] = np.percentile(cmap,99)
		ax2.scatter(spatial_meta['x'],spatial_meta['y'],s=0.1,c=cmap)
		ax2.set_title('spaGE ' + i, fontsize = 12)
		save_file = i + '_spaGE_train.jpg'
		plt.savefig(os.path.join(args.output_dir,save_file),dpi=500)
		
elif args.method == "gimVI":
	
	seq_gene_names = rna_adata.var_names
	n_genes = rna_adata.n_vars
	n_train_genes = int(n_genes*train_size)
	
	#remove cells with no counts
	scanpy.pp.filter_cells(spatial_adata, min_counts= 1)
	scanpy.pp.filter_cells(rna_adata, min_counts = 1)
	
	#Sets up AnnData object for scvi models
	GIMVI.setup_anndata(spatial_adata, labels_key='labels')
	GIMVI.setup_anndata(rna_adata)
	
	#create model for training
	model = GIMVI(rna_adata, spatial_adata)
	model.train(200)

	#plot imputation for genes
	for i in gene_set:
		
		fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(10,5))
		
		# Plot groundtruth
		ax1.axis('off')
		cmap = spatial_data.T[i]
		cmap[cmap > np.percentile(cmap,99)] = np.percentile(cmap,99)
		ax1.scatter(spatial_meta['x'],spatial_meta['y'],s=0.1,c=cmap)
		ax1.set_title('Measured ' + i, fontsize = 12)
		ax1.set_ylabel(i)
		
		ax2.axis('off')
		_, imputed = model.get_imputed_values(normalized=True)
		gene_id = get_gene_id(spatial_adata.var_names,i)
		cmap = imputed[:, gene_id]
		cmap[cmap > np.percentile(cmap,99)] = np.percentile(cmap,99)
		ax2.scatter(spatial_meta['x'],spatial_meta['y'],s=0.1,c=cmap)
		ax2.set_title('gimVI ' + i, fontsize = 12)
		save_file = i + '_gimVI.jpg'
		plt.savefig(os.path.join(args.output_dir,save_file),dpi=500)
		
elif args.method == "tangram":
	
	spatial_adata = spatial_adata[:,overlap_genes].copy()
	scanpy.pp.filter_cells(spatial_adata, min_counts= 1)
	scanpy.pp.filter_cells(rna_adata, min_counts = 1)
	
	tg.pp_adatas(rna_adata, spatial_adata, genes=overlap_genes)
	ad_map = tg.map_cells_to_space(rna_adata, spatial_adata, density_prior='uniform',num_epochs=500,device="cpu")
	ad_ge = tg.project_genes(adata_map=ad_map, adata_sc=rna_adata)
	
	adata_predicted = ad_ge

	for i in gene_set:
		fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(10,5))
		
		ax1.axis('off')
		cmap = spatial_data.T[i]
		cmap[cmap > np.percentile(cmap,99)] = np.percentile(cmap,99)
		ax1.scatter(spatial_meta['x'],spatial_meta['y'],s=0.1,c=cmap)
		ax1.set_title('Measured ' + i, fontsize = 12)
		ax1.set_ylabel(i)
		
		ax2.axis('off')
		cmap = np.array(adata_predicted[:, i].X).flatten()
		cmap[cmap > np.percentile(cmap,99)] = np.percentile(cmap,99)
		ax2.scatter(spatial_meta['x'],spatial_meta['y'],s=0.1,c=cmap)
		ax2.set_title('Tamgram ' + i, fontsize = 12)
		save_file = i + '_Tangram.jpg'
		plt.savefig(os.path.join("result",save_file),dpi=500)
	
	
