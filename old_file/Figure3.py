import os, sys
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns
import scanpy as sc
import torch
import scipy
sys.path.append('./')  # uncomment for local import
import tangram as tg
from functions_for_tangram import *
%load_ext autoreload
%autoreload 2
%matplotlib inline

tg.__version__

sc.settings.set_figure_params(dpi=300, facecolor='white')

ad_sp_sham = sc.read("ad_sp_sham.h5ad")
ad_ge_sham= sc.read("ad_ge_sham.h5ad")
ad_map_sham = sc.read("ad_map_sham.h5ad")
ad_ge_sham.obsm = ad_sp_sham.obsm
ad_ge_sham.uns = ad_sp_sham.uns
tg.project_cell_annotations(ad_map_sham, ad_sp_sham, annotation='subclass_label')

ad_sp_4h = sc.read("ad_sp_4h.h5ad")
ad_ge_4h= sc.read("ad_ge_4h.h5ad")
ad_map_4h = sc.read("ad_map_4h.h5ad")
ad_ge_4h.obsm = ad_sp_4h.obsm
ad_ge_4h.uns = ad_sp_4h.uns
tg.project_cell_annotations(ad_map_4h, ad_sp_4h, annotation='subclass_label')

ad_sp_12h = sc.read("ad_sp_12h.h5ad")
ad_ge_12h= sc.read("ad_ge_12h.h5ad")
ad_map_12h = sc.read("ad_map_12h.h5ad")
ad_ge_12h.obsm = ad_sp_12h.obsm
ad_ge_12h.uns = ad_sp_12h.uns
tg.project_cell_annotations(ad_map_12h, ad_sp_12h, annotation='subclass_label')

ad_sp_2d = sc.read("ad_sp_2d.h5ad")
ad_ge_2d= sc.read("ad_ge_2d.h5ad")
ad_map_2d = sc.read("ad_map_2d.h5ad")
ad_ge_2d.obsm = ad_sp_2d.obsm
ad_ge_2d.uns = ad_sp_2d.uns
tg.project_cell_annotations(ad_map_2d, ad_sp_2d, annotation='subclass_label')

ad_sp_6w = sc.read("ad_sp_6w.h5ad")
ad_ge_6w= sc.read("ad_ge_6w.h5ad")
ad_map_6w = sc.read("ad_map_6w.h5ad")
ad_ge_6w.obsm = ad_sp_6w.obsm
ad_ge_6w.uns = ad_sp_6w.uns
tg.project_cell_annotations(ad_map_6w, ad_sp_6w, annotation='subclass_label')

plot_cell_score_timepoints(adata_sp_list=[ad_sp_sham,ad_sp_4h,ad_sp_12h,ad_sp_2d,ad_sp_6w],
                           celltypes=['NewPT1','NewPT2'],
                           cmap="magma_r", 
                           order_cell=True,
                           perc = 0.0001
                          )

plot_cell_annotation_score(ad_sp_sham, ['NewPT1'],x='x', y='y',spot_size= 150, scale_factor=2,perc=0.02, cmap="magma_r")


plot_umap_celltype(ad_sp_2d, color='predicted_celltype', size=4, groups=('NewPT1'), palette="tab20_r")

plot_genes_timepoints_measured(genes=['krt20','spp1','hdc','fbxl13'], 
                      adata_list=[ad_sp_sham, ad_sp_4h, ad_sp_12h, ad_sp_2d, ad_sp_6w], 
                      perc=0.01, 
                      cmap='magma_r')

# genes that are included in the cartana probe list
plot_genes_timepoints(genes=['slc12a1','havcr1'], 
                      adata_list=[ad_ge_sham, ad_ge_4h, ad_ge_12h, ad_ge_2d, ad_ge_6w], 
                      perc=0.001, cmap='magma_r')

plot_genes_timepoints(genes=['foxm1','kcnh8'], 
                      adata_list=[ad_ge_sham, ad_ge_4h, ad_ge_12h, ad_ge_2d, ad_ge_6w], 
                      perc=0.0001, 
                      cmap='magma_r')

plot_genes_timepoints(genes=['cubn','slc34a1'], 
                      adata_list=[ad_ge_sham, ad_ge_4h, ad_ge_12h, ad_ge_2d, ad_ge_6w], 
                      perc=0.01, 
                      cmap='magma_r')

plot_genes_timepoints(genes=['krt20','spp1','hdc','fbxl13'], 
                      adata_list=[ad_ge_sham, ad_ge_4h, ad_ge_12h, ad_ge_2d, ad_ge_6w], 
                      perc=0.01, 
                      cmap='magma_r')

plot_genes_timepoints(genes=['gbe1','rtn4','krt20','prrg4'], 
                      adata_list=[ad_ge_sham, ad_ge_4h, ad_ge_12h, ad_ge_2d, ad_ge_6w], 
                      perc=0.01, 
                      cmap='magma_r')

plot_genes_timepoints(genes=['sult1a1','phyh','spink6','uap1'], 
                      adata_list=[ad_ge_sham, ad_ge_4h, ad_ge_12h, ad_ge_2d, ad_ge_6w], 
                      perc=0.01, 
                      cmap='magma_r')

plot_genes_predicted(["vcam1"], ad_ge_sham, spot_size=150, perc=0.002)
