import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import spateo as st
import glob
import scanpy as sc
import anndata
%matplotlib inline

## Read imputed spatial data
chunks = pd.read_csv(data_path,index_col=0,chunksize=1000000)
sp_count = pd.concat(chunks)
## Read metadata
sp_meta = pd.read_csv(data_path_meta)
cci = ["C3.cell", "FR.PT"]
sp_meta_cci = sp_meta[sp_meta["celltype3"].isin(cci)]

color_dict = {
 'FR.PT': "#01b0f0",
 'C3.cell': "#ca5d41"}
obs = pd.DataFrame()
obs["Annotation"] = sp_meta_cci.celltype3
sp_adata = anndata.AnnData(sp_count.values,obs = obs)
sp_adata.obs_names =  sp_count.index.astype(str)
sp_adata.var_names =  sp_count.columns 

df_loc = pd.concat([sp_meta_cci.x, sp_meta_cci.y], axis = 1)
sp_adata.obsm["spatial"] = np.array(df_loc)
sp_adata.uns["__type"] = 'UMI' 
sp_adata.uns["color_key"] = color_dict
## Find spatially adjacent celltypes, update sp_adata
weights_graph, distance_graph, sp_adata = st.tl.weighted_spatial_graph(
	sp_adata,
	n_neighbors=10,
	fixed='n_neighbors',
)
sender_ct = 'C3.cell'
receptor_ct = 'FR.PT'

result = st.tl.find_cci_two_group(sp_adata,
                               path='path/to/ligand-receptor-database',
                               species='mouse',
                               group='Annotation',
                               sender_group=sender_ct,
                               receiver_group=receptor_ct,
                               filter_lr='outer',
                               min_pairs=0,
                               min_pairs_ratio=0,
                               top=20,)
							
df = df.loc[df['lr_co_exp_num'] > 0].sort_values('lr_co_exp_ratio', ascending=False)[0:10]
data = df.iloc[:, [0, 1, 7]]
data = data.drop_duplicates(subset = ["from","to"]).reset_index(drop=True)
data = data[data["from"]!=data["to"]]
plot_data = data.pivot(index="from", columns="to",values="lr_co_exp_ratio").fillna(0)
## Plot ligand-receptor pairs
fig = plt.figure()
fig.set_size_inches(5.5, 3)
x_label = list(plot_data.columns.tolist())
y_label = list(plot_data.index)
ax = sns.heatmap(plot_data,
                 cmap="magma_r",
                 yticklabels=y_label,
                 linecolor='grey',
                 linewidths=0.3,
                 annot_kws={'size': 12, 'weight': 'bold', },
                 xticklabels=x_label,
                 mask=(plot_data < 0.01))
plt.gcf().subplots_adjust(bottom=0.3)
plt.xlabel("Receptor in FR-PT")
plt.ylabel("Ligand in C3+ cell")
ax.set_xticklabels(x_label, rotation=45, ha="right")
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('none')
plt.tight_layout()

## Plot adjacent cell types
chunks = pd.read_csv(data_path,index_col=0,chunksize=1000000)
sp_count_immune = pd.concat(chunks)
cells = ["C3.cell", "FR.PT","Macrophage","Ly6e.lymphocyte","Tcell","Cxcl10.cell","Plasma cell"]
sp_meta_immune = sp_meta[sp_meta["celltype3"].isin(cells)]
color_dict = {
 'Macrophage': "#70a845",
 'Ly6e.lymphocyte': "#b69141",
 'FR.PT': "#01b0f0",
 'Tcell': "#c8577d",
 'C3.cell': "#ca5d41",
 'Cxcl10.cell': '#994986',
 'Plasma cell': '#4bc3b7'}
obs = pd.DataFrame()
obs["Annotation"] = sp_meta_6w_immune.celltype3
sp_adata = anndata.AnnData(sp_count_immune.values,obs = obs)
sp_adata.obs_names =  sp_count_immune.index.astype(str)
sp_adata.var_names =  sp_count_immune.columns 
df_loc = pd.concat([sp_meta_6w_immune.x, sp_meta_6w_immune.y], axis = 1)
sp_adata.obsm["spatial"] = np.array(df_loc)
sp_adata.uns["__type"] = 'UMI'
sp_adata.uns["color_key"] = color_dict

weights_graph, distance_graph, sp_adata = st.tl.weighted_spatial_graph(
	sp_adata,
	n_neighbors=10,
	fixed='n_neighbors',
)
st.pl.plot_connections(
	sp_adata,
	cat_key='Annotation',
	save_show_or_return='all',
	title_str=" ",
	title_fontsize=10,
	label_fontsize=10,
	colormap=sp_adata.uns['color_key'],
	figsize=(5,5),
	save_kwargs = {"path": "ajd_ctp_immune", "dpi": 500, "ext":'jpg', "transparent": True, "close": True,
			"verbose": True},
)
