# function to run tangram on each sample
def run_tangram(sc_mat, sc_genes, sc_barcodes, sc_meta, sp_mat, sp_genes, sp_barcodes, sp_meta, markers):
    adata=sc.read_mtx(sc_mat)
    adata=adata.transpose()
    gene1 = pd.read_table(sc_genes, header=None)
    cell1 = pd.read_table(sc_barcodes, header=None)
    adata.var_names=gene1.iloc[:,0]
    adata.obs_names=cell1.iloc[:,0]
    metadata = pd.read_csv(sc_meta, index_col=0)
    metadata.dropna(inplace = True)
    adata.obs['subclass_label'] = metadata["name"]
    ad_sc = adata
    sc.pp.normalize_total(ad_sc)
    adata=sc.read_mtx(sp_mat)
    adata=adata.transpose()
    gene1 = pd.read_table(sp_genes, header=None)
    cell1 = pd.read_table(sp_barcodes, header=None)
    adata.var_names=gene1.iloc[:,0]
    adata.obs_names=cell1.iloc[:,0]
    metadata = pd.read_csv(sp_meta, index_col=0)
    metadata.dropna(inplace = True)
    adata.obs['x'] = metadata["x"]
    adata.obs['y'] = metadata["y"]
    ad_sp = adata
    df_genes = pd.read_csv(markers, index_col=0)
    markers = np.reshape(df_genes.values, (-1, ))
    markers = list(markers)
    tg.pp_adatas(ad_sc, ad_sp, genes=markers)
    assert ad_sc.uns['training_genes'] == ad_sp.uns['training_genes']
    ad_map = tg.map_cells_to_space(
    adata_sc=ad_sc,
    adata_sp=ad_sp,
    device='cpu',
    #device='cuda:0',
    )
    ad_map.write("ad_map.h5ad")
    ad_sc.write("ad_sc.h5ad")
    ad_sp.write("ad_sp.h5ad")

# function to scale the data locally or globally (adapted from tangram)
def construct_obs_plot(df_plot, adata, perc=0.001, suffix=None):
    df_plot = df_plot.clip(df_plot.quantile(perc), df_plot.quantile(1 - perc), axis=1)
    df_plot = (df_plot - df_plot.min()) / (df_plot.max() - df_plot.min())

    if suffix:
        df_plot = df_plot.add_suffix(" ({})".format(suffix))
    adata.obs = pd.concat([adata.obs, df_plot], axis=1)


# function for ploting multiple genes from the predicted gene expression on a single timepoint
def plot_genes_predicted(
    genes,
    adata_predicted,
    x = "x",
    y = "y",
    spot_size = 200,
    scale_factor=0.1, 
    cmap="inferno_r", 
    perc=0.001
):
    adata_predicted.obs.drop(
        ["{} (Inferred)".format(gene) for gene in genes],
        inplace=True,
        errors="ignore",
        axis=1,
    )
    adata_predicted.var.index = [g.lower() for g in adata_predicted.var.index]
    coords = [[x,y] for x,y in zip(adata_predicted.obs[x].values,adata_predicted.obs[y].values)]
    adata_predicted.obsm['spatial'] = np.array(coords)
    df = pd.DataFrame(
        data=np.array(adata_predicted[:, genes].X),
        columns=genes,
        index=adata_predicted.obs.index,
    )
    construct_obs_plot(df, adata_predicted, perc=perc, suffix="Inferred")
    genes = [s + " (Inferred)" for s in genes]
    sc.pl.spatial(
    adata_predicted,
    spot_size=spot_size,
    scale_factor=scale_factor,
    color=genes,
    frameon=False,
    show=False,
    cmap=cmap
        )
    

# function for ploting the predicted genes expression across IRI timepoints
def plot_genes_timepoints(
    genes,
    adata_list,
    order_cell=False,
    x = "x",
    y = "y",
    spot_size = 200,
    scale_factor=0.1, 
    cmap="inferno_r", 
    perc=0.001
):
    adata_predicted= adata_list[0].concatenate(adata_list[1:len(adata_list)])
    adata_predicted.var.index = [g.lower() for g in adata_predicted.var.index]
    coords = [[x,y] for x,y in zip(adata_predicted.obs[x].values,adata_predicted.obs[y].values)]
    adata_predicted.obsm['spatial'] = np.array(coords)
    df = pd.DataFrame(
        data=np.array(adata_predicted[:, genes].X),
        columns=genes,
        index=adata_predicted.obs.index,
    )
    construct_obs_plot(df, adata_predicted, perc=perc, suffix="Inferred")
    adata2=[]
    for hx, adata1 in enumerate(adata_list):
        new_adata=adata_predicted[adata_predicted.obs.index.isin([s + "-" + str(hx) for s in adata1.obs_names])]
        adata2.append(new_adata)
    fig = plt.figure(figsize=(len(genes) * 7, len(adata2) * 7))
    gs = GridSpec(len(adata2), len(genes), figure=fig)
    for hx, adata3 in enumerate(adata2):
        for vx, gene in enumerate(genes):
            ax = fig.add_subplot(gs[hx, vx])
            title_text = "" if hx!=0 else None
            clbr = "right" if vx==len(genes)-1 and hx==len(adata2)-1 else None
            sc.pl.spatial(
                adata3,
                spot_size=spot_size,
                scale_factor=scale_factor,
                color=["{} (Inferred)".format(gene)],
                frameon=False,
                show=False,
                ax=ax,
                cmap=cmap,
                title= title_text,
                sort_order=order_cell,
                vmin = 0,
                vmax = 1,
                colorbar_loc=clbr)

# function for ploting the original genes expression across IRI timepoints
def plot_genes_timepoints_measured(
    genes,
    adata_list,
    order_cell=False,
    x = "x",
    y = "y",
    spot_size = 200,
    scale_factor=0.1, 
    cmap="inferno_r", 
    perc=0.001
):
    adata_measured= adata_list[0].concatenate(adata_list[1:len(adata_list)])
    adata_measured.X = adata_measured.X.toarray()
    adata_measured.var.index = [g.lower() for g in adata_measured.var.index]
    coords = [[x,y] for x,y in zip(adata_measured.obs[x].values,adata_measured.obs[y].values)]
    adata_measured.obsm['spatial'] = np.array(coords)
    data = []
    for ix, gene in enumerate(genes):
        if gene not in adata_measured.var.index:
            data.append(np.zeros_like(np.array(adata_measured[:, 0].X).flatten()))
        else:
            data.append(np.array(adata_measured[:, gene].X).flatten())

    df = pd.DataFrame(
        data=np.array(data).T, columns=genes, index=adata_measured.obs.index,
    )
    construct_obs_plot(df, adata_measured, perc=perc, suffix="Measured")
    cmap1 = plt.cm.get_cmap(cmap)
    rgba = cmap1(0)
    na_color=matplotlib.colors.to_hex(rgba)
    adata2=[]
    for hx, adata1 in enumerate(adata_list):
        new_adata=adata_measured[adata_measured.obs.index.isin([s + "-" + str(hx) for s in adata1.obs_names])]
        adata2.append(new_adata)
    fig = plt.figure(figsize=(len(genes) * 7, len(adata2) * 7))
    gs = GridSpec(len(adata2), len(genes), figure=fig)
    for hx, adata3 in enumerate(adata2):
        for vx, gene in enumerate(genes):
            ax = fig.add_subplot(gs[hx, vx])
            title_text = "" if hx!=0 else None
            clbr = "right" if vx==len(genes)-1 and hx==len(adata2)-1 else None
            sc.pl.spatial(
                adata3,
                spot_size=spot_size,
                scale_factor=scale_factor,
                color=["{} (Measured)".format(gene)],
                frameon=False,
                show=False,
                ax=ax,
                cmap=cmap,
                title= title_text,
                sort_order=order_cell,
                vmin = 0,
                vmax = 1,
                colorbar_loc=clbr,
                na_color=na_color)

