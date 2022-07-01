def run_spagcn(
    adata,
    x,
    y,
    n_clusters=12,
    p=0.5,
    output_file="SpaGCN_results.h5ad"
):
    x_pixel=adata.obs[x].tolist()
    y_pixel=adata.obs[y].tolist()
    adj=spg.calculate_adj_matrix(x=x_pixel,y=y_pixel, histology=False)
    adata.var_names_make_unique()
    spg.prefilter_genes(adata,min_cells=3)
    spg.prefilter_specialgenes(adata)
    sc.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)
    l=spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)
    r_seed=t_seed=n_seed=100
    res=spg.search_res(adata, adj, l, n_clusters, start=0.7, step=0.1, tol=5e-3, lr=0.05, max_epochs=20, r_seed=r_seed, t_seed=t_seed, n_seed=n_seed)
    clf=spg.SpaGCN()
    clf.set_l(l)
    random.seed(r_seed)
    torch.manual_seed(t_seed)
    np.random.seed(n_seed)
    clf.train(adata,adj,init_spa=True,init="louvain",res=res, tol=5e-3, lr=0.05, max_epochs=200)
    y_pred, prob=clf.predict()
    adata.obs["pred"]= y_pred
    adata.obs["pred"]=adata.obs["pred"].astype('category')
    adj_2d=spg.calculate_adj_matrix(x=x_pixel,y=y_pixel, histology=False)
    refined_pred=spg.refine(sample_id=adata.obs.index.tolist(), pred=adata.obs["pred"].tolist(), dis=adj_2d, shape="hexagon")
    adata.obs["refined_pred"]=refined_pred
    adata.obs["refined_pred"]=adata.obs["refined_pred"].astype('category')
    adata.write_h5ad(output_file)
    
