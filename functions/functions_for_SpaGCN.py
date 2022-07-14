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
    refined_pred=spg.refine(sample_id=adata.obs.index.tolist(), pred=adata.obs["pred"].tolist(), dis=adj_2d, shape="square")
    adata.obs["refined_pred"]=refined_pred
    adata.obs["refined_pred"]=adata.obs["refined_pred"].astype('category')
    adata.write_h5ad(output_file)
    
def get_spaDEG(
    adata,
    target,
    x,
    y,
    min_in_group_fraction=0.8,
    min_in_out_group_ratio=1,
    min_fold_change=1.5,
    start = None,
    r=None,
    ratio=0.5
):
    if start == None:
        adj_2d=spg.calculate_adj_matrix(x=x_pixel, y=y_pixel, histology=False)
        start, end= np.quantile(adj_2d[adj_2d!=0],q=0.001), np.quantile(adj_2d[adj_2d!=0],q=0.1)
    if r == None:
        r=spg.search_radius(target_cluster=target, 
                    cell_id=adata.obs.index.tolist(), 
                    x=x_pixel, y=y_pixel, 
                    pred=adata.obs["pred"].tolist(), 
                    start=start, end=end, num_min=10, 
                    num_max=14,  max_run=100)
    nbr_domians=spg.find_neighbor_clusters(target_cluster=target,
                                   cell_id=adata.obs.index.tolist(), 
                                   x=adata.obs[x].tolist(), 
                                   y=adata.obs[y].tolist(), 
                                   pred=adata.obs["pred"].tolist(),
                                   radius=r,
                                   ratio=ratio)

    nbr_domians=nbr_domians[0:3]
    de_genes_info=spg.rank_genes_groups(input_adata=adata,
                                target_cluster=target,
                                nbr_list=nbr_domians, 
                                label_col="pred", 
                                adj_nbr=True, 
                                log=True)
    de_genes_info=de_genes_info[(de_genes_info["pvals_adj"]<0.05)]
    filtered_info=de_genes_info
    filtered_info=filtered_info[(filtered_info["pvals_adj"]<0.05) &
                            (filtered_info["in_out_group_ratio"]>min_in_out_group_ratio) &
                            (filtered_info["in_group_fraction"]>min_in_group_fraction) &
                            (filtered_info["fold_change"]>min_fold_change)]
    filtered_info=filtered_info.sort_values(by="in_group_fraction", ascending=False)
    filtered_info["target_dmain"]=target
    filtered_info["neighbors"]=str(nbr_domians)
    print("SVGs for domain ", str(target),":", filtered_info["genes"].tolist())
    return de_genes_info
