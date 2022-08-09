
function prepare_slideseq_data(h5ad_file, gene_anno_file)
        ctx=sc.read(h5ad_file)
        sc.pp.normalize_total(ctx, target_sum=1e4)
        sc.pp.log1p(ctx)
        coord=get(ctx.obsm,"X_spatial")
        coord=DataFrame(coord, :auto)
        rename!(coord, :x1=> :x, :x2=> :y);
        count_data=ctx.X
        count_data=count_data.toarray();
        count_data=DataFrame(count_data, :auto);
        genenames=collect(ctx.var_names);
        annot=DataFrame(CSV.File(gene_anno_file));
        annot.ensembl_gene_id=string.(annot.ensembl_gene_id);
        annot.mgi_symbol=string.(annot.mgi_symbol);
        annot=filter(:ensembl_gene_id => x -> x in genenames, annot)
        col_remove=[]
        for (i, gene) in enumerate(genenames)
            test=gene in annot.ensembl_gene_id
            if !test
                col_remove=append!(col_remove, i)
            end
        end
        col_n= "x" .* string.(col_remove)
        count_data=select!(count_data, Not(col_n))
        genenames=string.(genenames);
        genenames=genenames[Not(col_remove)]
        annot=spd.reorder(annot, :ensembl_gene_id, genenames)
        rename!(count_data, annot.mgi_symbol, makeunique=true)
        count_data.cells=collect(ctx.obs_names)
        coord.cells=collect(ctx.obs_names)
        coord.celltype=collect(get(ctx.obs,"cell_type"))
        #celltype1=unique(coord.celltype)
        #celltype2=["PCT","vEC","TAL","MC","Fib","gEC","CD-IC","CD-PC","DCT","vSMC","Macrophage","Leukocyte"]
        #coord=spd.mapvalues(coord, :celltype, :celltype2, celltype1, celltype2)
        slide_df=innerjoin(coord, count_data, on = :cells)
        return slide_df
end

function prepare_visium_data(visium_path)
        sham=sc.read_visium(visium_path)        
        sc.pp.normalize_total(sham, target_sum=1e4)
        sc.pp.log1p(sham)
        coord=get(sham.obsm,"spatial")
        coord=DataFrame(coord, :auto)
        rename!(coord, :x1=> :x, :x2=> :y);
        count_data=sham.X
        count_data=count_data.toarray();
        count_data=DataFrame(count_data, :auto);
        genenames=collect(sham.var_names);
        genenames=string.(genenames)
        rename!(count_data, genenames,makeunique=true)
        count_data.cells=collect(sham.obs_names)
        coord.cells=collect(sham.obs_names)
        visium_df=innerjoin(coord, count_data, on = :cells)
        visium_df=innerjoin(coord, count_data, on = :cells)
        return visium_df
end

function prepare_cartana_data(stat_file, count_file)
            coord =  DataFrame(CSV.File(stat_file))
            cm = DataFrame(CSV.File(count_file))
            genes=cm.gene
            cm=cm[!, 2:end]
            sham = sc.AnnData(Matrix(cm))
            sham = sham.transpose()
            sham.var_names = genes
            sham.obs_names = names(cm)
            sc.pp.normalize_total(sham, target_sum=1e4)
            sc.pp.log1p(sham)
            count_data=sham.X
            count_data=DataFrame(count_data, :auto)
            genenames=collect(sham.var_names)
            genenames=string.(genenames)
            rename!(count_data, genenames,makeunique=true)
            count_data.cells=collect(sham.obs_names)
            coord.cell=string.(coord.cell)
            rename!(coord, :cell => :cells)
            coord=coord[!, [:x, :y, :cells]]
            cartana_df=innerjoin(coord, count_data, on=:cells)
            return cartana_df
end

function plot_gene_expression(count_df, gene; 
        c_map=nothing, canvas_size=(600,600),x_lims=nothing, 
        y_lims=nothing, marker_size =2, x_reverse=false, y_reverse=false)
    if isa(x_lims, Nothing)
        x_lims=(minimum(count_df.x)-0.05*maximum(count_df.x),1.05*maximum(count_df.x))
    end
    if isa(y_lims, Nothing)
        y_lims=(minimum(count_df.y)-0.05*maximum(count_df.y),1.05*maximum(count_df.y))
    end
    if isa(c_map, Nothing)
        c_map = ColorSchemes.ColorScheme([colorant"gray94",colorant"pink",colorant"red"])
    end
    if x_reverse
       x_lims = reverse(x_lims)
    end
    if y_reverse
       y_lims = reverse(y_lims)
    end
    gene_expr = count_df[!, gene]
    colors = get.(Ref(c_map), (gene_expr .- minimum(gene_expr)) ./ maximum(gene_expr))
    plt_color="#" .* hex.(colors)
    fig = MK.Figure(resolution=canvas_size)
    fig[1, 1] = MK.Axis(fig; xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
                        xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false,
                        xgridvisible = false, ygridvisible = false)
    MK.scatter!(count_df.x, count_df.y; color=plt_color,
            strokewidth=0, markersize=marker_size)
    MK.xlims!(MK.current_axis(), x_lims)
    MK.ylims!(MK.current_axis(), y_lims)
    MK.current_figure()
end

function plot_annotation(sp_df::DataFrame, col::Union{String, Symbol}; 
        c_map=nothing, canvas_size=(600,600),x_lims=nothing, cell_order=nothing,
        y_lims=nothing, marker_size =2, x_reverse=false, y_reverse=false)
    if isa(cell_order, Nothing)
       celltype=unique(sp_df[!, col])
    end
    celltype=cell_order
    if isa(col, String)
        col=Symbol(col)
    end
    celltype=unique(sp_df[!, col])
    if isa(c_map, Nothing)
        c_map=hex.(Colors.distinguishable_colors(length(celltype), Colors.colorant"#007a10", lchoices=range(20, stop=70, length=15)))
        c_map="#" .* c_map
    end
    sp_df=spd.mapvalues(sp_df, col, :color, celltype, c_map)
    if isa(x_lims, Nothing)
        x_lims=(minimum(sp_df.x)-0.05*maximum(sp_df.x),1.05*maximum(sp_df.x))
    end
    if isa(y_lims, Nothing)
        y_lims=(minimum(sp_df.y)-0.05*maximum(sp_df.y),1.05*maximum(sp_df.y))
    end
    if x_reverse
       x_lims = reverse(x_lims)
    end
    if y_reverse
       y_lims = reverse(y_lims)
    end
    fig = MK.Figure(resolution=canvas_size)
    fig[1, 1] = MK.Axis(fig; xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
                        xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false,
                        xgridvisible = false, ygridvisible = false)
    MK.scatter!(sp_df.x, sp_df.y; color=sp_df.color,
            strokewidth=0, markersize=marker_size)
    MK.xlims!(MK.current_axis(), x_lims)
    MK.ylims!(MK.current_axis(), y_lims)
    MK.current_figure()
end

function dotplot_spatial(sp_df::DataFrame, genes::Union{Vector, String},
        cluster::Union{Symbol, String};expr_cutoff::Union{Float64, Int64}=0,
        x_title="Gene",y_title="Cell type", cell_order::Union{Vector, String, Nothing}=nothing,
        fontsize::Int64=12, color_scheme::String="yelloworangered",reverse_color::Bool=false,
        fig_height::Union{String, Int64}=400, fig_width::Union{String, Int64}=400)
    all_df=DataFrame()
    for (i, gene) in enumerate(genes)
        gene_expr=sp_df[!, gene]
        df = DataFrame()
        df.gene=gene_expr
        df.celltype=string.(sp_df[!, cluster])
        avg_expr=combine(groupby(df, :celltype), :gene => mean => :avg_exp);
        perc_expr=combine(groupby(df, :celltype), :gene => function(x) countmap(x.>expr_cutoff)[:1]*100/length(x) end => :perc_exp)
        df_plt=innerjoin(avg_expr, perc_expr, on = :celltype)
        df_plt.gene.=gene
        all_df=[all_df; df_plt]
    end
    p=all_df |> @vlplot(:circle,
        x={"gene:o", title="Gene", scale={
                domain=genes
            }, axis={labelFontSize=fontsize,titleFontSize=fontsize}},
        y={"celltype:o", title="Cell type",
           scale={
                domain=cell_order
            }, axis={labelFontSize=fontsize,titleFontSize=fontsize}},
        color={"avg_exp:q",
                scale={scheme=color_scheme,reverse=reverse_color}},
        size={"perc_exp:q", legend={symbolFillColor="transparent"}, scale={domain=(0,100)}},
        height= fig_height, width=fig_width
        )
    return p
end
