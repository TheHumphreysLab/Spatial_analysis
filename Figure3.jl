using Pkg, StatsBase, CSV, DataFrames, VegaLite, Colors, ColorSchemes, JLD2
import CairoMakie as MK
import SpaData as spd
using PyCall
ENV["PYTHON"]="/bin/miniconda3/bin/python"
Pkg.build("PyCall")
sc = pyimport("scanpy")

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

### data processing
slide_df=prepare_slideseq_data("/home/data/Jia/slideseqv2/czi/Puck_191204_15_healthy.h5ad",
                                "mouse_annot.csv");

visium_df=prepare_visium_data("/home/data/Jia/visium/fsham_137_processed/outs")

cartana_df=prepare_cartana_data("/home/data/Jia/spatial/Sham_Jia/segmentation_cell_stats.csv",
                                "/home/data/Jia/spatial/Sham_Jia/segmentation_counts.tsv")

### plot gene
plot_gene_expression(slide_df, "Umod")
plot_gene_expression(visium_df, "Umod")
plot_gene_expression(cartana_df, "Umod")

### Dot plot to show mixed gene expression in the visium cell type
order_by=["PEC","PT","TAL","MD","DCT","CD-PC","CD-IC","Fib","vSMC","JGA","MC",
            "Macrophage","Leukocyte","vEC"]
p2=dotplot_spatial(slide_df,genes,:celltype2; 
    fig_width=200,cell_order=order_by, fig_height=250, expr_cutoff=1, perc_cutoff=20)

order_by=["Pod","PT","PT-TAL","TAL","DCT","CD-PC","CD-IC","Uro","Fib","Immune"]
p1=dotplot_spatial(visium_df,genes,:celltype2; 
    fig_width=200,cell_order=order_by, fig_height=250,  expr_cutoff=1, perc_cutoff=20)

genes=["Podxl","Slc5a2","Havcr1","Umod","Aqp2","Slc26a4","Krt19","Col1a1","Ren1","Cd74","Flt1","Eln"];
order_by=["Podo","HealthyPT","InjPT","TAL","DCT","CD-PC","CD-IC","Uro","Fib","JGA","Immune", "gEC","vEC"]
p3=dotplot_spatial(cartana_df,genes,:celltype; 
    fig_width=200,cell_order=order_by, fig_height=250, expr_cutoff=1, perc_cutoff=20)

save("slideseq_dot.pdf", p2)
save("visium_dot.pdf",p1)
save("cartana_dot.pdf", p3)

### Plot cell type annotation
alpha_trans=1
anno2 = Dict("Podo" => ("magenta1",alpha_trans), 
                "HealthyPT"=>("green3",alpha_trans),
                "InjPT"=>("#f92874",alpha_trans),
                "TAL"=>("lightslateblue",alpha_trans),
                "DCT"=>("blue",alpha_trans),
                "CD-PC"=>("turquoise1",alpha_trans),
                "CD-IC"=>("#924cfa",alpha_trans),
                "vEC"=>("firebrick",alpha_trans),
                "gEC"=>("dodgerblue",alpha_trans),
                "Fib"=>("#edff4d",alpha_trans),
                "JGA"=>("sienna2",alpha_trans),
                "Immune"=>("darkgreen",alpha_trans),
                "Uro"=>("black",alpha_trans));
cellorder=collect(keys(anno2))
c_map=collect(values(anno2))
plot_annotation(cartana_df, :celltype; marker_size=1.5,cell_order=cellorder, c_map= c_map, y_reverse=false)

cellorder=["PEC","PT","TAL","MD","DCT","CD-PC","CD-IC","Fib","vSMC","JGA","MC",
            "Macrophage","Leukocyte","vEC"]

c_map=["gold","green3","lightslateblue","#339c72","blue","turquoise1","#924cfa","#edff4d",
    "yellow","sienna2","lightgoldenrod1","green","darkgreen","firebrick"];

plot_annotation(slide_df, :celltype2; marker_size=2,cell_order=cellorder, c_map= c_map, y_reverse=false)

cellorder=["Pod","PT","PT-TAL","TAL","DCT","CD-PC","CD-IC","Uro","Fib","Immune"];
c_map=["magenta1","green3","olivedrab1","lightslateblue","blue","turquoise1","#924cfa","black","#edff4d","darkgreen"];
plot_annotation(visium_df, :celltype2; marker_size=9,cell_order=cellorder, c_map= c_map, y_reverse=true )
