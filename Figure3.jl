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
        y_lims=nothing, marker_size =2)
    if isa(x_lims, Nothing)
        x_lims=(minimum(count_df.x)-0.05*maximum(count_df.x),1.05*maximum(count_df.x))
    end
    if isa(y_lims, Nothing)
        y_lims=(minimum(count_df.y)-0.05*maximum(count_df.y),1.05*maximum(count_df.y))
    end
    if isa(c_map, Nothing)
        c_map= cmap=ColorSchemes.ColorScheme([colorant"gray96",colorant"lemonchiffon1",colorant"red"])
    end
    gene_expr = count_df[!, gene]
    colors = get.(Ref(c_map), (gene_expr .- minimum(gene_expr)) ./ maximum(gene_expr))
    plt_color="#" .* hex.(colors)
    fig = MK.Figure(resolution=canvas_size)
    fig[1, 1] = MK.Axis(fig; xticklabelsize=12, yticklabelsize=12, xticksvisible=false, xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false,
        xgridvisible = false,ygridvisible = false);
    MK.scatter!(count_df.x, count_df.y; color=plt_color,
            strokewidth=0, markersize=marker_size)
    MK.xlims!(MK.current_axis(), x_lims)
    MK.ylims!(MK.current_axis(), y_lims)
    MK.current_figure()
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




