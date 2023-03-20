using Pkg, StatsBase, CSV, DataFrames, VegaLite, Colors, ColorSchemes, JLD2
import CairoMakie as MK
import SpaData as spd
using PyCall
ENV["PYTHON"]="/bin/miniconda3/bin/python"
Pkg.build("PyCall")
sc = pyimport("scanpy")



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
