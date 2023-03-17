### Figure 1C-E
using Pkg, StatsBase, CSV, DataFrames, VegaLite, Colors, ColorSchemes, JLD2
import CellScopes as cs
import CairoMakie as MK

sham = load("sham_final.jld2")
cmap = ColorSchemes.ColorScheme([colorant"gray98",colorant"red", colorant"red4"])
p1 = cs.sp_feature_plot(sham, ["Podxl"]; color_keys=cmap, x_lims=(21600,24400), y_lims=(4200,6200), width = 600, height = 500)
p2 = cs.plot_gene_polygons(sham, ["Podxl"]; color_keys= cmap, x_lims=(21600,24400), y_lims=(4200,6200), width = 600, height = 500)
alpha_trans=0.5
anno2 = Dict("Podo" => "green3", "HealthyPT"=>("gray98",alpha_trans),"InjPT"=>("gray98",alpha_trans),"TAL"=>("gray98",alpha_trans),"DCT"=>("gray98",alpha_trans),"CD-PC"=>("gray98",alpha_trans),
            "CD-IC"=>("gray98",alpha_trans),"aEC"=>("gray98",alpha_trans),"gEC"=>("gray98",alpha_trans),"Fib"=>("gray98",alpha_trans),"MC"=>("gray98",alpha_trans),"Immune"=>("gray98",alpha_trans),"Uro"=>("gray98",alpha_trans));
p3 = cs.plot_cell_polygons(sham, "celltype"; cell_colors = anno2, x_lims = (21600,24400), y_lims = (4200,6200), width = 600, height = 500)
