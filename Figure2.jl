## Here are the codes to reproduce figure2. Function for plotting are in the 'functions' folder.

include("functions/Julia_functions_for_ploting.jl") 
using Pkg
import Baysor
import Colors
import Images
import MAT
import MultivariateStats
import Plots
import CSV
using DataFrames
using DataFramesMeta
using NearestNeighbors
using ProgressMeter
using Statistics
using StatsBase
using JLD2
import Makie as MK
B = Baysor;

pre_data=load("day2_baysor_final.jld2")
polygons= pre_data["polygons"];
df_spatial= pre_data["df_spatial"];

alpha_trans=0.5
anno2 = Dict("Podo" => ("fuchsia",alpha_trans), "HealthyPT"=>("darkolivegreen1",alpha_trans),"InjPT"=>("darkolivegreen1",alpha_trans),"TAL"=>("cyan",alpha_trans),"DCT"=>("coral",alpha_trans),"CD-PC"=>("gray95",alpha_trans),
            "CD-IC"=>("gray95",alpha_trans),"aEC"=>("gray95",alpha_trans),"gEC"=>("blue",alpha_trans),"Fib"=>("gray95",alpha_trans),"MC"=>("darkgreen",alpha_trans),"Immune"=>("gray95",alpha_trans),"Uro"=>("gray95",alpha_trans))
@time plot_genes_overlaid_cells(df_spatial, polygons; 
    genes=["Podxl","Ehd3","Ren1"], colors=["fuchsia","blue","darkgreen"],
    markersize=2,annotation=:celltype, ann_colors=anno2, is_noise=:is_noise,noise_kwargs=(markersize=0, color="transparent"),
    show_legend=true
)
