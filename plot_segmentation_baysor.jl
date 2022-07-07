using Pkg
import Baysor as B
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

df_spatial, gene_names = B.load_df("day2/segmentation.csv")
scale=25;
min_pixels_per_cell = 15;
grid_step = scale / min_pixels_per_cell
bandwidth= scale / 10
df_spatial.ncv_color = convert(Vector{String}, df_spatial.ncv_color)
import Makie as MK
fig1=B.plot_molecules(df_spatial, polygons, color=:ncv_color, 
    markersize=1, size=(5000,6000)
);
MK.save("day2_plot_seg_whole.png", fig1)
fig1
fig2=B.plot_molecules(df_spatial, polygons, color=:ncv_color, 
    markersize=1, size=(1000,1000),
    xlims=(11000,16000), ylims=(8500,13500)
);
MK.save("day2_plot_seg_crop.png", fig2)
fig2
using JLD2
@save "segmentation_day2_ploting.jld"
