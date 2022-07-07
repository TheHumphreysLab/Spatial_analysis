# These are some Julia codes I made for clustering the cartana and merfish data after Baysor cell segmentation

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
using Pkg
B = Baysor;

# Set up python environment to use scanpy in Julia
ENV["PYTHON"]="/home/haojiawu/anaconda3/bin/python"
Pkg.build("PyCall")
using PyCall
using PyPlot
sc = pyimport("scanpy")

# Create a Julia function to take the input files from Baysor and cluster the cell with scanpy
function scanpy_clustering(stat_file::String, count_file::String, n_gene=20)
    cell_stat_df =  DataFrame(CSV.File(stat_file))
    cell_stat_df = @where(cell_stat_df, :n_transcripts .>n_gene)
    cm = DataFrame(CSV.File(count_file, delim = "\t", transpose=true))[cell_stat_df.cell, 2:end]
    adata = sc.AnnData(Matrix(cm))
    adata.var_names = names(cm)
    sc.pp.filter_genes(adata, min_cells=1)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=20, n_pcs=10)
    sc.tl.umap(adata, min_dist=0.2)
    sc.tl.leiden(adata, resolution=0.5)
    return adata
end

# Cell clustering for each IRI samples
sham = scanpy_clustering("sham_stats.csv", "sham.tsv", 20);
hour4 = scanpy_clustering("h4_stats.csv", "h4.tsv", 30);
hour12 = scanpy_clustering("h12_stats.csv", "h12.tsv", 30);
day2 = scanpy_clustering("d2_stats.csv", "d2.tsv", 30);
week6 = scanpy_clustering("wk6_stats.csv", "wk6.tsv", 20);

# Export the umap coordinates and clustering resutls
function write_output(pyObj, file_name::String)
  clusters = convert(Vector{String}, pyObj.obs["leiden"])
  clusters = DataFrame(cluster=clusters)
  CSV.write("cluster_" * file_name * ".csv",clusters)
  umap_coord = DataFrame(get(pyObj.obsm,"X_umap"), :auto)
  CSV.write("umap_" * file_name * ".csv", umap_coord)
end

write_output(sham, "sham")
write_output(hour4, "hour4")
write_output(hour12, "hour12")
write_output(day2, "day2")
write_output(week6, "week6")

# Save the final objects
using JLD2
@save "individual_clustering_obj.jld"

## plot segmentation polygons
