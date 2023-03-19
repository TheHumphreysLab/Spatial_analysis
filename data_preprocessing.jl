import CellScopes as cs
using CSV, DataFrames
using SparseArrays

### 1. sham ###
count_molecules =  DataFrame(CSV.File("/mnt/sdc/new_analysis_cellscopes/seg_files/Sham_Jia/segmentation.csv"))
count_df =  DataFrame(CSV.File("/mnt/sdc/new_analysis_cellscopes/seg_files/Sham_Jia/segmentation_counts.tsv"))
count_cells =  DataFrame(CSV.File("/mnt/sdc/new_analysis_cellscopes/seg_files/Sham_Jia/segmentation_cell_stats.csv"))
import Baysor as B
scale=30
min_pixels_per_cell = 15
grid_step = scale / min_pixels_per_cell
bandwidth= scale / 10
polygons = B.boundary_polygons(count_molecules, count_molecules.cell, grid_step=grid_step, bandwidth=bandwidth)
genes = count_df.gene
cells = string.(names(count_df)[2:end])
count_df = count_df[!, 2:end]
count_df = convert(SparseMatrixCSC{Int64, Int64},Matrix(count_df))
raw_count = cs.RawCountObject(count_df, cells, genes)
count_molecules.cell = string.(count_molecules.cell)
count_cells.cell = string.(count_cells.cell)
sham = cs.CartanaObject(count_molecules, count_cells, raw_count;
    prefix = "sham", min_gene = 0, min_cell = 3)
sham = cs.normalize_object(sham)
sham = cs.scale_object(sham)
sham = cs.find_variable_genes(sham; nFeatures=100)
sham = cs.run_pca(sham;  method=:svd, pratio = 1, maxoutdim = 20)
sham.metaData.cluster = string.(sham.spmetaData.cell.cluster)
sham.polygonData = polygons
sham = cs.polygons_cell_mapping(sham);
sham = cs.generate_polygon_counts(sham);
cs.save(sham; filename="sham.jld2")
spaGE_path = "/mnt/sdc/new_analysis_cellscopes/SpaGE"
data_path = "/mnt/sdc/new_analysis_cellscopes/for_imputate/IRI_sham/"
sham = cs.load(filename="sham.jld2")
sham = cs.run_spaGE(sham, data_path, spaGE_path)
center=[13806.73, 23484.30]
sham = cs.compute_kidney_coordinates(sham, center)
cs.save(sham; filename="sham_imputed.jld2")

### 2. IRI4h ###
count_molecules =  DataFrame(CSV.File("/mnt/sdc/new_analysis_cellscopes/seg_files/Hour4_Jia/segmentation.csv"))
count_df =  DataFrame(CSV.File("/mnt/sdc/new_analysis_cellscopes/seg_files/Hour4_Jia/segmentation_counts.tsv"))
count_cells =  DataFrame(CSV.File("/mnt/sdc/new_analysis_cellscopes/seg_files/Hour4_Jia/segmentation_cell_stats.csv"))
import Baysor as B
scale=30
min_pixels_per_cell = 15
grid_step = scale / min_pixels_per_cell
bandwidth= scale / 10
polygons = B.boundary_polygons(count_molecules, count_molecules.cell, grid_step=grid_step, bandwidth=bandwidth)
genes = count_df.gene
cells = string.(names(count_df)[2:end])
count_df = count_df[!, 2:end]
count_df = convert(SparseMatrixCSC{Int64, Int64},Matrix(count_df))
raw_count = cs.RawCountObject(count_df, cells, genes)
count_molecules.cell = string.(count_molecules.cell)
count_cells.cell = string.(count_cells.cell)
hour4 = cs.CartanaObject(count_molecules, count_cells, raw_count;
    prefix = "hour4", min_gene = 0, min_cell = 3)
hour4 = cs.normalize_object(hour4)
hour4 = cs.scale_object(hour4)
hour4 = cs.find_variable_genes(hour4; nFeatures=100)
hour4 = cs.run_pca(hour4;  method=:svd, pratio = 1, maxoutdim = 20)
hour4.metaData.cluster = string.(hour4.spmetaData.cell.cluster)
hour4.polygonData = polygons
hour4 = cs.polygons_cell_mapping(hour4)
hour4 = cs.generate_polygon_counts(hour4)
cs.save(hour4; filename="hour4.jld2")
spaGE_path = "/mnt/sdc/new_analysis_cellscopes/SpaGE"
data_path = "/mnt/sdc/new_analysis_cellscopes/for_imputate/IRI_4h/"
hour4 = cs.load(filename="hour4.jld2")
hour4 = cs.run_spaGE(hour4, data_path, spaGE_path)
center=[14708.91, 22491.68]
hour4 = cs.compute_kidney_coordinates(hour4, center)
cs.save(hour4; filename="hour4_imputed.jld2")

### 3. IRI12h ###
count_molecules =  DataFrame(CSV.File("/mnt/sdc/new_analysis_cellscopes/seg_files/Hour12_Jia/segmentation.csv"))
count_df =  DataFrame(CSV.File("/mnt/sdc/new_analysis_cellscopes/seg_files/Hour12_Jia/segmentation_counts.tsv"))
count_cells =  DataFrame(CSV.File("/mnt/sdc/new_analysis_cellscopes/seg_files/Hour12_Jia/segmentation_cell_stats.csv"))
import Baysor as B
scale=30
min_pixels_per_cell = 15
grid_step = scale / min_pixels_per_cell
bandwidth= scale / 10
polygons = B.boundary_polygons(count_molecules, count_molecules.cell, grid_step=grid_step, bandwidth=bandwidth)
genes = count_df.gene
cells = string.(names(count_df)[2:end])
count_df = count_df[!, 2:end]
count_df = convert(SparseMatrixCSC{Int64, Int64},Matrix(count_df))
raw_count = cs.RawCountObject(count_df, cells, genes)
count_molecules.cell = string.(count_molecules.cell)
count_cells.cell = string.(count_cells.cell)
hour12 = cs.CartanaObject(count_molecules, count_cells, raw_count;
    prefix = "hour12", min_gene = 0, min_cell = 3)
hour12 = cs.normalize_object(hour12)
hour12 = cs.scale_object(hour12)
hour12 = cs.find_variable_genes(hour12; nFeatures=100)
hour12 = cs.run_pca(hour12;  method=:svd, pratio = 1, maxoutdim = 20)
hour12.metaData.cluster = string.(hour12.spmetaData.cell.cluster)
hour12.polygonData = polygons
hour12 = cs.polygons_cell_mapping(hour12)
hour12 = cs.generate_polygon_counts(hour12)
center=[14685.08, 22769.85]
hour12 = cs.compute_kidney_coordinates(hour12, center)
cs.save(hour12; filename="hour12.jld2")

### 4. IRI2d ###
count_molecules =  DataFrame(CSV.File("/mnt/sdc/new_analysis_cellscopes/seg_files/Day2_Jia/segmentation.csv"))
count_df =  DataFrame(CSV.File("/mnt/sdc/new_analysis_cellscopes/seg_files/Day2_Jia/segmentation_counts.tsv"))
count_cells =  DataFrame(CSV.File("/mnt/sdc/new_analysis_cellscopes/seg_files/Day2_Jia/segmentation_cell_stats.csv"))
import Baysor as B
scale=30
min_pixels_per_cell = 15
grid_step = scale / min_pixels_per_cell
bandwidth= scale / 10
polygons = B.boundary_polygons(count_molecules, count_molecules.cell, grid_step=grid_step, bandwidth=bandwidth)
genes = count_df.gene
cells = string.(names(count_df)[2:end])
count_df = count_df[!, 2:end]
count_df = convert(SparseMatrixCSC{Int64, Int64},Matrix(count_df))
raw_count = cs.RawCountObject(count_df, cells, genes)
count_molecules.cell = string.(count_molecules.cell)
count_cells.cell = string.(count_cells.cell)
day2 = cs.CartanaObject(count_molecules, count_cells, raw_count;
    prefix = "day2", min_gene = 0, min_cell = 3)
day2 = cs.normalize_object(day2)
day2 = cs.scale_object(day2)
day2 = cs.find_variable_genes(day2; nFeatures=100)
day2 = cs.run_pca(day2;  method=:svd, pratio = 1, maxoutdim = 20)
day2.metaData.cluster = string.(day2.spmetaData.cell.cluster)
day2.polygonData = polygons
day2 = cs.polygons_cell_mapping(day2)
day2 = cs.generate_polygon_counts(day2)
cs.save(day2; filename="day2.jld2")
data_path = "/mnt/sdc/new_analysis_cellscopes/for_imputate/IRI_2d/"
spaGE_path = "/mnt/sdc/new_analysis_cellscopes/SpaGE"
day2 = cs.run_spaGE(day2, data_path, spaGE_path)
center=[15793.59,26864.04]
day2 = cs.compute_kidney_coordinates(day2, center)
cs.save(day2; filename="day2_imputed.jld2")

### 5. IRI6w ###
count_molecules =  DataFrame(CSV.File("/mnt/sdc/new_analysis_cellscopes/seg_files/Week6_Jia/segmentation.csv"))
count_df =  DataFrame(CSV.File("/mnt/sdc/new_analysis_cellscopes/seg_files/Week6_Jia/segmentation_counts.tsv"))
count_cells =  DataFrame(CSV.File("/mnt/sdc/new_analysis_cellscopes/seg_files/Week6_Jia/segmentation_cell_stats.csv"))
import Baysor as B
scale=30
min_pixels_per_cell = 15
grid_step = scale / min_pixels_per_cell
bandwidth= scale / 10
polygons = B.boundary_polygons(count_molecules, count_molecules.cell, grid_step=grid_step, bandwidth=bandwidth)
genes = count_df.gene
cells = string.(names(count_df)[2:end])
count_df = count_df[!, 2:end]
count_df = convert(SparseMatrixCSC{Int64, Int64},Matrix(count_df))
raw_count = cs.RawCountObject(count_df, cells, genes)
count_molecules.cell = string.(count_molecules.cell)
count_cells.cell = string.(count_cells.cell)
week6 = cs.CartanaObject(count_molecules, count_cells, raw_count;
    prefix = "week6", min_gene = 0, min_cell = 3)
week6 = cs.normalize_object(week6)
week6 = cs.scale_object(week6)
week6 = cs.find_variable_genes(week6; nFeatures=100)
week6 = cs.run_pca(week6;  method=:svd, pratio = 1, maxoutdim = 20)
week6.metaData.cluster = string.(week6.spmetaData.cell.cluster)
week6.polygonData = polygons
week6 = cs.polygons_cell_mapping(week6)
week6 = cs.generate_polygon_counts(week6)
cs.save(week6; filename="week6.jld2")
data_path = "/mnt/sdc/new_analysis_cellscopes/for_imputate/IRI_6w/"
spaGE_path = "/mnt/sdc/new_analysis_cellscopes/SpaGE"
week6 = cs.run_spaGE(week6, data_path, spaGE_path)
cell_sub=cs.subset_fov(week6, [503,509,743,749], 40,40);
cell_sub=filter(:celltype => ==("CD-PC"), cell_sub);
center=[mean(cell_sub.x),mean(cell_sub.y)]
week6 = cs.compute_kidney_coordinates(week6, center)
cs.save(week6; filename="week6_imputed.jld2")



