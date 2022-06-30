## These are the R functions I created for Cartana data analysis. 
## The goal is to create a R workspace to seamlessly process the results from scanpy
## so I can enhance the data analysis with R.

run_seurat<-function(stat_file, count_file, sample_name, min_gene){
  cartana_dat_stat<-read.csv(stat_file)
  cartana_dat_stat<-cartana_dat_stat[cartana_dat_stat$n_transcripts>min_gene,]
  cartana_dat<-read.table(count_file, header = T)
  rownames(cartana_dat)<-cartana_dat$gene
  cartana_dat<-cartana_dat[,-1]
  cartana_dat<-cartana_dat[, cartana_dat_stat$cell]
  colnames(cartana_dat)<-paste(sample_name,colnames(cartana_dat), sep = "_")
  pbmc <- CreateSeuratObject(counts = cartana_dat, project = sample_name, 
                             min.cells = 0, min.features = 0)
  pbmc <- SCTransform(pbmc)
  pbmc <- RunPCA(pbmc)
  pbmc <- FindNeighbors(pbmc, dims = 1:10)
  pbmc <- FindClusters(pbmc, resolution = 2)
  pbmc <- RunUMAP(pbmc, dims = 1:10, min.dist = 0.2)
  pbmc
}


convert_scanpy<-function(stat_file, count_file, sample_name, min_gene, 
                         umap_df, cluster_df){
  umap1<-read.csv(umap_df)
  colnames(umap1)<-c("UMAP_1","UMAP_2")
  cartana_dat_stat<-read.csv(stat_file)
  cartana_dat_stat<-cartana_dat_stat[cartana_dat_stat$n_transcripts>min_gene,]
  cartana_dat<-read.table(count_file, header = T)
  rownames(cartana_dat)<-cartana_dat$gene
  cartana_dat<-cartana_dat[,-1]
  cartana_dat<-cartana_dat[, cartana_dat_stat$cell]
  colnames(cartana_dat)<-paste(sample_name,colnames(cartana_dat), sep = "_")
  pbmc <- CreateSeuratObject(counts = cartana_dat, project = sample_name, 
                             min.cells = 0, min.features = 0)
  pbmc <- NormalizeData(pbmc)
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 50)
  all.genes <- rownames(pbmc)
  pbmc <- ScaleData(pbmc, features = all.genes)
  pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
  pbmc <- FindNeighbors(pbmc, dims = 1:10)
  pbmc <- FindClusters(pbmc, resolution = 0.5)
  pbmc <- RunUMAP(pbmc, dims = 1:10)
  rownames(umap1)<-colnames(pbmc)
  umap1<-as.matrix(umap1)
  pbmc@reductions$umap@cell.embeddings<-umap1
  clusters<-read.csv(cluster_df)
  clusters<-clusters$cluster
  pbmc@meta.data$cluster<-clusters
  pbmc<-SetIdent(pbmc, value = "cluster")
  pbmc@reductions$spatial<-pbmc@reductions$umap
  pbmc@reductions$spatial@key<-"Spatial_"
  spatial_coord<-cartana_dat_stat[,c("x","y")]
  rownames(spatial_coord)<-colnames(pbmc)
  colnames(spatial_coord)<-c("Spatial_1","Spatial_2")
  spatial_coord<-as.matrix(spatial_coord)
  pbmc@reductions$spatial@cell.embeddings<-spatial_coord
  pbmc
}

make_mock_object <- function(object_full, molecue_count, seg_file, append){
vizgen.obj <- LoadVizgen(data.dir = "/mnt/sdc/spatial/s2r1/", fov = "s2r1")
vizgen.obj@meta.data<-object_full@meta.data
colnames(vizgen.obj@meta.data)[2:3]<-c("nCount_Vizgen",   "nFeature_Vizgen")
vizgen.obj@assays$Vizgen@counts<-object_full@assays$RNA@counts
vizgen.obj@assays$Vizgen@data<-object_full@assays$RNA@data
cell_umap<-object_full@reductions$umap@cell.embeddings
colnames(cell_umap)<-c("x","y")
rownames(cell_umap)<-NULL
aa<-CreateCentroids(data.frame(object_full@reductions$umap@cell.embeddings))
vizgen.obj[["s2r1"]]@boundaries$centroids <- aa
aa<-read.csv(seg_file)
aa$cell<-paste0(append, aa$cell)
aa<-aa[aa$cell %in% colnames(vizgen.obj),]
bb <- CreateSegmentation(aa)
vizgen.obj[["s2r1"]]@boundaries$segmentation<-bb
aa <- read.csv(molecue_count)
colnames(aa)[2:3]<-c("x","y")
aa<-aa[,c(2,3,1)]
#aa<-add_fov(aa, n_fields_x = 10, n_fields_y = 10)
aa<-CreateMolecules(aa)
vizgen.obj[["s2r1"]]@molecules$molecules<-aa
vizgen.obj <- SCTransform(vizgen.obj, assay = "Vizgen", clip.range = c(-10, 10), )
vizgen.obj <- RunPCA(vizgen.obj, npcs = 30, features = rownames(vizgen.obj))
vizgen.obj <- RunUMAP(vizgen.obj, dims = 1:30)
vizgen.obj <- FindNeighbors(vizgen.obj, reduction = "pca", dims = 1:30)
vizgen.obj <- FindClusters(vizgen.obj, resolution = 0.3)
vizgen.obj@active.ident<-object_full@active.ident
vizgen.obj
}

plot_aoi <- function( vizgen.obj, 
                      molecule_df,
                      region, 
                      x_seg, 
                      y_seg, 
                      genes=c("Itga8","Nphs2","Ehd3"), 
                      mol_col=c("green", "yellow", "red")
                      ){
all_limits<-split_field(count_df = molecule_df, n_fields_x = x_seg, n_fields_y = y_seg)
if(length(region)==1){
  select_limit<-all_limits[[region]]
  cropped.coords <- Crop(vizgen.obj[["s2r1"]], x = c(select_limit[1], select_limit[2]),
                         y = c(select_limit[3], select_limit[4]), coords = "plot")
} else {
  select_limit<-matrix(data = NA, ncol = 4, nrow = length(region))
  for(i in 1:length(region)){
    select_limit[i,]<-all_limits[[region[i]]]
  }
  select_limit<-data.frame(select_limit)
  colnames(select_limit)<-c("x1","x2","y1","y2")
  cropped.coords <- Crop(vizgen.obj[["s2r1"]], x = c(min(c(select_limit$y1,select_limit$y2)), max(c(select_limit$y1,select_limit$y2))), 
                         y = c(min(c(select_limit$x1,select_limit$x2)), max(c(select_limit$x1,select_limit$x2))), coords = "plot")
}
vizgen.obj[["hippo"]] <- cropped.coords
vizgen.obj[["hippo"]][["simplified.segmentations"]] <- Simplify(coords = vizgen.obj[["hippo"]][["segmentation"]],
                                                                tol = 3)
DefaultBoundary(vizgen.obj[["hippo"]]) <- "simplified.segmentations"
ImageDimPlot(vizgen.obj, fov = "hippo", molecules = genes, cols = "polychrome",
             mols.size = 1, alpha = 0.6, mols.cols = mol_col,border.size = 0.1, border.color = 'gray50')
}

find_fov<-function(x1, x2, y1, y2, x, y){
  if(x >=x1 & x < x2 & y >= y1 & y < y2){
    TRUE
  } else {
    FALSE
  }
}

split_field <- function(count_df,n_fields_x, n_fields_y){
  x_min <- min(count_df$x)-1
  x_max <- max(count_df$x)+1
  y_min <- min(count_df$y)-1
  y_max <- max(count_df$y)+1
  x_seg<-(x_max-x_min)/n_fields_x
  y_seg<-(y_max-y_min)/n_fields_y
  x_collect<-c()
  x_collect[1]<-x_min
  for(i in 1:n_fields_x){
    x_collect[i+1]<- x_min + x_seg * i
  }
  y_collect<-c()
  y_collect[1]<-y_min
  for(i in 1:n_fields_y){
    y_collect[i+1]<- y_min + y_seg * i
  }
  x_collect2<-list()
  for(i in 1:n_fields_x){
    x_collect2[[i]] <-c(x_collect[i], x_collect[i+1])
  }
  y_collect2<-list()
  for(i in 1:n_fields_y){
    y_collect2[[i]] <-c(y_collect[i], y_collect[i+1])
  }
  coord_limits<-list()
  coord_limits2<-c()
  for(i in 1:length(x_collect2)){
    aa<-x_collect2[[i]]
    for(j in 1:length(y_collect2)){
      coord_limits[[j]]<-c(aa, y_collect2[[j]])
    }
    coord_limits2<-c(coord_limits2,coord_limits)
  }
  names(coord_limits2)<-paste0("fov",1:(n_fields_x*n_fields_y))
  coord_limits2
}

add_fov <- function(count_df, n_fields_x, n_fields_y){
  coord_limits<-split_field(count_df, n_fields_x, n_fields_y)
  fov_col <-c()
  pb <- progress_bar$new(
    format = "  Running [:bar] :percent eta: :eta",
    clear = FALSE, total = nrow(count_df), width = 100)
 for(i in 1:nrow(count_df)){
   x0=count_df$x[i]
   y0=count_df$y[i]
   fov_assign<-c()
   for(j in 1:length(coord_limits)){
     lim_xy<-coord_limits[[j]]
     fov_assign[j]<-find_fov(x1=lim_xy[1], x2=lim_xy[2],y1=lim_xy[3],y2=lim_xy[4], x=x0, y=y0)
    }
   fov_col[i]<-names(coord_limits)[fov_assign]
   pb$tick()
   Sys.sleep(1 / nrow(count_df))
 }
  count_df$fov<-fov_col
  count_df
}

plot_rect<-function(count_data, x_seg, y_seg, gene_id){
rect_num <- x_seg * y_seg
cc<-matrix(data = NA, ncol = 4, nrow = rect_num)
coord_limits<-split_field(count_data, x_seg, y_seg)
for(i in 1:rect_num){
  cc[i,]<-coord_limits[[i]]
}
cc<-data.frame(cc)
cc$fov<-paste0("fov",1:rect_num)
ggplot() + 
  geom_point(data = count_data, aes(x, y), size=0.01, color="gray80")+
  geom_point(data = count_data[count_data$gene==gene_id,], aes(x, y), size=0.01, color="red")+
  geom_rect(data=cc, 
            mapping=aes(xmin=X1, 
                        xmax=X2, ymin=X3, 
                        ymax=X4), color="black", fill=NA, 
            alpha=0.5) +
  geom_text(data=cc, aes(x=X1+(X2-X1)/2, y=X3+(X4-X3)/2, label=fov), size=4, color="black")
}

