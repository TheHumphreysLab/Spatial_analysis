library(Seurat)
library(plot1cell)
immune_sub <- CreateSeuratObject(counts = immune_counts, project = "immune_sub", min.cells = 0, min.features = 0)
immune_sub <- SCTransform(immune_sub,  verbose = FALSE, variable.features.n = 100)
immune_sub <- RunPCA(immune_sub, verbose = FALSE)
immune_sub <- RunUMAP(immune_sub, dims = 1:8, verbose = FALSE, min.dist = 0.2)
immune_sub <- FindNeighbors(immune_sub, dims = 1:8, verbose = FALSE)
immune_sub <- FindClusters(immune_sub, verbose = FALSE, resolution = 0.1)
DefaultAssay(immune_sub)<-"RNA"
immune_sub<-NormalizeData(immune_sub)
DimPlot(immune_sub, label = TRUE) + NoLegend()
immune_sub@meta.data$celltype=as.character(immune_sub@meta.data$seurat_clusters)
immune_sub@meta.data$celltype[immune_sub@meta.data$celltype %in% c(0, 2, 3)]<-"Macrophage"
immune_sub@meta.data$celltype[immune_sub@meta.data$celltype %in% c(4)]<-"Tcell"
immune_sub@meta.data$celltype[immune_sub@meta.data$celltype %in% c(1)]<-"Ly6e+ lymphocyte"
immune_sub@meta.data$celltype[immune_sub@meta.data$celltype %in% c(5)]<-"C3+ cell"
immune_sub@meta.data$celltype[immune_sub@meta.data$celltype %in% c(6)]<-"Cxcl10+ cell"
immune_sub@meta.data$celltype[immune_sub@meta.data$celltype %in% c(7)]<-"Plasma cell"
immune_sub<-SetIdent(immune_sub, value = "celltype")
save(immune_sub, file = "immune_subclustering.Rda")
levels(immune_sub)<-c("Macrophage","Ly6e+ lymphocyte","Tcell","Cxcl10+ cell","C3+ cell","Plasma cell")

tiff(filename = "immune_clustering.tiff", width = 5, height = 5, res = 300, units = "in")
DimPlot(immune_sub, label = TRUE, label.size = 5, repel = T) + NoLegend()+ NoAxes() + 
  scale_color_manual(values = c("seagreen2", "magenta","dodgerblue","firebrick1","sienna1","darkgreen"))
dev.off()
levels(immune_sub)<-rev(levels(immune_sub))
tiff(filename = "immune_markers.tiff", width = 5, height = 4, res = 300, units = "in")
DotPlot(immune_sub, features = c("Cd74","C1qa","Ly6e","Cd3e","Cxcl10","C3","Igkc"), dot.scale = 8, cols = c("white","red"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        axis.line = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.key = element_rect(colour = NA, fill = NA),
        axis.text = element_text(size = 12),axis.title=element_text(size=8),legend.text=element_text(size=10), 
        legend.title = element_text(size = 10),legend.position="right")+ylab("")+xlab("")+
  guides(size=guide_legend(title="% Cells"),colour = guide_colorbar(title="Avg Expr"))
dev.off()

##### plot annotation
coordinates<-read.csv("../spatial/juanru/juanrufinal/cartana/Baysorresult/Week6_Jia/segmentation_cell_stats.csv")
load("immune_subclustering.Rda")
immune_cells<-colnames(immune_sub)
metadata<-immune_sub@meta.data
coordinates$cell2<-paste("Week6", coordinates$cell, sep = "_")
immune_coord<-coordinates[coordinates$cell2 %in% immune_cells,]
immune_coord$celltype<-metadata$celltype
non_immune<-coordinates[coordinates$cell2 %in% setdiff(coordinates$cell2, immune_cells),]
levels(immune_sub)<-c("Macrophage","Ly6e+ lymphocyte","Tcell","Cxcl10+ cell","C3+ cell","Plasma cell")
immune_coord$celltype<-factor(immune_coord$celltype, levels = levels(immune_sub))
tiff(filename = "spatial_anno.tiff", width = 6, height = 5, res = 300, units = "in")
ggplot()+geom_point(data=non_immune, aes(x, y), size=0.05, color="gray90")+
  geom_point(data=immune_coord, aes(x, y, color=celltype), size=0.1)+
  theme(panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_blank(), axis.ticks = element_blank(),
        legend.key = element_rect(colour = NA, fill = NA),
        axis.text = element_blank(),legend.text=element_text(size=10), 
        legend.title = element_blank(),legend.position="right")+ylab("")+xlab("")+
  scale_color_manual(values = c("seagreen2", "magenta","dodgerblue","firebrick1","sienna1","darkgreen"))+
  guides(color=guide_legend(override.aes = list(size = 8)))
dev.off()

tiff(filename = "spatial_immune.tiff", width = 4.5, height = 5, res = 300, units = "in")
ggplot()+geom_point(data=non_immune, aes(x, y), size=0.05, color="gray90")+
  geom_point(data=immune_coord, aes(x, y), size=0.1, color="black")+
  theme(panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_blank(), axis.ticks = element_blank(),
        legend.key = element_rect(colour = NA, fill = NA),
        axis.text = element_blank(),legend.text=element_text(size=10), 
        legend.title = element_blank(),legend.position="right")+ylab("")+xlab("")+
  guides(color=guide_legend(override.aes = list(size = 8)))
dev.off()


tiff(filename = "immune_subtype_separate.tiff", width = 10, height = 8, res = 300, units = "in")
ggplot()+geom_point(data=non_immune, aes(x, y), size=0.05, color="gray90")+
  geom_point(data=immune_coord, aes(x, y, color=celltype), size=0.2)+
  theme(panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_blank(), axis.ticks = element_blank(),
        legend.key = element_rect(colour = NA, fill = NA),
        axis.text = element_blank(),legend.text=element_text(size=10), 
        legend.title = element_blank(),legend.position="none", 
        strip.background=element_blank(), strip.text = element_text(size=20))+ylab("")+xlab("")+
  scale_color_manual(values = c("seagreen2", "magenta","dodgerblue","firebrick1","sienna1","darkgreen"))+
  facet_wrap(~celltype)
dev.off()

#### analyzing the week6 snRNA-seq data -- Kirita etal 2020 PNAS
week6<-read.table("MouseIRI_6weeks.dge.txt.gz")
week6<-Matrix(as.matrix(week6), sparse = T)
immune_sub <- CreateSeuratObject(counts = week6, project = "week6", min.cells = 0, min.features = 0)
ifnb.list <- SplitObject(immune_sub, split.by = "orig.ident")
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = ifnb.list)
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
immune.combined <- IntegrateData(anchorset = immune.anchors)
DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

DimPlot(immune.combined, label = TRUE) + NoLegend()

DefaultAssay(immune.combined)<-"RNA"
immune.combined<-NormalizeData(immune.combined)
FeaturePlot(immune.combined, "Themis", order=T)

scRNA_immune<-subset(immune.combined, idents=8)
DefaultAssay(scRNA_immune)<-"integrated"

immune_sub1 <- CreateSeuratObject(counts = immune_uuo@assays$RNA@counts, project = "uuo", min.cells = 0, min.features = 0)
ifnb.list <- SplitObject(immune_sub1, split.by = "orig.ident")
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 500)
})
features <- SelectIntegrationFeatures(object.list = ifnb.list)
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
scRNA_immune <- IntegrateData(anchorset = immune.anchors)

DefaultAssay(scRNA_immune) <- "integrated"
scRNA_immune <- ScaleData(scRNA_immune, verbose = FALSE)

scRNA_immune <- RunPCA(scRNA_immune, verbose = FALSE)
scRNA_immune <- RunUMAP(scRNA_immune, reduction = "pca", dims = 1:8, min.dist = 0.5)
scRNA_immune <- FindNeighbors(scRNA_immune, reduction = "pca", dims = 1:8)
scRNA_immune <- FindClusters(scRNA_immune, resolution = 2)
DimPlot(scRNA_immune, label = TRUE) + NoLegend()
scRNA_immune@meta.data$celltype<-as.character(scRNA_immune@meta.data$seurat_clusters)
scRNA_immune@meta.data$celltype[scRNA_immune@meta.data$celltype %in% c(11)]<-"Mmp12+ Mac"
scRNA_immune@meta.data$celltype[scRNA_immune@meta.data$celltype %in% c(9)]<-"Proliferating"
scRNA_immune@meta.data$celltype[scRNA_immune@meta.data$celltype %in% c(17)]<-"cDC1"
scRNA_immune@meta.data$celltype[scRNA_immune@meta.data$celltype %in% c(3)]<-"cDC2"
scRNA_immune@meta.data$celltype[scRNA_immune@meta.data$celltype %in% c(13)]<-"Ccr7+ DC"
scRNA_immune@meta.data$celltype[scRNA_immune@meta.data$celltype %in% c(8)]<-"Arg1+ Mac"
scRNA_immune@meta.data$celltype[scRNA_immune@meta.data$celltype %in% c(10,12)]<-"Ly6c2+ Mono"
scRNA_immune@meta.data$celltype[scRNA_immune@meta.data$celltype %in% c(5)]<-"IFN induced Mac"
scRNA_immune@meta.data$celltype[scRNA_immune@meta.data$celltype %in% c(18,20)]<-"Plasma cell"
scRNA_immune@meta.data$celltype[scRNA_immune@meta.data$celltype %in% c(0,1,6)]<-"Ccr2+ Mac"
scRNA_immune@meta.data$celltype[scRNA_immune@meta.data$celltype %in% c(2)]<-"Egr1+ Mac"
scRNA_immune@meta.data$celltype[scRNA_immune@meta.data$celltype %in% c(14)]<-"Neutrophil"
scRNA_immune@meta.data$celltype[scRNA_immune@meta.data$celltype %in% c(15,16)]<-"NK cell"
scRNA_immune@meta.data$celltype[scRNA_immune@meta.data$celltype %in% c(4,7,19)]<-"T cell"

scRNA_immune<-SetIdent(scRNA_immune, value = "celltype")

tiff(filename = "immune_uuo.tiff", width = 6, height = 5, res = 300, units = "in")
DimPlot(scRNA_immune, label = TRUE, repel = T) + NoLegend()+NoAxes()
dev.off()
levels(scRNA_immune)<-c("Mmp12+ Mac", "Ccr2+ Mac", "Egr1+ Mac", "Arg1+ Mac", "Ly6c2+ Mono" ,
              "NK cell", "T cell","IFN induced Mac","Plasma cell","Proliferating"   ,"cDC1",
              "Ccr7+ DC"   ,   "cDC2",  "Neutrophil"  )
scRNA_immune<-ScaleData(scRNA_immune)
sc_immune<-AverageExpression(scRNA_immune, assays = "RNA", slot = "scale.data")
sc_immune<-sc_immune$RNA
sc_immune<-data.frame(sc_immune)
immune_sub<-ScaleData(immune_sub)

sp_immune<-AverageExpression(immune_sub, assays = "RNA", slot = "scale.data")
sp_immune<-sp_immune$RNA
sp_immune<-data.frame(sp_immune)
p=run_correlation(sc_immune, sp_immune, ngenes = 100, 
                color.use = colorRampPalette(c("white","lemonchiffon1", "red", "darkred"))(255))

tiff(filename = "correlation_score.tiff", width = 3, height = 4, res = 300, units = "in")
p
dev.off()

