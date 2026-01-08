library(Seurat)

obj <- readRDS("/media/singhlab/B684-19A6/CosMxJuly2022Experiment/data/seurat_processed_object/SeuratObject_sqrtNormalized_SCTransform_UMAP_annotated.RData")
obj_sub <- subset(obj, idents = c("Luminal", "Basal"))

# Clustering -------------------------------------------------------------------
obj_sub <- SCTransform(obj_sub, ncells = 15000)
# DefaultAssay(obj_sub) <- "Nanostring"
obj_sub <- NormalizeData(obj_sub)
obj_sub <- RunPCA(obj_sub, npcs = 50)
obj_sub <- FindNeighbors(obj_sub, dims = 1:50)
obj_sub <- FindClusters(obj_sub, algorithm = 4, method = "igraph")
obj_sub <- RunUMAP(obj_sub, dims = 1:50)
plotUMAP <- DimPlot(obj_sub, reduction = "umap",label = TRUE,pt.size=0.5, raster = FALSE)

# Celltype scoring -------------------------------------------------------------
# Genesets for each celltype are read into a list of dataframes, with rownames
# as gene names and a single column containing weights. The score for each
# geneset is calculated for each cell and stored in meta.data slot.

# Read marker genes
geneset <- read_xlsx(path = "/media/singhlab/B684-19A6/AACR 2024/geneset_epithelial.xlsx") %>% as.data.frame()
rownames(geneset) <- geneset$Gene
geneset$Gene <- NULL
data_df <- GetAssayData(obj_sub, assay = "Nanostring", layer = "counts")
common_genes <- intersect(rownames(geneset), rownames(data_df))
geneset <- geneset[common_genes, ]
data_df <- data_df[common_genes, ]
score_df <- as.matrix(t(data_df))%*%as.matrix(geneset) %>% as.data.frame()
obj_sub <- AddMetaData(obj_sub, metadata = score_df, col.name = paste0(colnames(geneset), "_Score"))
plt_score_cosmx_epi <- plot_prostate_features(obj_sub,
                                              variable_column = paste0(colnames(geneset), "_Score"),
                                              column_type = "score",
                                              reduction = "umap"
)

obj_ref <- readRDS(file = "/media/singhlab/B684-19A6/Single cell datasets/Prostate/integrated_dataset/Integrated_clustered_annotated_data_trim.RData")
obj_ref_sub <- subset(obj_ref, idents = c("Luminal", "Basal", "Hillock", "Club", "UK1", "UK2", "UK3", "MSMB+"))
obj_sub<- NormalizeData(obj_sub)
DefaultAssay(obj_sub) <- "Nanostring"
obj.anchors <- FindTransferAnchors(reference = obj_ref_sub, query = obj_sub, dims = 1:30,
                                   reference.reduction = "pca")
predictions <- TransferData(anchorset = obj.anchors, refdata = obj_ref_sub$CellAnnFine, dims = 1:30)
obj_sub_t <- AddMetaData(obj_sub, metadata = predictions)
Idents(obj_sub_t) <- "predicted.id"
DimPlot(obj_sub_t, reduction = "umap",label = TRUE, pt.size=0.5, raster = FALSE)
Idents(obj_sub_t) <- "seurat_clusters"
DimPlot(obj_sub_t, reduction = "umap",label = TRUE, pt.size=0.5, raster = FALSE)

obj_sub_t <- MapQuery(anchorset = obj.anchors,
                      reference = obj_ref,
                      query = obj_sub,
                      refdata = list(celltype = "CellAnnFine"),
                      reference.reduction = "pca",
                      reduction.model = "umap"
)



# Plot scores
score_expr_df <- matrix(nrow = dim(obj)[2], ncol = (length(names(common)) + 3)) %>% as.data.frame()
rownames(score_expr_df) <- colnames(obj)
colnames(score_expr_df) <- c("UMAP_1", "UMAP_2", "Cluster", names(common))
score_expr_df[, c("UMAP_1", "UMAP_2")] <- obj@reductions$umap@cell.embeddings
score_expr_df$Cluster <- obj$seurat_clusters
score_expr_df[, names(common)] <- obj@meta.data[, paste0(names(common), "_Score")]
score_expr_df$Density <- NA

colorlist <- c("#cccccc", "#00008b")
score_plots <- vector(mode = "list", length = length(names(common)))
names(score_plots) <- names(common)

x <- obj@reductions[["umap"]]@cell.embeddings[, 1]
y <- obj@reductions[["umap"]]@cell.embeddings[, 2]
for (celltype in names(common)) {

  w <- score_expr_df[, celltype] %>% as.vector(mode = "double")
  score_expr_df$Density <- plot_density(x, y, w)
  score_plots[[celltype]] <- ggplot(data = score_expr_df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = Density), size = 0.7, stroke = 0, shape = 16) +
    scale_colour_gradient2(low = "red", mid = "grey80", high = "blue") +
    labs(title = paste0(celltype, " Score"), color = "Score") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position = "none"
    )

}
design <- "
123
444
567
888
9ab
"
plt <- score_plots[[1]] + score_plots[[2]] + score_plots[[3]] +
  plot_spacer() +
  score_plots[[4]] + score_plots[[5]] + score_plots[[6]] +
  plot_spacer() +
  score_plots[[7]] + score_plots[[8]] + score_plots[[9]] +
  plot_layout(design = design,
              guides = 'collect',
              heights = c(10, 0, 10, 0, 10)
  ) &
  theme(plot.margin = unit(c(0, 0, 0, 0), units = "mm"))
plt
ggsave(filename = "Celltypes.tiff",
       plot = plt,
       device = "tiff",
       path = file.path(project_dir, "results", "figures"),
       scale = 1,
       width = 270,
       height = 270,
       units = "mm",
       dpi = 600
)
dev.off()
