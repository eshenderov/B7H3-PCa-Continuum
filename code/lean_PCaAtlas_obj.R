library(Seurat)

# Pick reductions that actually exist
keep_dr <- intersect(c("pca","umap_2D","tsne"), names(obj@reductions))

# Pick RNA layers that actually exist in this object
rna_layers <- intersect(c("counts","data"), names(obj[["RNA"]]@layers))

obj_lean <- DietSeurat(
  object   = obj,
  assays   = "RNA",
  layers   = list(RNA = rna_layers),  # <- v5 way to keep counts/data
  dimreducs= keep_dr,
  graphs   = NULL,   # drop all graphs
  misc     = FALSE   # drop misc slot
)

saveRDS(obj_lean, "/projects/eshenderov-hpc/Shivang/projects/B7-H3 ADC Manuscript/data/lean_seurat_object.rds", compress = TRUE)
