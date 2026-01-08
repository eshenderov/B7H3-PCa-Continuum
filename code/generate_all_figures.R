library(decontX)
library(recount3)
library(CellChat)
library(GSVA)
library(rshannonca)
library(conos)
library(leidenAlg)
library(tiledb)
library(TileDBArray)
library(reshape2)
library(Seurat)

library(Matrix)

# mat: dgCMatrix (genes x cells)
# targets, donors: character vectors of equal length (many-to-one allowed)
combine_keep_donors <- function(mat, targets, donors) {
  if (!inherits(mat, "dgCMatrix")) mat <- as(mat, "dgCMatrix")
  rn <- rownames(mat); cn <- colnames(mat)
  stopifnot(length(targets) == length(donors), !is.null(rn), !is.null(cn))
  
  # keep only pairs present in matrix
  ok <- (targets %in% rn) & (donors %in% rn)
  targets <- targets[ok]; donors <- donors[ok]
  if (!length(targets)) return(mat)
  
  # (optional) disallow chains to avoid multi-hop additions
  if (any(targets %in% donors))
    stop("Some targets also appear as donors; overlapping/chain merges not supported.")
  
  g  <- length(rn)
  ti <- match(targets, rn)
  di <- match(donors,  rn)
  
  # Build routing matrix G (g x g):
  # - identity routes each row to itself (keeps all rows)
  # - extra entries route donor rows to their target rows
  I_i <- seq_len(g); I_j <- seq_len(g)
  G <- sparseMatrix(
    i = c(I_i, di),
    j = c(I_j, ti),
    x = 1L,
    dims = c(g, g)
  )
  
  out <- t(G) %*% mat  # sums donors into targets; donors also kept
  dimnames(out) <- list(rn, cn)
  out
}

project_dir <- "/projects/eshenderov-hpc/Shivang/projects/B7-H3 ADC Manuscript"
obj <- readRDS(file = file.path("/projects/eshenderov-hpc/Shivang/projects/scRNA-seq_atlas/Prostate/2024", "Integrated", "data", "filtered_merged_data_annotated.rds"))
source("/projects/eshenderov-hpc/Shivang/projects/B7-H3 ADC Manuscript/code/.Rprofile.R")
source("/projects/eshenderov-hpc/Shivang/GitHub/shivangsharma/scRNAseq/.Rprofile")
# Project directory ------------------------------------------------------------


# Read integrated data obj_subect --------------------------------------------------
#data_obj_loc <- read.csv(file = file.path(project_dir, "data", "data_loc.txt"), header = FALSE)[1, 1]
#obj <- readRDS(data_obj_loc)
###obj_sub <- subset(obj_sub, cells = colnames(obj_sub)[obj_sub$`doublet_status_0.08` %in% "Singlet"])
###obj_sub <- subset(obj_sub, subset = (percent.mito < 25) & (doublet_status_0.08 == "Singlet"))
meta_data <- obj@meta.data
meta_data <- meta_data[!(meta_data$Dataset %in% c("GSE172357", "GSE143791", "GSE153892")), ]
obj_sub <- subset(obj, cells = rownames(meta_data))

## 1) Replace "Normal" -> "Healthy" in a chosen metadata column (likely "Subtype")
col <- "Subtype"   # <-- change this if needed
x <- obj_sub@meta.data[[col]]
x <- as.character(x)                 # ensure character before refactoring
x[x == "Normal"] <- "Healthy"

## 2) Re-factor with the exact order you want
lvl_order <- c("Healthy", "Benign", "BPH", "Adeno", "ICC/IDC", "CRPC", "NEPC", "Small Cell")
obj_sub@meta.data[[col]] <- factor(x, levels = lvl_order)

## 3) Convert Dataset column to a factor with alphabetical levels
ds_levels <- sort(unique(obj_sub@meta.data$Dataset))
obj_sub@meta.data$Dataset <- factor(obj_sub@meta.data$Dataset, levels = ds_levels)

# Example:
targets <- c("NECTIN1","NECTIN4","KLK3")
donors  <- c("PVRL1","PVRL4","HK3")

x <- LayerData(obj_sub, assay = "RNA", layer = "counts", with.dimnames = TRUE)
x_new <- combine_keep_donors(x, targets, donors)
x_new@Dimnames <- obj_sub@assays[["RNA"]]@layers[["counts"]]@Dimnames
obj_sub@assays[["RNA"]]@layers[["counts"]] <- x_new
obj_sub <- NormalizeData(obj_sub, normalization.method = "LogNormalize", scale.factor = 1000000)

genes <- c("STEAP1", "STEAP2", "PSCA", "FOLH1", "KLK2", "TACSTD2", "NECTIN1", "NECTIN4", "DLL3", "CD276", "CD274", "PDCD1", "PDCD1LG2", "TNFRSF4", "TNFRSF9", "CTLA4", "AR", "KLK3", "TMPRSS2", "FKBP5", "IGF1R")
gene_names <- c("STEAP1", "STEAP2", "PSCA", "PSMA", "KLK2", "TROP-2", "NECTIN1", "NECTIN4", "DLL3", "B7-H3", "PD-1", "PD-L1", "PD-L2", "OX40", "4-1BB", "CTLA4", "AR", "KLK3", "TMPRSS2", "FKBP5", "IGF1R")
heatmap_row_split <- c(9, 7, 5)

hla_genes <- c("HLA-A", "HLA-B", "HLA-C", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1")
hla_gene_names <- c("HLA-A", "HLA-B", "HLA-C", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1")

#remove(obj)
gc()

obj_sub@meta.data[obj_sub@meta.data['Sample_ID'] == 'GSE181294_HP4', 'Subtype'] <- 'Benign'

# Generate Fig 2a --------------------------------------------------------------
fig2a(obj_sub, reduction = "umap_2D", project_dir = project_dir, genes = genes[1:(length(genes) - 5)], gene_names = gene_names[1:(length(genes) - 5)])
gc()

genes <- c("STEAP1", "STEAP2", "PSCA", "FOLH1", "KLK2", "TACSTD2", "NECTIN1", "NECTIN4", "DLL3", "CD276", "CD274", "PDCD1LG2", "PDCD1", "LAG3", "TNFRSF4", "TNFRSF9", "CTLA4", "AR", "KLK3", "TMPRSS2", "FKBP5", "IGF1R")
gene_names <- c("STEAP1", "STEAP2", "PSCA", "PSMA", "KLK2", "TROP-2", "NECTIN1", "NECTIN4", "DLL3", "B7-H3", "PD-1", "PD-L2", "PD-L1", "LAG3", "OX40", "4-1BB", "CTLA4", "AR", "KLK3", "TMPRSS2", "FKBP5", "IGF1R")
heatmap_row_split <- c(9, 8, 5)

# Generate Fig 2b --------------------------------------------------------------
fig2b(obj_sub, project_dir, genes = genes, gene_names = gene_names, heatmap_row_split = heatmap_row_split)
gc()

# Generate Supplementary Figure S2 ---------------------------------------------
figS2(obj_sub, project_dir)
gc()

# Generate Supplementary Figure S3 ---------------------------------------------
figS3(obj_sub, project_dir, hla_genes = hla_genes, hla_gene_names = hla_gene_names)
gc()

# Generate Supplementary Figure S4 ---------------------------------------------
figS4(obj_sub, project_dir, genes = c(genes[(length(genes) - 5):length(genes)], "KRT8", "KRT18", "KRT5", "KRT14"), gene_names = c(gene_names[(length(genes) - 5):length(genes)], "KRT8", "KRT18", "KRT5", "KRT14"))
gc()



# Generate Supplementary Figure S6 ---------------------------------------------
figS6(obj_sub, project_dir)
gc()

# Generate Supplementary Figure S7 ---------------------------------------------
figS7(obj_sub, project_dir)
gc()

# Generate Supplementary Figure S8 ---------------------------------------------
figS8(obj_sub, project_dir, genes = genes, gene_names = gene_names, heatmap_row_split = heatmap_row_split)
gc()

# Generate Supplementary Figure S9 ---------------------------------------------
figS9(obj_sub, project_dir)
gc()

# Generate Supplementary Table S2 ----------------------------------------------
tableS2(obj_sub, project_dir)

# Generate Supplementary Table S3 ----------------------------------------------
tableS3(obj_sub, project_dir)

# Generate Supplementary Table S4 ----------------------------------------------
tableS4(obj_sub, project_dir)

# Generate Supplementary Table S5 ----------------------------------------------
tableS5(obj_sub, project_dir)

