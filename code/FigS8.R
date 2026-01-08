figS8 <- function (obj, project_dir, genes, gene_names, heatmap_row_split) {
  
  suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
    library(dplyr)
    library(tidyverse)
    library(ComplexHeatmap)
    library(gplots)
  })
  
  # ---------- Params ----------
  epithelial_types <- c("Luminal", "Neuroendocrine")
  genes <- genes  # your character vector of features
  
  # ---------- 1) Subset once & normalize once ----------
  obj_epithelial <- subset(
    obj_sub,
    subset = CoarseAnnotation %in% epithelial_types
  )
  
  # ---------- 2) Build grouping (Patient_ID_Subtype) without loops ----------
  md <- obj_epithelial@meta.data
  md$Sample_Name <- md$Sample_ID
  
  # Separate meta data of subsetted epithelial cells above -----------------------
  meta_data <- obj_epithelial@meta.data[, c("Subtype", "Dataset", "Patient_ID", "Sample_ID")]
  meta_data <- meta_data %>% unique() %>% remove_rownames()
  
  # Keep genes that exist
  genes <- intersect(genes, rownames(obj_epithelial))
  stopifnot(length(genes) > 0)
  
  # ---------- 3) Grab matrices once (sparse) ----------
  # data: normalized (RC)
  
  # counts: raw library sizes
  C <- LayerData(obj_epithelial, assay = "RNA", layer = "counts")[genes, , drop = FALSE]
  X <- LayerData(obj_epithelial, assay = "RNA", layer = "data")[genes, , drop = FALSE]
  # obj_epithelial2 <- NormalizeData(obj_epithelial, normalization.method = "RC", scale.factor = 1000000)
  # C <- GetAssayData(obj_epithelial, assay = "RNA", layer = "counts")
  # X <- t(t(C)/colSums2(C))*1e6
  # X_RC <- X
  # X_LN <- X %>% log1p()
  # X <- X_RC[genes, , drop = FALSE]
  # X <- X_LN[genes, , drop = FALSE]
  
  # ---------- 4) Make a sparse one-hot group matrix G (cells x groups) ----------
  grp <- factor(md$Sample_Name)
  G <- Matrix::Matrix(model.matrix(~ grp - 1), sparse = TRUE)     # columns in same order as levels(grp)
  colnames(G) <- levels(grp)
  n_cells_by_grp <- Matrix::colSums(G)                            # cell_count
  
  # ---------- 5) Per-group gene-wise mean & sd using sums and sums-of-squares ----------
  # Sums by group (genes x groups)
  S  <- X %*% G                               # sum(x)
  # Sums of squares (modify a copy's @x then multiply)
  X2 <- X
  slot(X2, "x") <- slot(X2, "x")^2   # same as: X2@x <- X2@x^2
  S2 <- X2 %*% G               # sum(x^2) without duplicating memory
  # (The trick above reuses X sparsity; if you prefer clarity, do: X2 <- X; X2@x <- X2@x^2; S2 <- X2 %*% G)
  
  # Broadcast n across genes
  N  <- matrix(rep(n_cells_by_grp, each = nrow(S)), nrow = nrow(S))
  
  # Means and (unbiased) variance
  mu  <- S / pmax(N, 1)
  var <- (S2 - (S^2) / pmax(N, 1)) / pmax(N - 1, 1)
  sdg <- sqrt(pmax(var, 0))
  
  # Coefficient of variation (%) per gene x group
  cv_pct <- 100 * (sdg / mu)
  cv_pct[!is.finite(cv_pct)] <- NA_real_
  
  # Arrange outputs
  dispersion_df <- as.data.frame(cv_pct)
  rownames(dispersion_df) <- genes
  colnames(dispersion_df) <- colnames(G)  # Sample_Name columns
  
  # ---------- 6) cell_count & median reads per cell (from counts) ----------
  # Per-cell library sizes (use counts for "reads")
  lib_sizes <- Matrix::colSums(C)  # length = number of cells
  
  # Median per group — fast, no per-sample Seurat subset
  # (median isn't linear, so we use tapply on the vector)
  median_reads_per_cell <- tapply(
    lib_sizes,
    grp,
    median,
    na.rm = TRUE
  ) %>% as.numeric()
  names(median_reads_per_cell) <- levels(grp)
  
  cell_count <- as.numeric(n_cells_by_grp)
  names(cell_count) <- levels(grp)
  
  # Quick peek
  cat(sprintf("genes x samples in dispersion: %d x %d\n",
              nrow(dispersion_df), ncol(dispersion_df)))
  
  # ----------------------- New method ends ------------------------------------
  # ----------------------------------------------------------------------------
  
  # Rename subtypes --------------------------------------------------------------
  meta_data$Subtype <- factor(meta_data$Subtype, levels = c("Healthy", "Benign", "Adeno", "ICC/IDC", "CRPC", "NEPC", "Small Cell") %>% rev())
  
  # Reorder data -----------------------------------------------------------------
  order <- c(which(meta_data$Subtype %in% c("Healthy")), 
             which(meta_data$Subtype %in% c("Benign")), 
             which(meta_data$Subtype %in% c("Adeno")), 
             which(meta_data$Subtype %in% c("ICC/IDC")), 
             which(meta_data$Subtype %in% "CRPC"), 
             which(meta_data$Subtype %in% "NEPC"), 
             which(meta_data$Subtype %in% "Small Cell")
  ) %>% rev()
  meta_data <- meta_data[order, ]
  dispersion_df <- dispersion_df[, meta_data$Sample_ID]
  
  luminal_basal_genes <- genes

  
  # Prepare data for heatmap -----------------------------------------------------
  anno <- meta_data %>% dplyr::select(-c("Dataset"))
  dispersion_df <- dispersion_df[, meta_data$`Sample_ID`]
  rownames(anno) <- anno$`Sample_ID`
  # AR_Score <- AR_Score[anno$`Sample Name`]
  # NE_Score <- NE_Score[anno$`Sample Name`]
  
  # Assign colors to subtypes of PCa, AR and NE Score in heatmap -----------------
  # colors <- c("#5E4FA2", "#5E4FA2", "#3288BD", "#ABDDA4", "#E6F598", "#E6F598", "#E6F598", "#E6F598", "#FEE08B", "#FEE08B", "#F46D43", "#D53E4F", "#9E0142")
  # names(colors) <- c("Healthy", "Healthy", "Benign", "Tumor", "Tumor", "Tumor", "Tumor", "ICC/IDC", "ICC/IDC", "CRPC", "mixed-NEPC", "NEPC", "Small Cell")
  colors <- c("#5E4FA2", "#3288BD", "#ABDDA4", "#E6F598", "#FEE08B", "#F46D43", "#9E0142")
  names(colors) <- c("Healthy", "Benign", "Adeno", "ICC/IDC", "CRPC", "NEPC", "Small Cell")
  anno_colors_col <- list(Subtype = rev(colors))
  
  # Create heatmap column section ------------------------------------------------
  col_ha <- HeatmapAnnotation(Subtype = anno$Subtype,
                              col = anno_colors_col,
                              annotation_label = c("Subtype"),
                              annotation_name_side = "right",
                              border = c(Subtype = FALSE),
                              simple_anno_size = unit(5, "mm")
  )
  
  # Assign colors to rows of genes based on their category -----------------------
  anno_colors_row <- list(Gene = c("#FF0000", "#49FF00", "#0092FF"))
  names(anno_colors_row$Gene) <- c("ADC", "Immune", "AR reg.")
  
  # Create heatmap row section ---------------------------------------------------
  row_ha <- HeatmapAnnotation(Gene = rep(x = c("ADC", "Immune", "AR reg."), times = heatmap_row_split) %>% 
                                factor(levels = c("ADC", "Immune", "AR reg.")),
                              col = anno_colors_row,
                              which = "row",
                              show_annotation_name = FALSE,
                              simple_anno_size = unit(2, "mm")
                              )
  
  # Create heatmap ---------------------------------------------------------------
  ht_plot <- Heatmap(matrix = dispersion_df %>% as.matrix() %>% log1p(),
                     name = "Dispersion",
                     col = colorpanel(n = 500,
                                      low = "#440154FF",
                                      mid = "#21908CFF",
                                      high = "#FDE725FF"),
                     row_split = rep(x = c("ADC", "Immune", "AR reg."), times = heatmap_row_split) %>% factor(levels = c("ADC", "Immune", "AR reg.")),
                     cluster_columns = FALSE,
                     cluster_rows = FALSE,
                     row_gap = unit(x = c(1, 2), units = "mm"),
                     border = TRUE,
                     left_annotation = row_ha,
                     top_annotation = col_ha,
                     row_labels = gene_names,
                     show_column_names = FALSE,
                     row_title = NULL,
                     row_names_side = "left",
                     width = (ncol(dispersion_df) * unit(2, "mm")) * 1.2,
                     height = (nrow(dispersion_df) + 3) * unit(3.5, "mm") * 1.2,
                     row_names_gp = gpar(fontsize = 10)
  )
  
  # Save heatmap into pdf --------------------------------------------------------
  pdf(file = file.path(project_dir, "figures", 'FigS8.pdf'), 
      width = convertUnit(x = (dim(dispersion_df)[2] * unit(2, "mm")) * 1.7, 
                          unitTo = "inches"), 
      height = convertUnit(x = (dim(dispersion_df)[1] + 3) * unit(3.5, "mm") * 1.7, 
                           unitTo = "inches"))
  ComplexHeatmap::draw(ht_plot)
  dev.off()
  
  tiff(file = file.path(project_dir, "figures", 'FigS8.tiff'), 
      width = convertUnit(x = (dim(dispersion_df)[2] * unit(2, "mm")) * 102, 
                          unitTo = "inches"), 
      height = convertUnit(x = (dim(dispersion_df)[2] + 3) * unit(3.5, "mm") * 17, 
                           unitTo = "inches"))
  ComplexHeatmap::draw(ht_plot)
  dev.off()
  
}