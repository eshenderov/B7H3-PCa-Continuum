# ---- Dependencies -------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(matrixStats)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(limma)
  library(readxl)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# ---- Small helpers ------------------------------------------------------------
df_transpose <- function(x) {
  as.data.frame(t(as.data.frame(x, check.names = FALSE)), check.names = FALSE)
}

RLE_plot <- function(df, main = "RLE (median-centered)") {
  # Simple, fast RLE-like plot: subtract per-sample median, boxplot the deltas
  z <- as.matrix(df)
  z <- sweep(z, 2, matrixStats::colMedians(z, na.rm = TRUE), FUN = "-")
  oldpar <- par(no.readonly = TRUE); on.exit(par(oldpar), add = TRUE)
  par(mar = c(6, 4, 2, 1))
  boxplot(as.data.frame(z), outline = FALSE, las = 2, main = main, ylab = "Δ from median")
  invisible(NULL)
}

dispersion_safe <- function(x, eps = 1e-8) {
  x <- x[!is.na(x)]
  if (length(x) < 2) return(NA_real_)
  v <- stats::var(x)
  m <- mean(x)
  if (abs(m) < eps) return(NA_real_)
  v / m
}

# ---- Main function ------------------------------------------------------------
# Arguments:
#   obj:                 Seurat object
#   project_dir:         base project directory (plots saved to project_dir/figures)
#   genes:               character vector (order defines row order & group splits)
#   gene_names:          character vector of pretty labels (same length/order as genes)
#   heatmap_row_split:   integer vector of length 3, e.g. c(n_ADC, n_Immune, n_AR)
make_dispersion_heatmap <- function(obj, project_dir, genes, gene_names, heatmap_row_split) {
  
  stopifnot(length(genes) == length(gene_names))
  stopifnot(length(heatmap_row_split) == 3)
  stopifnot(sum(heatmap_row_split) == length(genes))
  
  # ---- Part 1: Subset, pseudobulk, batch-correct -----------------------------
  obj_epithelial <- subset(
    obj,
    subset = CoarseAnnotation %in% c("Luminal", "Neuroendocrine")
  )
  #obj_epithelial <- NormalizeData(obj_epithelial, normalization.method = "RC")
  data_epithelial <- LayerData(obj_epithelial, assay = "RNA", layer = "counts")
  n_cells_epithelial <- ncol(obj_epithelial)
  
  meta_data <- obj_epithelial@meta.data[, c("Subtype", "Dataset", "Patient_ID")]
  meta_data <- meta_data %>% distinct() %>% tibble::remove_rownames()
  meta_data$`Sample Name` <- paste0(meta_data$Patient_ID, "_", meta_data$Subtype)
  
  data <- data.frame(
    matrix(nrow = nrow(obj_epithelial), ncol = nrow(meta_data)),
    row.names = rownames(obj_epithelial),
    check.names = FALSE
  )
  colnames(data) <- meta_data$`Sample Name`
  
  obj_meta_data <- obj_epithelial@meta.data[, c("Subtype", "Dataset", "Patient_ID")]
  
  for (i in seq_len(nrow(meta_data))) {
    subtype     <- meta_data$Subtype[i]
    patient_id  <- meta_data$Patient_ID[i]
    sample_name <- paste0(patient_id, "_", subtype)
    which_cells <- which((obj_meta_data$Patient_ID %in% patient_id) &
                           (obj_meta_data$Subtype %in% subtype))
    data[, sample_name] <- rowSums2(as.matrix(data_epithelial[, which_cells, drop = FALSE]))
  }
  
  # Normalize to relative counts (per-billion scale, then transpose back)
  data <- df_transpose(df_transpose(data) / colSums(data) * 1e9)
  
  # Rename/relevel Subtype factor for consistent ordering (rev of listed order)
  meta_data$Subtype <- factor(
    meta_data$Subtype,
    levels = rev(c("Healthy", "Benign", "Adeno", "ICC/IDC", "CRPC", "NEPC", "Small Cell"))
  )
  
  # Reorder columns by subtype groups (Healthy -> ... -> Small Cell)
  order_idx <- c(
    which(meta_data$Subtype %in% "Healthy"),
    which(meta_data$Subtype %in% "Benign"),
    which(meta_data$Subtype %in% "Adeno"),
    which(meta_data$Subtype %in% "ICC/IDC"),
    which(meta_data$Subtype %in% "CRPC"),
    which(meta_data$Subtype %in% "NEPC"),
    which(meta_data$Subtype %in% "Small Cell")
  )
  # keep only those present; then reverse to Small Cell ... Healthy (to match example)
  order_idx <- rev(order_idx)
  meta_data <- meta_data[order_idx, , drop = FALSE]
  data <- data[, meta_data$`Sample Name`, drop = FALSE]
  
  # Batch correction on log1p scale
  vsd <- log1p(as.matrix(data[, meta_data$`Sample Name`, drop = FALSE]))
  vsd_df_uncorrected <- as.data.frame(vsd)
  # (Optional diagnostic)
  try(RLE_plot(vsd_df_uncorrected, main = "RLE before removeBatchEffect"), silent = TRUE)
  
  vsd_df_corrected <- limma::removeBatchEffect(
    vsd,
    batch = meta_data$Dataset,
    group = meta_data$Subtype
  ) %>% as.data.frame()
  
  try(RLE_plot(vsd_df_corrected, main = "RLE after removeBatchEffect"), silent = TRUE)
  vsd <- vsd_df_corrected
  
  # (Optional) AR/NE scores if the Excel file exists (not used downstream)
  try({
    geneset_path <- file.path(
      "/projects/eshenderov-hpc/Shivang/projects/scRNA-seq_atlas/Prostate/integrated_dataset",
      "Genesets_Emir.xlsx"
    )
    if (file.exists(geneset_path)) {
      features <- read_xlsx(path = geneset_path) %>% as.list()
      features$`AR Score` <- features$`AR Score`[!is.na(features$`AR Score`)]
      features$`NE Score` <- features$`NE Score`[!is.na(features$`NE Score`)]
      AR_Score <- (colSums2(as.matrix(vsd[features$`AR Score`, , drop = FALSE])) /
                     length(features$`AR Score`))
      NE_Score <- (colSums2(as.matrix(vsd[features$`NE Score`, , drop = FALSE])) /
                     length(features$`NE Score`))
    }
  }, silent = TRUE)
  
  # Subset expression matrix to requested genes
  # (Only keep genes that exist in vsd; warn for any missing)
  found_genes <- intersect(genes, rownames(vsd))
  missing_genes <- setdiff(genes, found_genes)
  if (length(missing_genes) > 0) {
    warning("These genes were not found and will be dropped: ",
            paste(missing_genes, collapse = ", "))
  }
  vsd_sub <- vsd[found_genes, , drop = FALSE]
  
  # Annotation dataframe (drop Dataset)
  anno <- meta_data %>% dplyr::select(-any_of("Dataset"))
  vsd_sub <- vsd_sub[, anno$`Sample Name`, drop = FALSE]
  rownames(anno) <- anno$`Sample Name`
  
  # ---- Part 2: Inter-patient, intra-subtype dispersion ------------------------
  df_IPH <- t(vsd_sub) %>% as.data.frame(check.names = FALSE)
  df_IPH$Subtype <- anno$Subtype
  
  gene_cols <- setdiff(names(df_IPH), "Subtype")
  disp_long <- df_IPH %>%
    pivot_longer(all_of(gene_cols), names_to = "gene", values_to = "value") %>%
    group_by(Subtype, gene) %>%
    summarise(
      n          = sum(!is.na(value)),
      mean       = mean(value, na.rm = TRUE),
      variance   = if (n() > 1) var(value, na.rm = TRUE) else NA_real_,
      dispersion = dispersion_safe(value),
      .groups = "drop"
    )
  
  disp_wide <- disp_long %>%
    dplyr::select(Subtype, gene, dispersion) %>%
    tidyr::pivot_wider(names_from = gene, values_from = dispersion)
  
  # tidy to matrix-like with genes as rows, Subtype as columns
  disp_wide <- disp_wide %>% t() %>% as.data.frame(check.names = FALSE)
  colnames(disp_wide) <- disp_wide[1, ]
  disp_wide <- disp_wide[-1, , drop = FALSE]
  
  # Keep only requested genes & keep their order
  disp_wide <- disp_wide[found_genes, , drop = FALSE]
  
  # ---- Part 3: Heatmap assembly & saving -------------------------------------
  # Column order (keep only columns that exist)
  col_order <- c("Small Cell","NEPC","CRPC","ICC/IDC","Adeno","Benign","Healthy")
  col_order <- intersect(col_order, colnames(disp_wide))
  mat <- disp_wide[, col_order, drop = FALSE]
  
  # Derive gene groups from heatmap_row_split on the (found) genes
  n_adc    <- heatmap_row_split[1]
  n_immune <- heatmap_row_split[2]
  n_ar     <- heatmap_row_split[3]
  
  # Reconcile splits with possibly missing genes
  # Choose indices within found_genes (current row order = found_genes)
  idx_all <- seq_len(nrow(mat))
  idx_adc    <- idx_all[seq_len(min(n_adc, nrow(mat)))]
  idx_immune <- setdiff(idx_all, idx_adc)
  idx_immune <- idx_immune[seq_len(min(n_immune, length(idx_immune)))]
  idx_ar     <- setdiff(idx_all, c(idx_adc, idx_immune))
  
  row_order <- c(idx_adc, idx_immune, idx_ar)
  mat <- mat[row_order, , drop = FALSE]
  
  gene_group <- rep(NA_character_, nrow(mat))
  gene_group[seq_along(idx_adc)] <- "ADC"
  gene_group[seq(length(idx_adc) + 1, length(idx_adc) + length(idx_immune))] <- "Immune"
  gene_group[(length(idx_adc) + length(idx_immune) + 1):nrow(mat)] <- "AR reg."
  gene_group <- factor(gene_group, levels = c("ADC","Immune","AR reg."))
  
  # Replace rownames with pretty gene_names (aligned to current rows)
  # Map found_genes -> gene_names by original position in 'genes'
  gene_to_pretty <- setNames(gene_names, genes)
  pretty_names   <- gene_to_pretty[rownames(mat)]
  # Fallback to original symbol if pretty missing
  pretty_names[is.na(pretty_names)] <- rownames(mat)[is.na(pretty_names)]
  rownames(mat) <- pretty_names
  
  # Color strip for groups
  grp_cols <- c("ADC"="#FF0000","Immune"="#49FF00","AR reg."="#0092FF")
  row_anno <- rowAnnotation(
    Group = gene_group,
    col   = list(Group = grp_cols),
    width = unit(3, "mm"),
    show_annotation_name = FALSE,
    annotation_legend_param = list(title = "Gene")
  )
  
  # Numeric coercion + log1p
  rn <- rownames(mat); cn <- colnames(mat)
  for (j in seq_len(ncol(mat))) mat[, j] <- as.numeric(mat[, j])
  dimnames(mat) <- list(rn, cn)
  mat <- log1p(as.matrix(mat))
  
  # Heatmap palette (auto-scaled)
  rng <- range(as.numeric(mat), na.rm = TRUE)
  col_fun <- colorRamp2(c(rng[1], mean(rng), rng[2]), c("#1f2a44", "white", "#b40426"))
  
  # In-cell labels
  cell_labeller <- function(j, i, x, y, w, h, fill) {
    v <- mat[i, j]
    lab <- if (is.na(v)) "" else sprintf("%.2f", v)
    grid::grid.text(lab, x, y, gp = grid::gpar(fontsize = 7))
  }
  
  ht <- Heatmap(
    mat, name = "value",
    col = col_fun,
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_names_side = "left",
    row_gap = unit(x = c(1, 1), units = "mm"),
    column_names_side = "top",
    left_annotation = row_anno,
    row_split = gene_group,
    cell_fun = cell_labeller,
    row_title = NULL
  )
  
  # ---- Save to files ----------------------------------------------------------
  fig_dir <- file.path(project_dir, "figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  
  ts   <- format(Sys.time(), "%Y%m%d_%H%M%S")
  base <- file.path(fig_dir, paste0("dispersion_heatmap_", ts))
  
  # TIFF
  tiff(filename = paste0(base, ".tiff"), width = 2400, height = 1600, res = 300, compression = "lzw")
  ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  
  # PNG
  png(filename = paste0(base, ".png"), width = 2400, height = 1600, res = 300)
  ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  
  # PDF
  pdf(file = paste0(base, ".pdf"), width = 10, height = 7)
  ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  
  message("Saved: ",
          paste0(base, c(".tiff", ".png", ".pdf"), collapse = " , "))
  
  invisible(list(
    heatmap = ht,
    matrix  = mat,
    annotation = anno
  ))
}


# ---- Example (commented) ------------------------------------------------------
result <- make_dispersion_heatmap(
  obj = obj_sub,
  project_dir = "/projects/eshenderov-hpc/Shivang/projects/B7-H3 ADC Manuscript",
  genes = c("STEAP1", "STEAP2", "PSCA", "FOLH1", "KLK2", "TACSTD2", "NECTIN1", "NECTIN4", "DLL3",
            "CD276", "CD274", "PDCD1LG2", "PDCD1", "LAG3", "TNFRSF4", "TNFRSF9", "CTLA4",
            "AR", "KLK3", "TMPRSS2", "FKBP5", "IGF1R"),
  gene_names = c("STEAP1", "STEAP2", "PSCA", "PSMA", "KLK2", "TROP-2", "NECTIN1", "NECTIN4", "DLL3",
                 "B7-H3", "PD-1", "PD-L2", "PD-L1", "LAG3", "OX40", "4-1BB", "CTLA4", 
                 "AR", "KLK3", "TMPRSS2", "FKBP5", "IGF1R"),
  heatmap_row_split = c(9, 8, 5)   # ADC, Immune, AR reg.
)

