library(readxl)
library(tidyverse)
library(limma)
library(RColorBrewer)
library(ggsci)
library(ggpubr)
library(ggplotify)
library(data.table)
library(matrixStats)
library(circlize)
library(patchwork)
library(ComplexHeatmap)
library(gplots)

project_dir <- "/projects/eshenderov-hpc/Shivang/projects/B7-H3 ADC Manuscript"

combine_duplicate_genes <- function(counts_data) {
  # Provide a data table with first column containing gene names and the name of
  # the first column should be "id". The data table is supplied as an argument
  # to counts_data variable. The output will be a data table with first column
  # containing gene names and named as "gene". The rows with same gene names
  # will be summed together.
  simplified_counts_data <- counts_data %>%
    group_by(c(id)) %>%
    summarize(across(.cols = -id, .fns = sum))
  colnames(simplified_counts_data)[1] <- "gene"
  simplified_counts_data <-
    simplified_counts_data[!(simplified_counts_data$gene %in% c("")), ]
  simplified_counts_data <-
    simplified_counts_data[!is.na(simplified_counts_data$gene), ]
  return(simplified_counts_data %>% as.data.table())
  
}

fpkm_to_logTPM <- function(log1p_fpkm_df) {
  
  fpkm_df <- 2^log1p_fpkm_df
  
  # Convert to matrix to ensure numeric operations
  fpkm <- as.matrix(fpkm_df)
  
  # Step 1. FPKM → TPM
  tpm <- apply(fpkm, 2, function(x) {
    (x / sum(x, na.rm = TRUE)) * 1e6
  })
  
  # Step 2. Library size normalization (optional — TPM already normalized)
  lib_size <- colSums(tpm, na.rm = TRUE)
  tpm_norm <- t( t(tpm) / lib_size * mean(lib_size) )
  
  # Step 3. log1p transform
  tpm_log <- log1p(tpm_norm)
  
  # Return as data frame with same row/col names
  tpm_log_df <- as.data.frame(tpm_log)
  rownames(tpm_log_df) <- rownames(log1p_fpkm_df)
  colnames(tpm_log_df) <- colnames(log1p_fpkm_df)
  
  return(tpm_log_df)
}

RLE_plot <- function(data, show_legend = FALSE, show_outliers = FALSE) {
  # Provie a matrix of genesx samples, with gene names as rownames and
  # sample names as column names
  library(matrixStats)
  if (is.matrix(data)) {
    
    RLE_mat <- data
    
  } else {
    
    RLE_mat <- as.matrix(data)
    
  }
  RLE_mat <- sweep(x = RLE_mat, MARGIN = 1, STATS = matrixStats::rowMedians(RLE_mat), FUN = "-")
  RLE_df <- as.data.frame(RLE_mat)
  RLE_df <- melt(RLE_df, variable.name = "Patient_ID", value.name = "RLE")
  
  whiskers <- RLE_df %>%
    group_by(Patient_ID) %>%
    summarize(lower = quantile(RLE, 0.25) - 1.5*IQR(RLE),
              upper = quantile(RLE, 0.75) + 1.5*IQR(RLE),
              median = median(RLE),
              RLE = mean(RLE))
  whiskers$discrete_x <- c(1:dim(whiskers)[1])
  ylim_abs <- max(abs(whiskers$lower), whiskers$upper)
  
  plot_output <- ggplot(RLE_df, aes(x = Patient_ID, y = RLE, fill = Patient_ID)) +
    scale_x_discrete() +
    geom_segment(data = whiskers, aes(x = discrete_x - 0.2, xend = discrete_x + 0.2, y = upper, yend = upper)) +
    geom_segment(data = whiskers, aes(x = discrete_x - 0.2, xend = discrete_x + 0.2, y = lower, yend = lower)) +
    geom_linerange(data = whiskers, aes(x = Patient_ID, ymin = lower, ymax = upper),
                   size = 0.2, linetype = "dashed") +
    geom_boxplot(width = 0.4, color = "black", show.legend = show_legend,
                 outlier.shape = if (show_outliers) "circle" else NA, coef = 0) +
    scale_y_continuous(limits = c(-ylim_abs, ylim_abs)) +
    theme(axis.line = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", color = "black", linewidth = 2),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank()
    ) +
    labs(x = "Samples", y = "Relative expression")
  print(plot_output)
  
}

expression_map <- function (expr, meta, genes, gene_names, heatmap_row_split, project_dir, file_name) {
  
  # Calculate AR and NE Score ----------------------------------------------------
  vsd <- expr
  features <- read_xlsx(path = file.path("/projects/eshenderov-hpc/Shivang/projects/scRNA-seq_atlas/Prostate/integrated_dataset", "Genesets_Emir.xlsx")) %>%
    as.list()
  features$`AR Score` <- features$`AR Score`[!(features$`AR Score` %in% NA)]
  features$`NE Score` <- features$`NE Score`[!(features$`NE Score` %in% NA)]
  AR_Score <- (vsd[features$`AR Score`, ] %>% as.matrix() %>% colSums2())/(features$`AR Score` %>% length())
  NE_Score <- (vsd[features$`NE Score`, ] %>% as.matrix() %>% colSums2())/(features$`NE Score` %>% length())
  
  # Subset expression data matrix for genes of interest --------------------------
  luminal_basal_genes <- genes
  vsd_sub <- vsd[luminal_basal_genes, ]
  
  anno <- meta
  
  # Assign colors to subtypes of PCa, AR and NE Score in heatmap -----------------
  #colors <- brewer.pal(9, "Set1")
  colors <- pal_npg("nrc")(meta$Subtype %>% levels() %>% length())
  names(colors) <- meta$Subtype %>% levels()
  anno_colors_col <- list(Subtype = rev(colors),
                          "AR_Score" = colorRamp2(breaks = seq(min(AR_Score), max(AR_Score), length.out = 100),
                                                  hcl_palette = "SunsetDark",
                                                  transparency = 0,
                                                  reverse = TRUE),
                          "NE_Score" = colorRamp2(breaks = seq(min(NE_Score), max(NE_Score), length.out = 100),
                                                  hcl_palette = "Purples 3",
                                                  transparency = 0,
                                                  reverse = TRUE)
  )
  
  # Assign colors to rows of genes based on their category -----------------------
  gene_group <- factor(rep(c("ADC", "AR reg."), times = heatmap_row_split),
                       levels = c("ADC", "AR reg."))
  
  # colors for the row group bar
  anno_colors_row <- list(Gene = c("ADC" = "#FF0000", "AR reg." = "#0092FF"))
  
  # --- row-side annotation (use rowAnnotation, not HeatmapAnnotation) ---
  row_ha <- rowAnnotation(
    Gene = gene_group,
    col  = anno_colors_row,
    show_annotation_name = FALSE,     # <-- hide "Gene"
    width = unit(3, "mm"),
    annotation_legend_param = list(
      Gene = list(
        title      = "Gene",                 # (or set to "" for no title)
        title_gp   = gpar(fontsize = 30, fontface = "bold"),
        labels_gp  = gpar(fontsize = 30),
        grid_width = unit(8, "mm"),                  # bigger legend keys
        grid_height= unit(12, "mm"),
        title_gap  = unit(3, "mm"),                  # <-- more space title ↔ body
        row_gap    = unit(3, "mm")                   # space between legend rows
      )
    )
  )
  
  # Create heatmap column section ------------------------------------------------
  
  at_ar     <- pretty(range(AR_Score, na.rm = TRUE), n = 5)
  lab_ar    <- as.character(at_ar)
  lab_ar[length(lab_ar)] <- ""
  lab_ar[1] <- ""
  
  at_ne     <- pretty(range(NE_Score, na.rm = TRUE), n = 5)
  lab_ne    <- as.character(at_ne)
  lab_ne[length(lab_ne)] <- ""
  lab_ne[1] <- ""
  
  col_ha <- HeatmapAnnotation(
    Subtype   = anno$Subtype,
    "AR_Score" = AR_Score %>% as.numeric() %>% unname(),
    "NE_Score" = NE_Score %>% as.numeric() %>% unname(),
    col = anno_colors_col,
    annotation_label = c("Subtype", "AR Score", "NE Score"),
    annotation_name_side = "right",
    annotation_name_gp   = gpar(fontsize = 30, fontface = "bold"),
    border = c(Subtype = FALSE, AR_Score = FALSE, NE_Score = FALSE),
    simple_anno_size = unit(16, "mm"),
    
    ## >>> EDIT HERE: legend title/text sizes for the TOP annotations <<<
    annotation_legend_param = list(
      Subtype = list(
        title_gp  = gpar(fontsize = 30, fontface = "bold"),   # legend title
        labels_gp = gpar(fontsize = 30),                       # legend labels
        grid_width  = unit(8, "mm"),   # smaller color keys
        grid_height = unit(12, "mm"),
        title_gap = unit(3, "mm"), 
        row_gap    = unit(3, "mm")
        #ncol = 2
        # grid_width = unit(4, "mm"), grid_height = unit(4, "mm")  # (optional) key size
      ),
      AR_Score = list(
        title = "AR Score",                                   # pretty title (optional)
        title_gp  = gpar(fontsize = 30, fontface = "bold"),
        labels_gp = gpar(fontsize = 30),
        legend_height = unit(60, "mm"),                        # (optional) gradient height
        legend_width  = unit(8, "mm"),   # thinner bar (reduces right-side width)
        at = at_ar, 
        labels = lab_ar,
        title_gap = unit(100, "mm")  # fewer ticks/labels
      ),
      NE_Score = list(
        title = "NE Score",
        title_gp  = gpar(fontsize = 30, fontface = "bold"),
        labels_gp = gpar(fontsize = 30),
        legend_height = unit(60, "mm"),
        legend_width  = unit(8, "mm"),
        at = at_ne,
        labels = lab_ne,
        title_gap = unit(100, "mm")
      )
    )
  )
  
  min_v <- min(vsd_sub, na.rm = TRUE)
  max_v <- max(vsd_sub, na.rm = TRUE)
  
  #at_expr <- pretty(range(vsd_sub, na.rm = TRUE), n = 3)
  #at_expr <- head(at_expr, -1)   # <-- remove the max tick (prevents crowding at top)
  
  at_expr   <- pretty(range(vsd_sub,  na.rm = TRUE), n = 5)
  lab_expr  <- as.character(at_expr)
  lab_expr[length(lab_expr)] <- ""  # blank top
  lab_expr[1] <- ""

  n_rows  <- nrow(vsd_sub)
  n_cols  <- ncol(vsd_sub)
  
  cell_h_mm <- 18        # mm per row (you already set height via this)
  cell_w_mm <- 8        # mm per column
  marg_mm   <- 20    # extra margin mm
  
  ht_plot <- Heatmap(
    vsd_sub %>% as.matrix(),
    name = "Expression",
    height = unit(n_rows * cell_h_mm, "mm"),   # per-row height
    width  = unit(n_cols * cell_w_mm, "mm"),   # ← per-column width
    col  = colorpanel(500, low = "#440154FF", mid = "#21908CFF", high = "#FDE725FF"),
    cluster_rows = FALSE, cluster_columns = FALSE,
    top_annotation  = col_ha,
    left_annotation = row_ha,
    row_labels = gene_names,
    row_names_side = "left",
    row_names_gp   = gpar(fontsize = 30),
    row_split = gene_group,                 # keep the grouping
    row_gap   = unit(1, "mm"),              # visible separation between groups
    row_title = NULL,                       # <-- no row category labels
    rect_gp   = gpar(col = "grey90", lwd = 0.5),
    show_column_names = FALSE,
    column_names_gp   = gpar(fontsize = 0),
    column_title      = NULL,
    heatmap_legend_param = list(
      title_gp  = gpar(fontsize = 30, fontface = "bold"),
      labels_gp = gpar(fontsize = 30),
      legend_width  = unit(8, "mm"),
      legend_height = unit(60, "mm"),
      at = at_expr,
      labels = lab_expr,
      title_gap = unit(100, "mm")
    )
  )
  
  # Save heatmap into pdf --------------------------------------------------------
  
  # png(file = file.path(project_dir, file_name), 
  #     width = convertUnit(x = (ncol(vsd_sub) * unit(1200*0.7, "mm")) * 1.7, 
  #                         unitTo = "inches"), 
  #     height = convertUnit(x = (nrow(vsd_sub) + 3) * unit(450*0.7, "mm") * 1.7, 
  #                          unitTo = "inches"))
  # ComplexHeatmap::draw(ht_plot)
  # dev.off()
  
  height_in <- (n_rows * cell_h_mm + marg_mm) / 20
  width_in  <- (n_cols * cell_w_mm + marg_mm) / 12
  dpi <- 600
  max_px <- 32000
  width_in_bitmap  <- min(width_in,  max_px / dpi)
  height_in_bitmap <- min(height_in, max_px / dpi)
  
  tiff(paste0(file.path(project_dir, file_name), ".tiff"),
       width = width_in_bitmap, height = height_in_bitmap,
       units = "in", res = dpi, compression = "lzw")
  ComplexHeatmap::draw(ht_plot)
  dev.off()
  
  pdf(paste0(file.path(project_dir, file_name), ".pdf"),
      width = width_in, height = height_in)
  ComplexHeatmap::draw(ht_plot)
  dev.off()
  
  # pdf(file = file.path(project_dir, file_name), 
  #     width = convertUnit(x = (ncol(vsd_sub) * unit(2, "mm")) * 1.7, 
  #                         unitTo = "inches"), 
  #     height = convertUnit(x = (nrow(vsd_sub) + 3) * unit(3.5, "mm") * 1.7, 
  #                          unitTo = "inches"))
  # ComplexHeatmap::draw(ht_plot,
  #                      heatmap_legend_side    = "bottom",
  #                      annotation_legend_side = "bottom")
  # dev.off()
  
  # width_in <- convertUnit(x = (ncol(vsd_sub) * unit(3, "mm")) * 1.7, unitTo = "inches", valueOnly = TRUE)
  # height_in <- convertUnit(x = (nrow(vsd_sub) + 3) * unit(9, "mm") * 1.7, unitTo = "inches", valueOnly = TRUE)
  # 
  # # Set resolution (dpi)
  # dpi <- 600
  # 
  # # Convert to pixels
  # width_px <- width_in * dpi
  # height_px <- height_in * 1.7 * dpi
  # 
  # # Save heatmap into TIFF -------------------------------------------------------
  # tiff(filename = file.path(project_dir, file_name),
  #      width = width_px, height = height_px, res = dpi, units = "px", compression = "lzw")
  # 
  # ComplexHeatmap::draw(ht_plot)
  # 
  # dev.off()
  
  invisible(ht_plot)
  
}


corr_heatmap_genes <- function(expr_df, genes, project_dir,
                               method = c("pearson","spearman"),
                               file = NULL, width = 7, height = 7,
                               show_gene_labels = TRUE) {
  method <- match.arg(method)
  
  # --- Coerce to numeric matrix (genes as rows, samples as columns)
  stopifnot(!is.null(rownames(expr_df)))
  mat <- as.matrix(expr_df)
  storage.mode(mat) <- "numeric"
  
  # --- Subset to requested genes (keep user order)
  genes <- unique(genes)
  genes_in <- genes[genes %in% rownames(mat)]
  missing <- setdiff(genes, genes_in)
  if (length(genes_in) < 2L) stop("Fewer than 2 requested genes were found in the matrix.")
  
  M <- mat[genes_in, , drop = FALSE]
  
  # --- Drop zero-variance genes (they yield NA correlations)
  keep <- apply(M, 1, function(x) sd(x, na.rm = TRUE) > 0)
  dropped_const <- rownames(M)[!keep]
  M <- M[keep, , drop = FALSE]
  if (nrow(M) < 2L) stop("After removing constant genes, fewer than 2 genes remain.")
  
  # --- Gene–gene correlation across samples
  cor_mat <- suppressWarnings(cor(t(M), use = "pairwise.complete.obs", method = method))
  
  # --- Heatmap (diverging palette, symmetric around 0)
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE) ||
      !requireNamespace("circlize", quietly = TRUE)) {
    stop("Please install packages: ComplexHeatmap and circlize")
  }
  library(ComplexHeatmap)
  library(circlize)
  
  col_fun <- circlize::colorRamp2(c(-1, 0, 1),
                                  c("#3B4CC0", "#F7F7F7", "#B40426"))  # crisp RdBu-like
  
  ht_plot <- Heatmap(cor_mat,
                name = paste0(method, " r"),
                col = col_fun,
                cluster_rows = TRUE, cluster_columns = TRUE,
                show_row_names = show_gene_labels,
                show_column_names = show_gene_labels,
                row_names_gp = grid::gpar(fontsize = 8),
                column_names_gp = grid::gpar(fontsize = 8),
                heatmap_legend_param = list(at = c(-1, -0.5, 0, 0.5, 1)))
  
  # Save heatmap into pdf --------------------------------------------------------
  square_dim <- convertUnit(x = (50 * unit(200*0.7, "mm")) * 1.7, 
                            unitTo = "inches")
  
  png(file = file.path(project_dir, "figures", 'FigPDX_c.png'), 
      width = square_dim, 
      height = square_dim)
  ComplexHeatmap::draw(ht_plot)
  dev.off()
  
  square_dim <- convertUnit(x = (length(genes) * unit(4, "mm")) * 1.7, unitTo = "inches", valueOnly = TRUE)
  
  pdf(file = file.path(project_dir, "figures", 'FigPDX_c.pdf'), 
      width = square_dim, 
      height = square_dim)
  ComplexHeatmap::draw(ht_plot)
  dev.off()
  
  square_dim <- convertUnit(x = (length(genes) * unit(4, "mm")) * 1.7, unitTo = "inches", valueOnly = TRUE)
  
  # Set resolution (dpi)
  dpi <- 1200
  
  # Convert to pixels
  width_px <- square_dim * dpi
  height_px <- square_dim * dpi
  
  # Save heatmap into TIFF -------------------------------------------------------
  tiff(filename = file.path(project_dir, "figures", 'FigPDX_c.tif'),
       width = width_px, height = height_px, res = dpi, units = "px", compression = "lzw")
  
  ComplexHeatmap::draw(ht_plot)
  
  dev.off()
  
  invisible(list(cor = cor_mat,
                 heatmap = ht_plot,
                 missing_genes = missing,
                 dropped_constant_genes = dropped_const))
}


expr_line_strips <- function(expr_df, genes = NULL, project_dir,
                             order = c("hclust","heatmap","none"),
                             heatmap_ht = NULL,
                             scale_method = c("zscore","minmax","none"),
                             cluster_dist = c("correlation","euclidean"),
                             cluster_method = "complete",
                             zlim = c(-3, 3),
                             point_size = 12, line_size = 6,
                             show_x_labels = TRUE, 
                             hide_y_labels_at = c(-2, 0, 2), 
                             plot_name = "Fig_PDXb", 
                             gene_names) {
  # --- checks & args
  stopifnot(!is.null(rownames(expr_df)), !is.null(colnames(expr_df)))
  order <- match.arg(order)
  scale_method <- match.arg(scale_method)
  cluster_dist <- match.arg(cluster_dist)
  
  # --- coerce numeric matrix, keep sample order
  M <- as.matrix(expr_df); storage.mode(M) <- "numeric"
  sample_order <- colnames(M)
  
  # --- subset genes (keep user order for "none")
  if (is.null(genes)) genes <- rownames(M)
  genes <- unique(genes)
  genes_in <- genes[genes %in% rownames(M)]
  gene_names_in <- gene_names[genes %in% rownames(M)]
  missing <- setdiff(genes, genes_in)
  if (length(genes_in) < 1L) stop("None of the requested genes are present.")
  M <- M[genes_in, , drop = FALSE]
  
  # --- per-gene scaling for comparable y ranges
  if (scale_method == "zscore") {
    Msc <- t(scale(t(M)))      # z-score by gene
    Msc[is.na(Msc)] <- 0       # constant genes -> 0
    Y <- Msc; ylab <- "Z-Score"
  } else if (scale_method == "minmax") {
    rng <- apply(M, 1, function(x) diff(range(x, na.rm = TRUE)))
    minv <- apply(M, 1, min, na.rm = TRUE)
    Y <- (M - minv) / ifelse(rng == 0, 1, rng)
    ylab <- "Scaled expression (min–max)"
  } else {
    Y <- M; ylab <- "Expression"
  }
  
  rownames(Y) <- gene_names_in
  
  # --- gene order (optionally from ComplexHeatmap object)
  if (order == "heatmap" && !is.null(heatmap_ht)) {
    if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
      stop("ComplexHeatmap needed for 'heatmap' order.")
    }
    ord <- ComplexHeatmap::row_order(heatmap_ht)
    ord <- unlist(ord, use.names = FALSE)
    gene_order <- rownames(Y)[ord[ord <= nrow(Y)]]
    gene_order <- gene_order[gene_order %in% rownames(Y)]
  } else if (order == "hclust") {
    if (cluster_dist == "correlation") {
      D <- as.dist(1 - stats::cor(t(Y), use = "pairwise.complete.obs"))
    } else {
      D <- dist(Y)
    }
    hc <- hclust(D, method = cluster_method)
    gene_order <- rownames(Y)[hc$order]
  } else {
    gene_order <- rownames(Y)
  }
  
  # --- long table for ggplot
  df_long <- as.data.frame(Y)
  df_long$Gene <- rownames(Y)
  df_long <- tidyr::pivot_longer(df_long, -Gene,
                                 names_to = "Sample", values_to = "y")
  df_long$Gene   <- factor(df_long$Gene,   levels = gene_order)
  df_long$Sample <- factor(df_long$Sample, levels = sample_order)
  
  # --- plot
  suppressPackageStartupMessages(library(ggplot2))
  p <- ggplot(df_long, aes(x = Sample, y = y)) +
    geom_line(aes(group = Gene), linewidth = line_size) +
    geom_point(size = point_size) +
    facet_grid(Gene ~ ., scales = "fixed", switch = "y") +
    labs(x = "Samples", y = ylab) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y      = element_blank(),
      axis.ticks.y     = element_blank(),
      strip.placement = "outside",
      strip.background.y = element_rect(fill = "grey95", colour = NA), 
      axis.title.x = element_text(size = 40, face = "bold", margin = margin(t = 6)),
      axis.title.y = element_text(size = 40, face = "bold", margin = margin(r = 6)),
      strip.text.y.left = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 40),
      axis.text.x = if (show_x_labels) element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14) else element_blank(),
      axis.ticks.x = if (show_x_labels) element_line() else element_blank(),
      panel.spacing.y = unit(0.15, "lines")
    ) +
    coord_cartesian(ylim = if (!is.null(zlim) && scale_method == "zscore") zlim else NULL)
  
  # Save heatmap into pdf --------------------------------------------------------
  p
  ggsave(file.path(project_dir, "figures", paste0(plot_name, '.tiff')), device = "tiff", width = 540, height = 120, units = "mm", dpi = 600)
  
  invisible(list(plot = p,
       gene_order = gene_order,
       missing_genes = missing,
       scaled_matrix = Y))
}

# expr: genes x samples (rows = genes, cols = samples)
# genes: character vector of genes to keep (in this order)
gene_row_derivative <- function(expr, genes = NULL) {
  stopifnot(!is.null(rownames(expr)), !is.null(colnames(expr)))
  
  # ensure numeric matrix
  M <- as.matrix(expr)
  storage.mode(M) <- "numeric"
  
  # subset to requested genes, keep user order
  if (is.null(genes)) {
    genes_in <- rownames(M)
  } else {
    genes <- unique(genes)
    genes_in <- genes[genes %in% rownames(M)]
    missing <- setdiff(genes, genes_in)
    if (length(missing)) {
      warning("These genes were not found and will be skipped: ",
              paste(missing, collapse = ", "))
    }
  }
  M <- M[genes_in, , drop = FALSE]
  
  # discrete derivative along samples with baseline sample(0)=0:
  # D[,1] = M[,1] - 0;  D[,i] = M[,i] - M[,i-1] for i>=2
  ns <- ncol(M)
  D <- matrix(NA_real_, nrow = nrow(M), ncol = ns,
              dimnames = dimnames(M))
  if (ns >= 1) D[, 1] <- M[, 1]
  if (ns >= 2) D[, 2:ns] <- M[, 2:ns, drop = FALSE] - M[, 1:(ns - 1), drop = FALSE]
  
  as.data.frame(D)
}

plot_intact_vs_castrated <- function(
    df, project_dir, sample_col = "Sample_Names", y_col = "B7-H3",
    test = c("t.test","wilcox.test"),
    x_title_size = 20, y_title_size = 20,
    x_text_size  = 20, y_text_size  = 20,
    file_name = "figPDX_f.tiff"
) {
  test <- match.arg(test)
  
  df1 <- df %>%
    dplyr::mutate(
      sample = .data[[sample_col]],
      y      = as.numeric(.data[[y_col]]),
      is_cr  = grepl("CR\\d*$", sample),                                   # CR, CR1, CR2, ...
      status = factor(ifelse(is_cr, "Castrated", "Intact"),
                      levels = c("Intact","Castrated")),                   # Intact first
      base   = sub("CR\\d*$", "", sample)                                  # strip trailing CR*
    ) %>%
    dplyr::filter(!is.na(y))
  
  # keep only bases that have both states
  bases_with_both <- df1 %>%
    dplyr::distinct(base, status) %>% dplyr::count(base) %>%
    dplyr::filter(n == 2) %>% dplyr::pull(base)
  
  df2 <- df1 %>% dplyr::filter(base %in% bases_with_both)
  
  # per-patient means for paired testing / connecting line
  patient_means <- df2 %>%
    dplyr::group_by(base, status) %>%
    dplyr::summarise(y = median(y, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = status, values_from = y) %>%
    tidyr::drop_na(Intact, Castrated)
  
  # paired p-value (means per patient)
  pval <- if (nrow(patient_means) > 1) {
    if (test == "t.test") {
      t.test(patient_means$Intact, patient_means$Castrated, paired = TRUE)$p.value
    } else {
      wilcox.test(patient_means$Intact, patient_means$Castrated, paired = TRUE, exact = FALSE)$p.value
    }
  } else NA_real_
  
  # compute per-group means for robust red lines
  means_df <- df2 %>%
    dplyr::group_by(status) %>%
    dplyr::summarise(y_mean = mean(y, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(xn = as.numeric(status))
  
  # y range + bracket position kept inside limits
  y_limits <- c(0, 9)
  #y_top <- min(y_limits[2] - 0.2, max(df2$y, na.rm = TRUE) + 0.2)
  y_top <- 8
  lbl   <- if (is.na(pval)) NULL else paste0("p = ", format.pval(pval, digits = 3, eps = 1e-3))
  
  p <- ggplot(df2, aes(x = status, y = y)) +
    # one line per patient (means)
    geom_segment(data = patient_means,
                 aes(x = "Intact", xend = "Castrated",
                     y = Intact,   yend = Castrated,
                     group = base),
                 inherit.aes = FALSE, linewidth = 0.6, alpha = 0.8) +
    # all replicate dots
    geom_point(size = 2, position = position_jitter(width = 0.06, height = 0)) +
    # robust horizontal red mean lines per group
    geom_segment(
      data = means_df,
      aes(x = xn - 0.15, xend = xn + 0.15, y = y_mean, yend = y_mean),
      inherit.aes = FALSE, linewidth = 2, color = "red", lineend = "round"
    ) +
    labs(
      title = y_col,
      x = NULL,
      y = expression(log[2] * " TPM")
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, size = x_title_size, face = "bold"), 
      axis.title.x = element_text(size = x_title_size, face = "bold", margin = margin(t = 6), color = "black"),
      axis.title.y = element_text(size = y_title_size, face = "bold", margin = margin(r = 6), color = "black"),
      axis.text.x  = element_text(size = x_text_size,  color = "black"),
      axis.text.y  = element_text(size = y_text_size,  color = "black")
    ) +
    coord_cartesian(ylim = y_limits) +
    scale_y_continuous(breaks = seq(y_limits[1], y_limits[2], by = 1))
  
  if (!is.null(lbl)) {
    p <- p + ggpubr::stat_pvalue_manual(
      data.frame(group1 = "Intact", group2 = "Castrated",
                 y.position = y_top, label = lbl),
      label = "label",
      xmin = "group1", xmax = "group2", y.position = "y.position",
      size = 8,       # optional: bold label
      family = "sans",         # optional: font family
      color = "black",         # color for label & bracket
      bracket.size = 1.4,      # ↑ bracket line thickness
      tip.length = 0.02        # optional: slightly longer tips
    )
  }
  
  dir.create(file.path(project_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(project_dir, "figures", paste0(file_name, ".tiff")),
         plot = p, device = "tiff", width = 120, height = 120, units = "mm", dpi = 600)
  ggsave(file.path(project_dir, "figures", paste0(file_name, ".pdf")),
         plot = p, device = "pdf", width = 120, height = 120, units = "mm", dpi = 600)
  
  invisible(list(
    plot = p,
    patient_means = patient_means,
    group_means = means_df,
    paired_p = pval,
    used_bases = bases_with_both
  ))
}

plot_ar_status <- function(
    df, project_dir,
    status_col = "AR Status",     # column with "AR+" / "AR-"
    y_col      = "B7-H3",         # expression column (log2 TPM)
    method     = c("t.test","wilcox.test"),
    title      = "AR+ vs AR-",
    width_mm   = 120, height_mm = 120, dpi = 600,
    file_base  = "FigPDX_d",
    x_title_size = 20, y_title_size = 20,
    x_text_size  = 20, y_text_size  = 20
) {
  method <- match.arg(method)
  
  dat_ar <- df %>%
    dplyr::mutate(
      B7_H3     = as.numeric(.data[[y_col]]),
      AR_Status = factor(.data[[status_col]], levels = c("AR+","AR-"))
    ) %>%
    dplyr::filter(!is.na(B7_H3), !is.na(AR_Status))
  
  # y position for p-value bracket
  y_top <- max(dat_ar$B7_H3, na.rm = TRUE) + 0.6
  
  # per-group means for red lines
  means_df <- dat_ar %>%
    dplyr::group_by(AR_Status) %>%
    dplyr::summarise(y_mean = median(B7_H3, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(xn = as.numeric(AR_Status))
  
  y_top <- 8
  
  p <- ggplot(dat_ar, aes(x = AR_Status, y = B7_H3)) +
    geom_point(position = position_jitter(width = 0.08, height = 0), size = 2) +
    geom_segment(
      data = means_df,
      aes(x = xn - 0.15, xend = xn + 0.15, y = y_mean, yend = y_mean),
      inherit.aes = FALSE, linewidth = 2, color = "red", lineend = "round"
    ) +
    ggpubr::stat_compare_means(
      comparisons = list(c("AR+","AR-")),
      method = method, label = "p.format",
      size = 8,                 # ↑ p-value text size
      color = "black",          # label + bracket color
      bracket.size = 1.4,       # ↑ bracket thickness
      tip.length = 0.02,        # a touch longer tips
      vjust = -0.2,             # nudge up from the bracket
      label.y = y_top
    ) +
    labs(title = y_col, x = NULL, y = expression("log"[2]*" TPM")) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, size = 16,face = "bold"),
      axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 6), color = "black"),
      axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 6), color = "black"),
      axis.text.x  = element_text(size = 20, color = "black"),
      axis.text.y  = element_text(size = 20, color = "black")
    ) +
    coord_cartesian(ylim = c(0, 9)) +
    scale_y_continuous(breaks = seq(0, 9, by = 1))
  
  # save
  outdir <- file.path(project_dir, "figures")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(outdir, paste0(file_base, ".tiff")),
         plot = p, device = "tiff",
         width = width_mm, height = height_mm, units = "mm", dpi = dpi)
  ggsave(file.path(outdir, paste0(file_base, ".pdf")),
         plot = p, device = "pdf",
         width = width_mm, height = height_mm, units = "mm")
  
  invisible(list(plot = p, data = dat_ar, means = means_df))
}

make_coexp_pie <- function(
    expr, gene_x, gene_y,
    threshold = 1,
    out_dir = file.path(getwd(), "figures"),
    out_basename = NULL,
    legend_title = NULL,
    gene_x_color = "#4DBBD5FF",
    gene_y_color = "#00A087FF",
    both_color   = "#E64B35FF",
    none_color   = "#999999FF",
    # --- sizing ---
    lock_panel = TRUE,
    panel_mm = 90,            # fixed pie diameter
    legend_side = "right",    # "right" or "bottom"
    legend_space_mm = 80,     # reserved width (right) or height (bottom) for legend
    margins_mm = 10,          # outer plot margins on all sides
    width_mm = NULL,          # if NULL, computed from panel_mm + legend_space_mm + margins
    height_mm = NULL,
    dpi = 600
) {
  stopifnot(is.matrix(expr) || is.data.frame(expr))
  expr <- as.data.frame(expr, stringsAsFactors = FALSE)
  if (is.null(rownames(expr))) stop("`expr` must have rownames = gene symbols.")
  if (!all(c(gene_x, gene_y) %in% rownames(expr))) {
    stop("Missing genes: ", paste(setdiff(c(gene_x, gene_y), rownames(expr)), collapse = ", "))
  }
  
  pie_df <- t(as.matrix(expr[c(gene_x, gene_y), , drop = FALSE])) |> as.data.frame()
  colnames(pie_df) <- c(gene_x, gene_y)
  if (!is.null(colnames(expr))) rownames(pie_df) <- colnames(expr)
  
  pie_df$CoExp <- dplyr::case_when(
    pie_df[[gene_x]] > threshold & pie_df[[gene_y]] > threshold ~ "Both",
    pie_df[[gene_x]] > threshold & pie_df[[gene_y]] <= threshold ~ gene_x,
    pie_df[[gene_x]] <= threshold & pie_df[[gene_y]] > threshold ~ gene_y,
    TRUE ~ "None"
  )
  lvl <- c("Both", gene_x, gene_y, "None")
  pie_df$CoExp <- factor(pie_df$CoExp, levels = lvl)
  
  plot_df <- dplyr::count(pie_df, CoExp, .drop = FALSE) |>
    dplyr::mutate(percent = n / sum(n) * 100)
  
  pal <- c("Both" = both_color, setNames(gene_x_color, gene_x),
           setNames(gene_y_color, gene_y), "None" = none_color)
  
  if (is.null(legend_title)) legend_title <- paste0(gene_x, "/", gene_y)
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = "", y = percent, fill = CoExp)) +
    ggplot2::geom_col(width = 1, color = "white") +
    ggplot2::coord_polar(theta = "y", clip = "off") +                 # <- no clipping in panel
    ggplot2::scale_fill_manual(values = pal, limits = lvl, drop = FALSE, name = legend_title) +
    ggplot2::theme_void(base_size = 14) +
    ggplot2::theme(
      legend.position   = legend_side,
      legend.title      = ggplot2::element_text(size = 20, face = "bold", color = "black"),
      legend.text       = ggplot2::element_text(size = 20, color = "black"),
      legend.key.size   = grid::unit(1.0, "cm"),
      plot.margin       = ggplot2::margin(margins_mm, margins_mm, margins_mm, margins_mm)
    )
  
  to_save <- p
  if (lock_panel && requireNamespace("egg", quietly = TRUE)) {
    to_save <- egg::set_panel_size(
      to_save,
      width  = grid::unit(panel_mm, "mm"),
      height = grid::unit(panel_mm, "mm")
    )
  }
  
  # Auto-compute canvas size if not supplied
  if (is.null(width_mm) || is.null(height_mm)) {
    if (legend_side == "right") {
      if (is.null(width_mm))  width_mm  <- panel_mm + legend_space_mm + 2 * margins_mm
      if (is.null(height_mm)) height_mm <- panel_mm + 2 * margins_mm
    } else { # legend at bottom
      if (is.null(width_mm))  width_mm  <- panel_mm + 2 * margins_mm
      if (is.null(height_mm)) height_mm <- panel_mm + legend_space_mm + 2 * margins_mm
    }
  }
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(out_basename)) out_basename <- sprintf("pie_%s_%s.tiff", gene_x, gene_y)
  out_path <- file.path(out_dir, out_basename)
  
  ggplot2::ggsave(out_path, plot = to_save, device = "tiff",
                  width = width_mm, height = height_mm, units = "mm", dpi = dpi)
  
  invisible(list(file = out_path, plot = p, locked_plot = to_save, counts = plot_df))
}





jhu_pdx_data <- read_xlsx(path = "/projects/eshenderov-hpc/Shivang/projects/B7-H3 ADC Manuscript/data/JHU_PDX_PCa.xlsx", col_names = FALSE)
jhu_repeat_names <-jhu_pdx_data[3, c(2:dim(jhu_pdx_data)[2])] %>% as.character()
jhu_sample_names <- jhu_repeat_names
jhu_genes <- jhu_pdx_data[c(4:dim(jhu_pdx_data)[1]), 1][["...1"]]
jhu_expr <- jhu_pdx_data[c(4:dim(jhu_pdx_data)[1]), c(2:dim(jhu_pdx_data)[2])] %>% as.data.frame()
jhu_expr <- jhu_expr %>% mutate_if(is.character, as.numeric)
colnames(jhu_expr) <- jhu_repeat_names
jhu_expr <- cbind(jhu_genes, jhu_expr)
colnames(jhu_expr)[1] <- "id"
jhu_expr <- combine_duplicate_genes(jhu_expr) %>% as.data.frame()
rownames(jhu_expr) <- jhu_expr$gene
jhu_expr$gene <- NULL
jhu_expr_norm <- fpkm_to_logTPM(jhu_expr)


uw_pdx_data <- read_xlsx(path = "/projects/eshenderov-hpc/Shivang/projects/B7-H3 ADC Manuscript/data/UW_PDX_PCa.xlsx", col_names = FALSE)
uw_repeat_names <- uw_pdx_data[3, c(2:dim(uw_pdx_data)[2])] %>% as.character()
uw_sample_names <- uw_pdx_data[2, c(2:dim(uw_pdx_data)[2])] %>% as.character()
uw_genes <- uw_pdx_data[c(4:dim(uw_pdx_data)[1]), 1][["...1"]]
uw_expr <- uw_pdx_data[c(4:dim(uw_pdx_data)[1]), c(2:dim(uw_pdx_data)[2])] %>% as.data.frame()
uw_expr <- uw_expr %>% mutate_if(is.character, as.numeric)
colnames(uw_expr) <- uw_repeat_names
uw_expr <- cbind(uw_genes, uw_expr)
colnames(uw_expr)[1] <- "id"
uw_expr <- combine_duplicate_genes(uw_expr) %>% as.data.frame()
rownames(uw_expr) <- uw_expr$gene
uw_expr$gene <- NULL
uw_expr_norm <- fpkm_to_logTPM(uw_expr)

annotations <- read_xlsx(path = "/projects/eshenderov-hpc/Shivang/projects/B7-H3 ADC Manuscript/data/combined_meta_data_PDX.xlsx", col_names = FALSE) %>% distinct()
annotations <- annotations[c(1, 2), c(3:dim(annotations)[2])] %>% t() %>% as.data.frame() %>% distinct()
anno_dict <- annotations$V1
names(anno_dict) <- annotations$V2

meta <- data.frame(Repeat_Names = c(jhu_repeat_names, uw_repeat_names), 
                   Sample_Names = c(jhu_sample_names, uw_sample_names), 
                   Subtype = anno_dict[c(jhu_sample_names, uw_sample_names)])
rownames(meta) <- meta$Repeat_Names
meta$batch <- c(rep("JHU", times = length(jhu_repeat_names)), rep("UW", times = length(uw_repeat_names)))
meta$Subtype <- anno_dict[meta$Sample_Names]

merged <- merge(jhu_expr_norm, uw_expr_norm, by = "row.names", all = TRUE)
merged[is.na(merged)] <- 0
rownames(merged) <- merged$Row.names
merged$Row.names <- NULL
expr_mat <- as.matrix(sapply(merged, as.numeric))
rownames(expr_mat) <- rownames(merged)
RLE_pre_batch_correction <- RLE_plot(as.matrix(expr_mat))
ggsave(file.path(file.path(project_dir, "figures", "PDX_Supp_Data"), paste0("RLE_pre_batch_correction", ".tiff")),
       plot = RLE_pre_batch_correction, device = "tiff",
       width = 360, height = 120, units = "mm", dpi = 600)

meta <- meta[colnames(expr_mat), ]
batch <- factor(meta$batch)  # batch variable
expr_corrected <- removeBatchEffect(expr_mat, batch = batch)
RLE_post_batch_correction <- RLE_plot(as.matrix(expr_corrected))
ggsave(file.path(file.path(project_dir, "figures", "PDX_Supp_Data"), paste0("RLE_post_batch_correction", ".tiff")),
       plot = RLE_post_batch_correction, device = "tiff",
       width = 360, height = 120, units = "mm", dpi = 600)

intersecting_genes <- read_xlsx("/projects/eshenderov-hpc/Shivang/projects/B7-H3 ADC Manuscript/data/Hu_Genes.xlsx")
intersecting_genes <- intersect(intersecting_genes$JHU_Hu_Genes, intersecting_genes$UW_Hu_Genes)

expr <- as.data.frame(expr_corrected[intersecting_genes, ])


# Rearrange samples in expr and meta -------------------------------------------
# your global subtype order
desired <- c("DNPC","ARneg_NEneg","ARneg_NElow","ARneg_NEpos",
             "ARpos_NEpos","AMPHICRINE","ARlow_NEneg","ARpos_NEneg","ARPC")

joint <- c("ARpos_NEneg","ARPC")  # treat as one block

meta <- meta %>%
  mutate(
    Subtype = factor(Subtype, levels = desired),
    # ends with CR / CRdigits / R / Rdigits
    cr_or_r = grepl("(CR\\d*|R\\d*)$", Repeat_Names),
    
    # rank for the big blocks; joint subtypes share the same rank
    block_rank = if_else(Subtype %in% joint,
                         match("ARpos_NEneg", desired),      # use the earlier of the two
                         match(as.character(Subtype), desired)),
    
    # within the joint block: CR/R first (0), then non-CR/R (1)
    joint_priority = case_when(
      Subtype %in% joint & cr_or_r ~ 0L,
      Subtype %in% joint & !cr_or_r ~ 1L,
      TRUE ~ 0L  # irrelevant for non-joint blocks
    ),
    
    .row_id = row_number()  # to preserve original order elsewhere
  ) %>%
  arrange(block_rank, joint_priority, .row_id) %>%
  dplyr::select(-block_rank, -joint_priority, -.row_id)

orig_subtypes <- c("ARPC","AMPHICRINE","DNPC","ARpos_NEpos","ARpos_NEneg",
                   "ARneg_NEpos","ARneg_NEneg","ARneg_NElow","ARlow_NEneg")
new_subtypes  <- c("AR+/NE-","AR+/NE+","AR-/NE-","AR+/NE+","AR+/NE-",
                   "AR-/NE+","AR-/NE-","AR-/NE_l","AR_l/NE-")
map <- setNames(new_subtypes, orig_subtypes)

meta <- meta %>%
  mutate(
    Subtype_old = Subtype,
    Subtype = recode(as.character(Subtype), !!!map)  # no .default needed
  )

# (optional) set plotting order for the new labels
desired_orig <- c("DNPC","ARneg_NEneg","ARneg_NElow","ARneg_NEpos",
                  "ARpos_NEpos","AMPHICRINE","ARlow_NEneg","ARpos_NEneg","ARPC")
desired_new  <- unique(map[desired_orig])
meta$Subtype <- factor(meta$Subtype, levels = desired_new)

expr <- expr[, meta$Repeat_Names]

# ------------------------------------------------------------------------------


# Correlation plot and horizontal strips with hclust on genes
genes <- c("STEAP1", "STEAP2",  "PSCA", "FOLH1", "KLK2", "TACSTD2", "NECTIN1", "NECTIN4", "DLL3", "CD276", "AR", "KLK3", "TMPRSS2", "FKBP5", "IGF1R")
corr_ht <- corr_heatmap_genes(expr, genes, project_dir, method = "pearson",
                              file = "gene_correlation_heatmap.pdf", width = 7, height = 7)
ht_drawn <- ComplexHeatmap::draw(corr_ht$heatmap)
library(ggplotify)   # install.packages("ggplotify")

ht <- corr_ht$heatmap
p <- as.ggplot(function() ComplexHeatmap::draw(ht))

ggsave(
  file.path(project_dir, "figures", "PDX_Supp_Data", "gene_correlation_plot.tiff"),
  plot = p,
  device = "tiff",
  width = 120, height = 120, units = "mm", dpi = 600, compression = "lzw"
)


# Expression variation strip
res_strip_PSCA <- expr_line_strips(expr,
                                   genes = c("PSCA", "CD276", "NECTIN1"), 
                                   project_dir, 
                                   order = "heatmap",          # or "heatmap" with heatmap_ht=...
                                   heatmap_ht = ht_drawn, 
                                   cluster_dist = c("euclidean"), 
                                   scale_method = "zscore",    # comparable y-ranges
                                   zlim = c(-3, 3),           # same y-limits for all strips
                                   point_size = 0.6, line_size = 0.3, 
                                   show_x_labels = FALSE, 
                                   plot_name = "FigPDX_b", 
                                   gene_names = c("PSCA", "B7-H3", "NECTIN1")
)
print(res_strip_PSCA$plot)
ggsave(file.path(file.path(project_dir, "figures", "PDX_Supp_Data"), paste0("PSCA_CD276_NECTIN1_expr_var_strip", ".tiff")),
       plot = res_strip_PSCA$plot, device = "tiff",
       width = 600, height = 120, units = "mm", dpi = 600)


# Expression variation strip
res_strip_PSCA <- expr_line_strips(expr,
                                   genes = c("CD276", "PSCA"), 
                                   project_dir, 
                                   order = "heatmap",          # or "heatmap" with heatmap_ht=...
                                   heatmap_ht = ht_drawn, 
                                   cluster_dist = c("euclidean"), 
                                   scale_method = "zscore",    # comparable y-ranges
                                   zlim = c(-3, 3),           # same y-limits for all strips
                                   point_size = 0.6, line_size = 0.3, 
                                   show_x_labels = FALSE, 
                                   plot_name = "FigPDX_c", 
                                   gene_names = c("B7-H3", "PSCA")
)
print(res_strip_PSCA$plot)
ggsave(file.path(file.path(project_dir, "figures", "PDX_Supp_Data"), paste0("PSCA_CD276_expr_var_strip", ".tiff")),
       plot = res_strip_PSCA$plot, device = "tiff",
       width = 600, height = 120, units = "mm", dpi = 600)



expr_der <- gene_row_derivative(expr, genes)
corr_ht_der <- corr_heatmap_genes(expr_der, genes, project_dir = project_dir, method = "pearson",
                                  file = "gene_derivative_correlation_heatmap.pdf", width = 7, height = 7)
ht_drawn_der <- ComplexHeatmap::draw(corr_ht$heatmap)
p <- as.ggplot(function() ComplexHeatmap::draw(ht_drawn_der))

ggsave(
  file.path(project_dir, "figures", "PDX_Supp_Data", "gene_derivative_correlation_heatmap.tiff"),
  plot = p,
  device = "tiff",
  width = 120, height = 120, units = "mm", dpi = 600, compression = "lzw"
)



#######################################################
# -------------------- Main figures ------------------#
#######################################################

# -------------------- Figure a -----------------------
genes <- c("STEAP1", "STEAP2",  "PSCA", "FOLH1", "KLK2", "TACSTD2", "NECTIN1", "NECTIN4", "DLL3", "CD276", "AR", "KLK3", "TMPRSS2", "FKBP5", "IGF1R")
gene_names <- c("STEAP1", "STEAP2", "PSCA", "PSMA", "KLK2", "TROP-2", "NECTIN1", "NECTIN4", "DLL3", "B7-H3", "AR", "KLK3", "TMPRSS2", "FKBP5", "IGF1R")
heatmap_row_split <- c(10, 5)

expr_ht <- expression_map(expr, meta, genes, gene_names, heatmap_row_split, project_dir, file_name = "figures/PDX_PSCA_NECTIN1/FigPDX_a")


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# ---- Required packages ----
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(grDevices)
})

# ---- Inputs expected ----
# expr:           numeric matrix (rows = genes, cols = samples)
# genes:          character vector of row names to plot (in desired order)
# gene_names:     row label vector, same length/order as 'genes'
# heatmap_row_split: integer vector of lengths for c("ADC","Immune","AR reg.")
# project_dir:    base output directory (e.g., "results")
# file_name:      relative path under project_dir without extension
#                  e.g., "figures/PDX_PSCA_NECTIN1/FigPDX_a"

# -------- Align metadata to expression columns --------
if (!is.null(rownames(meta)) && all(colnames(expr) %in% rownames(meta))) {
  meta_use <- meta[colnames(expr), , drop = FALSE]
} else if ("Sample" %in% colnames(meta) && all(colnames(expr) %in% meta$Sample)) {
  meta_use <- meta[match(colnames(expr), meta$Sample), , drop = FALSE]
  rownames(meta_use) <- meta_use$Sample
} else {
  stop("Cannot align `meta` to `expr` columns. Ensure rownames(meta) or meta$Sample match colnames(expr).")
}

stopifnot("Subtype" %in% colnames(meta_use))

# Optional: enforce a canonical subtype order if present
canonical_order <- c("AR-/NE-", "AR-/NE_l", "AR-/NE+", "AR+/NE+", "AR_l/NE-", "AR+/NE-")
subtypes_present <- unique(meta_use$Subtype)
if (all(subtypes_present %in% canonical_order)) {
  meta_use$Subtype <- factor(meta_use$Subtype, levels = canonical_order)
} else {
  meta_use$Subtype <- factor(meta_use$Subtype)  # keep observed order
}

# -------- Subset genes --------
expr_sub <- expr[genes, , drop = FALSE]


# -------- Compute dispersion (sd / mean) per Subtype for each gene --------
# guard against division-by-zero and NAs
eps <- 1e-12
st_levels <- levels(meta_use$Subtype)
disp_mat <- sapply(st_levels, function(st) {
  idx <- which(meta_use$Subtype == st)
  if (length(idx) == 0) return(rep(NA_real_, length(genes)))
  m <- rowMeans(expr_sub[, idx, drop = FALSE], na.rm = TRUE)
  s <- apply(expr_sub[, idx, drop = FALSE], 1, sd, na.rm = TRUE)
  d <- s / pmax(abs(m), eps)
  # replace non-finite with NA
  d[!is.finite(d)] <- NA_real_
  d
}, simplify = TRUE)
rownames(disp_mat) <- genes
colnames(disp_mat) <- st_levels

# -------- Optionally cap extreme values to 99th percentile for visualization --------
finite_vals <- as.numeric(disp_mat[is.finite(disp_mat)])
if (length(finite_vals) >= 10) {
  cap <- quantile(finite_vals, 0.99, na.rm = TRUE)
  disp_mat[disp_mat > cap] <- cap
}

# -------- Row grouping (ADC | Immune | AR reg.) --------
row_groups <- rep(c("ADC","AR reg."), times = heatmap_row_split)
stopifnot(length(row_groups) == nrow(disp_mat))

# -------- Row annotation --------
row_group_cols <- c("ADC"="#FF0000","AR reg."="#0092FF")
row_anno <- rowAnnotation(
  Gene = factor(row_groups, levels = c("ADC","AR reg.")),
  col  = list(Gene = row_group_cols),
  annotation_legend_param = list(Gene = list(title = "Group"))
)

# -------- Heatmap color function (robust to NA / narrow range) --------
rng <- range(disp_mat, na.rm = TRUE, finite = TRUE)
if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) rng <- rng + c(-1e-6, 1e-6)
col_fun <- circlize::colorRamp2(
  seq(rng[1], rng[2], length.out = 256),
  hcl.colors(256, "Viridis", rev = TRUE)
)

# -------- Heatmap (summary statistic = dispersion) --------
ht <- Heatmap(
  matrix = disp_mat,
  name = "Dispersion\n(sd/mean)",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = factor(row_groups, levels = c("ADC","AR reg.")),
  left_annotation = row_anno,
  row_names_side = "left",
  row_labels = gene_names,
  show_column_names = TRUE,
  row_title = NULL,
  border = TRUE,
  row_gap = unit(2, "mm"),
  width  = (ncol(disp_mat) * unit(6, "mm")) * 1.1,
  height = (nrow(disp_mat) * unit(3.5, "mm")) * 1.2,
  row_names_gp = gpar(fontsize = 10)
)

file_name <- "dispersion_heatmap_PDX"

# -------- Output sizes --------
w_in <- as.numeric(convertUnit((ncol(disp_mat) * unit(6, "mm")) * 2.5, "in", valueOnly = TRUE))
h_in <- as.numeric(convertUnit((nrow(disp_mat) * unit(3.5, "mm")) * 2, "in", valueOnly = TRUE))

# -------- Ensure output directory --------
out_base <- file.path(project_dir, "figures")
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

# -------- Save PNG --------
png(file = file.path(project_dir, "figures", paste0(file_name, ".png")),
    width = w_in, height = h_in, units = "in", res = 300)
ComplexHeatmap::draw(ht)
dev.off()

# -------- Save PDF --------
pdf(file = file.path(project_dir, "figures", paste0(file_name, ".pdf")),
    width = w_in, height = h_in, onefile = FALSE, useDingbats = FALSE)
ComplexHeatmap::draw(ht)
dev.off()

# -------- Save TIFF (high-res) --------
dpi <- 600
tiff(file = file.path(project_dir, "figures", paste0(file_name, ".tif")),
     width = w_in, height = h_in, units = "in", res = dpi, compression = "lzw")
ComplexHeatmap::draw(ht)
dev.off()


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

ar_pos_cat <- c("AR+/NE-", "AR+/NE+", "AR+/NE+", "AR+/NE-", "AR_l/NE-")
ar_neg_cat <- c("AR-/NE-","AR-/NE+", "AR-/NE-", "AR-/NE_l")
ar_status <- c(rep("AR+", times = length(ar_pos_cat)), rep("AR-", times = length(ar_neg_cat)))
names(ar_status) <- c(ar_pos_cat, ar_neg_cat)
meta$`AR Status` <- ar_status[meta$Subtype]

intact_castration_annotation <- read_xlsx(path = "/projects/eshenderov-hpc/Shivang/projects/B7-H3 ADC Manuscript/data/Intact_vs_Castration_anntoation.xlsx", col_names = TRUE)
intact_samples <- intact_castration_annotation$Intact[!is.na(intact_castration_annotation$Intact)]
castration_samples <- intact_castration_annotation$Castration[!is.na(intact_castration_annotation$Castration)]
castration_status <- c(rep("Intact", times = length(intact_samples)), rep("Castrated", times = length(castration_samples)))
names(castration_status) <- c(intact_samples, castration_samples)
meta$`Castration Status` <- castration_status[meta$Sample_Names]

meta$`B7-H3` <- expr["CD276", ] %>% as.numeric()
meta$PSCA <- expr["PSCA", ] %>% as.numeric()
meta$NECTIN1 <- expr["NECTIN1", ] %>% as.numeric()


# ----- 1) Tidy your columns (handle spaces/hyphen) -----
# df = your data frame shown above
df_clean <- meta %>%
  mutate(
    `B7-H3`           = as.numeric(`B7-H3`), 
    PSCA              = as.numeric(PSCA),
    NECTIN1           = as.numeric(NECTIN1),
    AR_Status         = factor(`AR Status`, levels = c("AR+","AR-")),
    Castration_Status = `Castration Status`,
    # patient id = everything before the final underscore in Sample_Names
    patient           = Sample_Names
  )

# -----------------------------------------------------------------------------
# Plot A: AR+ vs AR- (all rows used; no NA issue here)
# -----------------------------------------------------------------------------

# -------------------- Figure b, c, d -----------------------
res_ar <- plot_ar_status(
  df = df_clean,                # or your raw df; just ensure the column names match
  project_dir = project_dir,    # path to your project
  status_col = "AR Status",
  y_col = "B7-H3",
  method = "wilcox.test",            # or "wilcox.test"
  file_base = "PDX_PSCA_NECTIN1/Fig2_b"
)
res_ar$plot

res_ar <- plot_ar_status(
  df = df_clean,                # or your raw df; just ensure the column names match
  project_dir = project_dir,    # path to your project
  status_col = "AR Status",
  y_col = "PSCA",
  method = "wilcox.test",            # or "wilcox.test"
  file_base = "PDX_PSCA_NECTIN1/Fig2_c"
)
res_ar$plot

res_ar <- plot_ar_status(
  df = df_clean,                # or your raw df; just ensure the column names match
  project_dir = project_dir,    # path to your project
  status_col = "AR Status",
  y_col = "NECTIN1",
  method = "wilcox.test",            # or "wilcox.test"
  file_base = "PDX_PSCA_NECTIN1/Fig2_d"
)
res_ar$plot


# -------------------- Figure e, f, g -----------------------
res_cast <- plot_intact_vs_castrated(df_clean, project_dir, sample_col = "Sample_Names", y_col = "B7-H3",
                                    test = "wilcox.test", file_name = "PDX_PSCA_NECTIN1/Fig2_e")
res_cast$plot
res_cast <- plot_intact_vs_castrated(df_clean, project_dir, sample_col = "Sample_Names", y_col = "PSCA",
                                     test = "wilcox.test", file_name = "PDX_PSCA_NECTIN1/Fig2_f")
res_cast$plot
res_cast <- plot_intact_vs_castrated(df_clean, project_dir, sample_col = "Sample_Names", y_col = "NECTIN1",
                                     test = "wilcox.test", file_name = "PDX_PSCA_NECTIN1/Fig2_g")
res_cast$plot

# (Optional) Save
# ggsave("ARpos_vs_ARneg_B7H3.pdf", p_ar, width = 3.5, height = 4.5)
# ggsave("Intact_vs_Castrate_B7H3_paired.pdf", p_cast, width = 3.5, height = 5)

# -------------------- Figure h, i -----------------------
pie_df <- merged[c("PSCA", "CD276", "NECTIN1"), ] %>% 
  as.matrix() %>% 
  t() %>% 
  as.data.frame()
colnames(pie_df) <- c("PSCA", "CD276", "NECTIN1")
rownames(pie_df) <- colnames(merged)

pie_df$CD276_NECTIN1_CoExp <- dplyr::case_when(
  pie_df$CD276 > 1 & pie_df$NECTIN1 > 1 ~ "Both",
  pie_df$CD276 > 1 & pie_df$NECTIN1 <= 1 ~ "CD276",
  pie_df$CD276 <= 1 & pie_df$NECTIN1 > 1 ~ "NECTIN1",
  TRUE ~ "None"
)

pie_df$CD276_PSCA_CoExp <- dplyr::case_when(
  pie_df$CD276 > 1 & pie_df$PSCA > 1 ~ "Both",
  pie_df$CD276 > 1 & pie_df$PSCA <= 1 ~ "CD276",
  pie_df$CD276 <= 1 & pie_df$PSCA > 1 ~ "PSCA",
  TRUE ~ "None"
)

# Set factor levels in desired order
pie_df$CD276_NECTIN1_CoExp <- factor(
  pie_df$CD276_NECTIN1_CoExp,
  levels = c("Both", "CD276", "NECTIN1", "None")
)

# Set factor levels in desired order
pie_df$CD276_PSCA_CoExp <- factor(
  pie_df$CD276_PSCA_CoExp,
  levels = c("Both", "CD276", "PSCA", "None")
)

pie_df <- as.data.frame(pie_df)

# Count proportions
plot_CD276_NECTIN1_df <- pie_df %>%
  dplyr::count(CD276_NECTIN1_CoExp) %>%
  mutate(percent = n / sum(n) * 100)

plot_CD276_PSCA_df <- pie_df %>%
  dplyr::count(CD276_PSCA_CoExp) %>%
  mutate(percent = n / sum(n) * 100)

# Pie chart

# Pie Chart colors
pie_colors <- c(
  "Both"    = "#9E6A9D",  # red
  "CD276"   = "#3F7FBF",  # blue
  "NECTIN1" = "#2FA37B",  # green
  "None"    = "#B5B5B5"   # gray
)

# assume plot_CD276_NECTIN1_df already completed to include all 4 categories
df <- plot_CD276_NECTIN1_df %>%
  mutate(CD276_NECTIN1_CoExp = factor(CD276_NECTIN1_CoExp, levels = names(pie_colors))) %>%
  complete(CD276_NECTIN1_CoExp = names(pie_colors), fill = list(percent = 0)) %>%
  arrange(CD276_NECTIN1_CoExp) %>%
  mutate(
    ymax = cumsum(percent),
    ymin = lag(ymax, default = 0),
    ymid = (ymin + ymax) / 2,
    lab  = ifelse(percent > 0, paste0(CD276_NECTIN1_CoExp, " (", round(percent, 1), "%)"), "")
  )

p_CD276_NECTIN1 <- ggplot(df, aes(x = 1, y = percent, fill = CD276_NECTIN1_CoExp)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y", clip = "off") +   # allow labels outside
  # Less horizontal span -> bigger pie; still some room for labels
  xlim(c(0.5, 1.5)) +
  scale_fill_manual(
    name   = "Co-expression",
    values = pie_colors,
    breaks = names(pie_colors),
    drop   = FALSE
  ) +
  theme_void(base_size = 14) +
  theme(
    legend.position  = "right",
    legend.title     = element_text(size = 16, face = "bold"),
    legend.text      = element_text(size = 14),
    legend.key.size  = unit(1.2, "cm"),
    plot.margin      = margin(10, 40, 10, 10)  # extra right margin for longer labels
  )
p_CD276_NECTIN1

ggsave(file.path(file.path(project_dir, "figures"), paste0("PDX_PSCA_NECTIN1/FigPDX_h", ".tiff")),
       plot = p_CD276_NECTIN1, device = "tiff",
       width = 120, height = 120, units = "mm", dpi = 600)
ggsave(file.path(file.path(project_dir, "figures"), paste0("PDX_PSCA_NECTIN1/FigPDX_h", ".pdf")),
       plot = p_CD276_NECTIN1, device = "pdf",
       width = 120, height = 120, units = "mm", dpi = 600)

pie_colors <- c(
  "Both"    = "#9E6A9D",  # red
  "CD276"   = "#3F7FBF",  # blue
  "PSCA"    = "#2FA37B",  # green
  "None"    = "#B5B5B5"   # gray
)

# Pie chart
df <- plot_CD276_PSCA_df %>%
  mutate(CD276_PSCA_CoExp = factor(CD276_PSCA_CoExp, levels = names(pie_colors))) %>%
  complete(CD276_PSCA_CoExp = names(pie_colors), fill = list(percent = 0)) %>%
  arrange(CD276_PSCA_CoExp) %>%
  mutate(
    ymax = cumsum(percent),
    ymin = lag(ymax, default = 0),
    ymid = (ymin + ymax) / 2,
    lab  = ifelse(percent > 0, paste0(CD276_PSCA_CoExp, " (", round(percent, 1), "%)"), "")
  )

p_CD276_PSCA <- ggplot(df, aes(x = 1, y = percent, fill = CD276_PSCA_CoExp)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y", clip = "off") +   # allow labels outside
  # Less horizontal span -> bigger pie; still some room for labels
  xlim(c(0.5, 1.5)) +
  scale_fill_manual(
    name   = "Co-expression",
    values = pie_colors,
    breaks = names(pie_colors),
    drop   = FALSE
  ) +
  theme_void(base_size = 14) +
  theme(
    legend.position  = "right",
    legend.title     = element_text(size = 16, face = "bold"),
    legend.text      = element_text(size = 14),
    legend.key.size  = unit(1.2, "cm"),
    plot.margin      = margin(10, 40, 10, 10)  # extra right margin for longer labels
  )
p_CD276_PSCA

ggsave(file.path(file.path(project_dir, "figures"), paste0("PDX_PSCA_NECTIN1/FigPDX_i", ".tiff")),
       plot = p_CD276_PSCA, device = "tiff",
       width = 120, height = 120, units = "mm", dpi = 600)
ggsave(file.path(file.path(project_dir, "figures"), paste0("PDX_PSCA_NECTIN1/FigPDX_i", ".pdf")),
       plot = p_CD276_PSCA, device = "pdf",
       width = 120, height = 120, units = "mm", dpi = 600)

# ------------------- Main text plots end here ---------------------


make_coexp_pie(
  expr = merged, gene_x = "CD276", gene_y = "NECTIN1",
  out_dir = file.path(project_dir, "figures"),
  out_basename = "FigPDX_c.tiff"
)

# PSCA pair (saves like your FigPDX_b.tiff)
make_coexp_pie(
  expr = merged, gene_x = "CD276", gene_y = "PSCA",
  out_dir = file.path(project_dir, "figures"),
  out_basename = "FigPDX_b.tiff"
)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
#                        Plots for CD276, NECTIN1 & STEAP2                      
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------


pie_df <- merged[c("STEAP2", "CD276", "NECTIN1"), ] %>% 
  as.matrix() %>% 
  t() %>% 
  as.data.frame()
colnames(pie_df) <- c("STEAP2", "CD276", "NECTIN1")
rownames(pie_df) <- colnames(merged)

pie_df$CD276_NECTIN1_CoExp <- dplyr::case_when(
  pie_df$CD276 > 1 & pie_df$NECTIN1 > 1 ~ "Both",
  pie_df$CD276 > 1 & pie_df$NECTIN1 <= 1 ~ "CD276",
  pie_df$CD276 <= 1 & pie_df$NECTIN1 > 1 ~ "NECTIN1",
  TRUE ~ "None"
)

pie_df$CD276_STEAP2_CoExp <- dplyr::case_when(
  pie_df$CD276 > 1 & pie_df$STEAP2 > 1 ~ "Both",
  pie_df$CD276 > 1 & pie_df$STEAP2 <= 1 ~ "CD276",
  pie_df$CD276 <= 1 & pie_df$STEAP2 > 1 ~ "STEAP2",
  TRUE ~ "None"
)

# Set factor levels in desired order
pie_df$CD276_NECTIN1_CoExp <- factor(
  pie_df$CD276_NECTIN1_CoExp,
  levels = c("Both", "CD276", "NECTIN1", "None")
)

# Set factor levels in desired order
pie_df$CD276_STEAP2_CoExp <- factor(
  pie_df$CD276_STEAP2_CoExp,
  levels = c("Both", "CD276", "STEAP2", "None")
)

pie_df <- as.data.frame(pie_df)

# Count proportions
plot_CD276_NECTIN1_df <- pie_df %>%
  dplyr::count(CD276_NECTIN1_CoExp) %>%
  mutate(percent = n / sum(n) * 100)

plot_CD276_STEAP2_df <- pie_df %>%
  dplyr::count(CD276_STEAP2_CoExp) %>%
  mutate(percent = n / sum(n) * 100)

# Pie chart
p_CD276_NECTIN1 <- ggplot(plot_CD276_NECTIN1_df, aes(x = "", y = percent, fill = CD276_NECTIN1_CoExp)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(
    values = c(
      "Both"    = "#E64B35FF",  # red
      "CD276"   = "#4DBBD5FF",  # blue
      "NECTIN1" = "#00A087FF",  # green
      "None"    = "#999999FF"   # gray
    )
  ) +
  theme_void(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text  = element_text(size = 14),
    legend.key.size = unit(1.2, "cm"),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(file.path(file.path(project_dir, "figures"), paste0("FigPDX_c", ".tiff")),
       plot = p_CD276_NECTIN1, device = "tiff",
       width = 120, height = 120, units = "mm", dpi = 600)

# Pie chart
p_CD276_PSCA <- ggplot(plot_CD276_STEAP2_df, aes(x = "", y = percent, fill = CD276_STEAP2_CoExp)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(
    values = c(
      "Both"    = "#E64B35FF",  # red
      "CD276"   = "#4DBBD5FF",  # blue
      "STEAP2" = "#00A087FF",  # green
      "None"    = "#999999FF"   # gray
    )
  ) +
  theme_void(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text  = element_text(size = 14),
    legend.key.size = unit(1.2, "cm"),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(file.path(file.path(project_dir, "figures"), paste0("FigPDX_b", ".tiff")),
       plot = p_CD276_PSCA, device = "tiff",
       width = 120, height = 120, units = "mm", dpi = 600)


make_coexp_pie(
  expr = merged, gene_x = "CD276", gene_y = "NECTIN1",
  out_dir = file.path(project_dir, "figures"),
  out_basename = "FigPDX_c.tiff"
)

# PSCA pair (saves like your FigPDX_b.tiff)
make_coexp_pie(
  expr = merged, gene_x = "CD276", gene_y = "STEAP2",
  out_dir = file.path(project_dir, "figures"),
  out_basename = "FigPDX_b.tiff"
)

# -----------------------------------------------------------------------------
# Plot A: AR+ vs AR- (all rows used; no NA issue here)
# -----------------------------------------------------------------------------

meta$STEAP2 <- expr["STEAP2", ] %>% as.numeric()

# ----- 1) Tidy your columns (handle spaces/hyphen) -----
# df = your data frame shown above
df_clean <- meta %>%
  mutate(
    B7_H3             = as.numeric(`B7-H3`), 
    STEAP2              = as.numeric(STEAP2),
    AR_Status         = factor(`AR Status`, levels = c("AR+","AR-")),
    Castration_Status = `Castration Status`,
    # patient id = everything before the final underscore in Sample_Names
    patient           = Sample_Names
  )


res_ar <- plot_ar_status(
  df = df_clean,                # or your raw df; just ensure the column names match
  project_dir = project_dir,    # path to your project
  status_col = "AR Status",
  y_col = "B7-H3",
  method = "t.test",            # or "wilcox.test"
  file_base = "FigPDX_d"
)
res_ar <- plot_ar_status(
  df = df_clean,                # or your raw df; just ensure the column names match
  project_dir = project_dir,    # path to your project
  status_col = "AR Status",
  y_col = "STEAP2",
  method = "t.test",            # or "wilcox.test"
  file_base = "FigPDX_e"
)

res_ar$plot


res_cast <- plot_intact_vs_castrated(df_clean, project_dir, sample_col = "Sample_Names", y_col = "B7-H3",
                                     test = "t.test", file_name = "FigPDX_f.tiff")
res_cast <- plot_intact_vs_castrated(df_clean, project_dir, sample_col = "Sample_Names", y_col = "STEAP2",
                                     test = "t.test", file_name = "FigPDX_g.tiff")
res_cast$plot