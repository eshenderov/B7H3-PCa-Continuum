fig2a <- function (obj, reduction, project_dir, genes, gene_names) {
  
  features <- genes
  features_title <- gene_names
  axes <- list(x = -12, y = -15, x_len = 5, y_len = 5)
  stallion = c("1"="#D4477D","2"="#272E6A","3"="#208A42","4"="#89288F", 
               "5"="#F47D2B", "6"="#FEE500","7"="#8A9FD1","8"="#C06CAB", 
               "19"="#E6C2DC", "10"="#90D5E4", "11"="#89C75F","12"="#F37B7D",
               "13"="#9983BD","14"="#D24B27","15"="#3BBCA8", "16"="#6E4B9E",
               "17"="#0C727C", "18"="#7E1416","9"="#D8A767","20"="#3D3D3D"
  )
  stallion_cols <- c("#D4477D", "#2F8AC4", "#7DD06F", "#89288F", "#FEE500", "#8A9FD1", "#F47D2B", "#7B6FD0", "#90D5E4", "#60824f", "#88CFA4", "#D19EC4", "#D8A767", "#D0CD47", "#484125")
  names(stallion_cols) <- levels(Idents(obj_sub))
  umap_df <- cbind(obj_sub@reductions[[reduction]]@cell.embeddings, 
                   data.frame(Label = obj_sub@meta.data[["CoarseAnnotation"]])
  )
  colnames(umap_df) <- c("UMAP1", "UMAP2", "Label")
  
  centroids <- umap_df %>%
    group_by(Label) %>%
    summarize(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2), .groups = 'drop')
  
  # UMAP with legend and Labels
  umap_plot_label_legend <- ggplot(data = umap_df, mapping = aes(x = UMAP1, y = UMAP2, color = Label)) + 
    geom_point(alpha = 0.1, size = 0.0001) + 
    scale_color_manual(values = stallion_cols) + 
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      axis.ticks = element_blank(),  # Remove axis ticks
      axis.title = element_blank(),  # Increase axis title size
      axis.text = element_blank(),   # Increase axis text size
      legend.title = element_text(size = 14),  # Increase legend title size
      legend.text = element_text(size = 11)    # Increase legend text size
    ) + 
    geom_text_repel(data = centroids, 
                    aes(label = Label), 
                    size = 3, 
                    nudge_y = 0.1, 
                    color = "black") + 
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +  # Increase dot size in legend
    annotate("segment", 
             x = axes$x, xend = axes$x + c(axes$x_len, 0), 
             y = axes$y, yend = axes$y + c(0, axes$y_len), 
             arrow = arrow(type = "closed", length = unit(10, 'pt')), 
             linewidth = 1
    ) + 
    annotate(geom = "text", x = axes$x + axes$x_len / 2, y = axes$y - 1.5, label = "UMAP1", size = 3) + 
    annotate(geom = "text", x = axes$x - 1.5, y = axes$y + axes$y_len / 2, label = "UMAP2", angle = "90", size = 3)
  
  # UMAP without labels and with legends for standalone UMAP to extract legend
  umap_plot_nolabel_legend <- ggplot(data = umap_df, mapping = aes(x = UMAP1, y = UMAP2, color = Label)) + 
    geom_point(alpha = 0.1, size = 0.0001) + 
    scale_color_manual(values = stallion_cols) + 
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      axis.ticks = element_blank(),  # Remove axis ticks
      axis.title = element_blank(),  # Increase axis title size
      axis.text = element_blank(),   # Increase axis text size
      legend.title = element_text(size = 14),  # Increase legend title size
      legend.text = element_text(size = 11, margin = margin(r = 20))    # Increase legend text size 
    ) +
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 3, alpha = 1))) +  # Increase dot size in legend
    annotate("segment", 
             x = axes$x, xend = axes$x + c(axes$x_len, 0), 
             y = axes$y, yend = axes$y + c(0, axes$y_len), 
             arrow = arrow(type = "closed", length = unit(10, 'pt')), 
             linewidth = 1
    ) + 
    annotate(geom = "text", x = axes$x + axes$x_len / 2, y = axes$y - 1.5, label = "UMAP1", size = 3) + 
    annotate(geom = "text", x = axes$x - 1.5, y = axes$y + axes$y_len / 2, label = "UMAP2", angle = "90", size = 3)
  legend_UMAP <- cowplot::get_legend(umap_plot_nolabel_legend)
  
  # UMAP without legend and without label for Fig 3a panel
  umap_plot_nolabel_nolegend <- ggplot(data = umap_df, mapping = aes(x = UMAP1, y = UMAP2, color = Label)) + 
    geom_point(alpha = 0.1, size = 0.0001) + 
    scale_color_manual(values = stallion_cols) + 
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      axis.ticks = element_blank(),        # Remove axis ticks
      axis.title = element_blank(),        # Increase axis title size
      axis.text = element_blank(),         # Increase axis text size
      legend.title = element_text(size = 14),  # Increase legend title size
      legend.text = element_text(size = 11)    # Increase legend text size
    ) + 
    NoLegend() +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +  # Increase dot size in legend
    annotate("segment", 
             x = axes$x, xend = axes$x + c(axes$x_len, 0), 
             y = axes$y, yend = axes$y + c(0, axes$y_len), 
             arrow = arrow(type = "closed", length = unit(10, 'pt')), 
             linewidth = 1
    ) + 
    annotate(geom = "text", x = axes$x + axes$x_len / 2, y = axes$y - 1.5, label = "UMAP1", size = 3) + 
    annotate(geom = "text", x = axes$x - 1.5, y = axes$y + axes$y_len / 2, label = "UMAP2", angle = 90, size = 3)
  
  FP <- vector(mode = "list", length = length(features))
  names(FP) <- features
  for (i in c(1:length(features))) {
    
    if (features[i] %in% rownames(obj)) {
      
      feature <- features[i]
      x <- obj@reductions[[reduction]]@cell.embeddings[, 1]
      y <- obj@reductions[[reduction]]@cell.embeddings[, 2]
      z <- FetchData(obj, assay = "RNA", layer = "data", vars = feature)[, feature]
      obj$z <- calc_density(x, y, z)
      FP[[features[i]]] <- FeaturePlot(obj, 
                                       features = c("z"),
                                       raster = FALSE,
                                       pt.size = 0.0001,
                                       label = FALSE,
                                       label.color = "black", 
                                       reduction = reduction, 
                                       dims = c(1, 2)
      ) + 
        ggtitle(label = features_title[i]) +
        theme(plot.title = element_text(vjust = - 3.5),
              plot.margin = unit(c(-5, -5, -5, -5), units = "mm")
        ) +
        NoLegend() +
        NoAxes()
      
    }
    obj$z <- NULL
    
  }
  
  design <- "
  1123456
  1177777
  1189abc
  ddddddd
  efghijk
  "
  
  empty_plot <- ggplot() + theme_void()
  
  plt <-
    (
      umap_plot_nolabel_nolegend +
        FP[[1]] + FP[[2]] + FP[[3]] + FP[[4]] + FP[[5]] + 
        plot_spacer() + 
        FP[[6]] + FP[[7]] + FP[[8]] + FP[[9]] + FP[[10]] + 
        plot_spacer() + 
        plot_spacer() + FP[[11]] + FP[[12]] + FP[[13]] + FP[[14]] + FP[[15]] + FP[[16]]
    ) +
    plot_layout(
      design  = design,
      guides  = "collect",
      heights = c(10, -1.5, 10, -1.5, 10)  # ensure these are valid for your layout
    ) &
    theme(plot.margin = margin(0, 0, 0, 0))
  
  
  gc()
  # rasterize(geom_point(), layers='Point', dpi=600)
  ggsave(plot = plt, filename = file.path(project_dir, "figures", "Fig2a.tiff"), device = "tiff", width = 350, height = 162, units = "mm", dpi = 250)
  ggsave(plot = umap_plot_label_legend, filename = file.path(project_dir, "figures", 'Annotated UMAP.pdf'), device = "pdf", width = 350, height = 162, units = "mm")
  ggsave(plot = umap_plot_label_legend, filename = file.path(project_dir, "figures", "Annotated UMAP.tiff"), device = "tiff", width = 350, height = 162, units = "mm", dpi = 250)
  ggsave(plot = legend_UMAP, filename = file.path(project_dir, "figures", 'UMAP Legend.pdf'), device = "pdf", width = 350, height = 162, units = "mm")
  ggsave(plot = legend_UMAP, filename = file.path(project_dir, "figures", "UMAP Legend.tiff"), device = "tiff", width = 350, height = 162, units = "mm", dpi = 250)
  
}


