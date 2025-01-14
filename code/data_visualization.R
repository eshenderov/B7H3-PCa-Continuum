plot_annotated_umap <- function(obj, reduction, cols, label = FALSE, legend = TRUE, cell_annotation_column) {

  # obj is a R Seurat object
  # cols is a named vector of colors, with cell types as names
  umap_df <- cbind(obj@reductions[[reduction]]@cell.embeddings,
                   data.frame(Label = obj@meta.data[[cell_annotation_column]])
  )
  colnames(umap_df) <- c("UMAP1", "UMAP2", "Label")

  centroids <- umap_df %>%
    group_by(Label) %>%
    summarize(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2), .groups = 'drop')

  axes <- list(x = -12, y = -15, x_len = 5, y_len = 5)

  umap_plot <- ggplot(data = umap_df, mapping = aes(x = UMAP1, y = UMAP2, color = Label)) +
    geom_point(alpha = 0.1, size = 0.0001) +
    scale_color_manual(values = cols) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      axis.ticks = element_blank(),  # Remove axis ticks
      axis.title = element_blank(),  # Increase axis title size
      axis.text = element_blank(),   # Increase axis text size
      legend.title = element_text(size = 14),  # Increase legend title size
      legend.text = element_text(size = 11)    # Increase legend text size
      )

  if (legend == FALSE) {

    umap_plot <- umap_plot + NoLegend()

  }

  if (label == TRUE) {

    umap_plot <- umap_plot +
      geom_text_repel(data = centroids,
                      aes(label = Label),
                      size = 3,
                      nudge_y = 0.1,
                      color = "black"
                      )

  }

  umap_plot <- umap_plot +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +  # Increase dot size in legend
    annotate("segment",
             x = axes$x, xend = axes$x + c(axes$x_len, 0),
             y = axes$y, yend = axes$y + c(0, axes$y_len),
             arrow = arrow(type = "closed", length = unit(10, 'pt')),
             linewidth = 1
             ) +
    annotate(geom = "text", x = axes$x + axes$x_len / 2, y = axes$y - 1.5, label = "UMAP1", size = 3) +
    annotate(geom = "text", x = axes$x - 1.5, y = axes$y + axes$y_len / 2, label = "UMAP2", angle = "90", size = 3)

}

plot_density_over_UMAP <- function(obj, features, feature_titles, umap_slot_name) {

  features <- features[features %in% c(rownames(obj), colnames(obj@meta.data))]
  not_features <- features[!(features %in% c(rownames(obj), colnames(obj@meta.data)))]
  if (!(length(not_features) == 0)) {

    paste0("These features were not found: ", paste(not_features, sep = ", ")) %>% message()

  }

  plt <- vector(mode = "list", length = length(features))
  names(plt) <- features
  df <- FetchData(obj, assay = "RNA", layer = "data", vars = features) %>% as.data.frame()
  for (i in c(1:length(features))) {

    feature <- features[i]
    feature_title <- feature_titles[i]
    paste0("Plotting feature: ", feature) %>% message()
    x <- obj@reductions[[umap_slot_name]]@cell.embeddings[, 1]
    y <- obj@reductions[[umap_slot_name]]@cell.embeddings[, 2]
    data <- df[[feature]]
    if (is.factor(data)) {

      data <- as.integer(data)

    }

    z <- calc_density(x, y, data)
    names(z) <- colnames(obj)
    obj$z <- z
    plt[[feature]] <- FeaturePlot(obj,
                                  features = c("z"),
                                  reduction = umap_slot_name,
                                  raster = FALSE,
                                  pt.size = 0.0005,
                                  label = FALSE,
                                  label.color = "black"
                                    ) +
      ggtitle(label = feature) +
      theme(plot.title = element_text(vjust = - 3.5),
            plot.margin = unit(c(-5, -5, -5, -5), units = "mm")
            ) +
      NoLegend() +
      NoAxes()

  }
  return(plt)

}

plot_prostate_features <- function(obj,
                                    celltypes = "none",
                                    variable_column = NA,
                                    column_type = c("identity", "score"),
                                    reduction = "umap"
                                    ) {

  celltype_markers <- list("Luminal" = c("KLK2", "KLK3", "ACPP", "MSMB"),
                           "Hillock" = c("SCGB3A1", "PIGR", "RARRES1"),
                           "Club" = c("KRT13", "GDF15", "GPX2"),
                           "Basal" = c("KRT5", "KRT15", "DST"),
                           "Neuroendocrine" = c("CALCA", "CHGA", "GRP", "SCG2", "TPH1", "PCSK1N"),
                           "Endothelial" = c("VWF", "FLT1", "SELE", "ENG", "PECAM1"),
                           "Fibroblast" = c("COL1A1", "COL1A2", "FBLN1", "DCN"),
                           "SMC" = c("MYH11", "RGS5", "ACTA2", "CALD1", "TAGLN", "MYL9"),
                           "Myeloid" = c("PTPRC", "CD68", "CD14", "C1QA", "C1QB", "LYZ", "CD163"),
                           "TCell" = c("PTPRC", "CD3D", "CD3E", "CD3G", "IL7R"),
                           "BCell" = c("PTPRC", "CD19", "CD79A", "LY9", "CXCR5", "BANK1"),
                           "NK" = c("PTPRC", "NCAM1", "NKG7", "CD3E"),
                           "Plasma" = c("PTPRC", "IGKC", "IGHA1", "IGHA2", "IGHG3", "IGHG4"),
                           "Mast" = c("TPSAB1", "MS4A2")
                           )
  features <- NA

  if (celltypes == "all") {

    celltypes <- names(celltype_markers)
    features <- celltype_markers[celltypes]

  } else if (celltypes == "none") {

    if (((variable_column %in% colnames(obj@meta.data)) %>% sum()) > 0) {

      if ((column_type == "identity") & (length(variable_column) == 1)) {

        celltypes <- obj@meta.data[, variable_column] %>% unique()
        features <- celltype_markers[celltypes]

      } else if (column_type == "score") {

        features <- list(variable_column)
        names(features) <- "Scores"

      } else {

        writeLines("Provide only one column for celltypes/identities.")

      }

    } else {

      print("No celltype selected and celltype column not found.")
      return(NA)

    }

  } else {

    features <- celltype_markers[celltypes]

  }

  patched_plot <- vector(mode = "list", length = length(features))
  names(patched_plot) <- names(features)
  for (feature in names(features)) {

    plt <- plot_density_over_UMAP(obj,
                                  features = features[[feature]],
                                  reduction
                                  )

    umap_plot <-  DimPlot(obj,
                          reduction = reduction,
                          label = TRUE,
                          repel = TRUE,
                          raster = FALSE,
                          pt.size = 0.0005
    )
    design <- "
      1112
      1112
      3456
      789a
      "
    patched_plot[[feature]] <- umap_plot + guide_area()
    for (plot in plt) {

      patched_plot[[feature]] <- patched_plot[[feature]] + plot

    }
    patched_plot[[feature]] + plot_layout(design = design,
                                          guides = 'collect',
                                          heights = c(10, 10, 10, 10)
    )

  }

  return(patched_plot)

}

explore_unknown_subtypes <- function(obj, known_celltypes, known_celltypes_column = NA) {

  if (!is.na(known_celltypes_column)) {

    obj_sub <- subset(obj, subset = !!sym(known_celltypes_column) == known_celltypes)

  } else {

    obj_sub <- subset(obj, idents = known_celltypes)

  }
  obj_sub <- obj_sub %>% NormalizeData() %>% ScaleData(features = rownames(obj_sub))
  obj_sub.markers <- FindAllMarkers(obj_sub, only.pos = TRUE, logfc.threshold = 1.4)
  obj_sub.markers <- obj_sub.markers[obj_sub.markers$p_val_adj < 0.05, ]
  obj_sub.markers <- obj_sub.markers[!duplicated(obj_sub.markers$gene), ]
  geneset <- dcast(data = obj_sub.markers[ ,c("cluster", "gene", "avg_log2FC")],
                   gene ~ cluster,
                   value.var = "avg_log2FC") %>% as.data.frame()
  rownames(geneset) <- geneset$gene
  geneset$gene <- NULL
  geneset[is.na(geneset)] <- 0
  obj_data <- FetchData(object = obj, layer = "data", vars = rownames(geneset))
  obj_data <- obj_data[, rownames(geneset)]
  score_df <- data.frame(as.matrix(obj_data)%*%as.matrix(geneset))
  colnames(score_df) <- paste0(colnames(score_df), "_Score")
  obj <- AddMetaData(obj, metadata = score_df)
  plt_score <- plot_prostate_features(obj,
                                     variable_column = colnames(score_df),
                                     column_type = "score",
                                     reduction = "umap"
                                     )

  return(list(plt_score, geneset))

}
