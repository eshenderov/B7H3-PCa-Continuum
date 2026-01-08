figS7 <- function (obj, project_dir) {
  
  # Subset counts data for epithelial cells minus Basal cells --------------------
  obj_epithelial <- subset(obj, 
                           subset = CoarseAnnotation %in% 
                             c("Luminal", "Club", "Hillock", "Neuroendocrine")
  )
  obj_epithelial <- NormalizeData(obj_epithelial, normalization.method = "RC")
  data_epithelial <- LayerData(obj_epithelial, assay = "RNA", layer = "counts")
  n_cells_epithelial <- dim(obj_epithelial)[2]
  
  # Separate meta data of subsetted epithelial cells above -----------------------
  meta_data <- obj_epithelial@meta.data[, c("Subtype", "Dataset", "Patient_ID")]
  meta_data <- meta_data %>% unique() %>% remove_rownames()
  meta_data$"Sample Name" <- paste0(meta_data$Patient_ID, "_", meta_data$Subtype)
  
  # Create a matrix to hold pseudobulk data of subsetted epithelial cells for all patients
  data <- data.frame(matrix(nrow = dim(obj_epithelial)[1], ncol = nrow(meta_data)), row.names = rownames(obj_epithelial))
  colnames(data) <- meta_data$"Sample Name"
  
  # Create pseudobulk data -------------------------------------------------------
  obj_meta_data <- obj_epithelial@meta.data[, c("Subtype", "Dataset", "Patient_ID")]
  for (i in c(1:nrow(meta_data))) {
    
    subtype <- meta_data$Subtype[i]
    patient_id <- meta_data$Patient_ID[i]
    sample_name <- paste0(patient_id, "_", subtype)
    sample_name %>% green() %>% cat(., "\n")
    which_cells <- ((obj_meta_data$Patient_ID %in% patient_id) & 
                      (obj_meta_data$Subtype %in% subtype)) %>% 
      which()
    data[, sample_name] <- 
      data_epithelial[, which_cells] %>% 
      as.matrix() %>% 
      rowSums2()
    
  }
  
  # Normalize counts data to relative counts -------------------------------------
  data <- df_transpose(df_transpose(data) / colSums(data) * 1e9)
  
  # Remove spurious samples ------------------------------------------------------
  meta_data <- meta_data[!(meta_data$Dataset %in% c("GSE172357")), ]
  meta_data$Dataset <- factor(meta_data$Dataset, levels = unique(meta_data$Dataset))
  data <- data[, meta_data$`Sample Name`]
  
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
  data <- data[, meta_data$`Sample Name`]
  
  # Batch correction -------------------------------------------------------------
  vsd <- log1p(data[, meta_data$"Sample Name" %>% as.character()] %>% as.matrix())
  vsd_df_uncorrected <- as.data.frame(vsd)
  RLE_plot(vsd_df_uncorrected)
  vsd_df_corrected <- 
    limma::removeBatchEffect(vsd, 
                             batch = meta_data$Dataset, 
                             group = meta_data$Subtype) %>% 
    as.data.frame()
  
  p1 <- RLE_plot(vsd_df_uncorrected) +
    ggtitle(label = "Uncorrected") +
    theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))
  p2 <- RLE_plot(vsd_df_corrected) +
    ggtitle(label = "Corrected") +
    theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))
  p3 <- p1 / p2
  ggsave(filename = "FigS7.pdf", path = file.path(project_dir, "figures"), plot = p3, device = "pdf", width = 350, height = 300, units = "mm", dpi = 600)
  ggsave(filename = "FigS7.tiff", path = file.path(project_dir, "figures"), plot = p3, device = "tiff", width = 350, height = 300, units = "mm", dpi = 600)
  
}
