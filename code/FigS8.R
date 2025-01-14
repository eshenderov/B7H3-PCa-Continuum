fig2b <- function (obj, project_dir) {
  
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
  
  # Create gene x patient matrix of intrasample cov ----------------------------
  for (patient in patient_ID) {
    
    print(patient)
    
    obj_lum_sub <- subset(obj_lum, Patient_ID == patient)
    obj_lum_sub <- NormalizeData(obj_lum_sub, normalization.method = "RC")
    data <- FetchData(obj_lum_sub, vars = genes, layer = "data", assay = "RNA") %>% as.matrix()
    # data_mat <- GetAssayData(obj_lum_sub, layer = "counts", assay = "RNA")
    # lib_size <- colSums(data_mat)
    # cts <- FetchData(obj_lum_sub, vars = "CD276", layer = "counts") %>% as.matrix() %>% as.vector()
    # data <- (cts + 1) / lib_size
    # H_patient <- vector(mode = "numeric", length = 10000)
    # for (i in c(1:10000)) {
    #   
    #   data_mc <- sample(x = data, size = 0.2 * length(data), replace = FALSE)
    #   H_patient[i] <- ((- log(data_mc)) * data_mc) %>% sum()
    #   
    # }
    # H[patient] <- median(H_patient) / length(data) * 10000
    
    for (gene in genes) {
      
      
      df[gene, patient] <- sd(data[, gene] %>% as.vector()) / mean(data[, gene] %>% as.vector()) * 100
      
    }
    
    data <- GetAssayData(obj_lum_sub, layer = "counts", assay = "RNA")
    med_reads_per_cell[patient] <- colSums(data) %>% as.matrix() %>% as.vector() %>% median()
    cell_count[patient] <- dim(obj_lum_sub)[2]
    
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
  RLE_plot(vsd_df_corrected)
  
  # Calculate AR and NE Score ----------------------------------------------------
  vsd <- vsd_df_corrected
  features <- read_xlsx(path = file.path("/projects/eshenderov-hpc/Shivang/projects/scRNA-seq_atlas/Prostate/integrated_dataset", "Genesets_Emir.xlsx")) %>%
    as.list()
  features$`AR Score` <- features$`AR Score`[!(features$`AR Score` %in% NA)]
  features$`NE Score` <- features$`NE Score`[!(features$`NE Score` %in% NA)]
  AR_Score <- (vsd[features$`AR Score`, ] %>% as.matrix() %>% colSums2())/(features$`AR Score` %>% length())
  NE_Score <- (vsd[features$`NE Score`, ] %>% as.matrix() %>% colSums2())/(features$`NE Score` %>% length())
  
  # Subset expression data matrix for genes of interest --------------------------
  luminal_basal_genes <- c("STEAP1", "STEAP2", "PSCA", "FOLH1", "TACSTD2", "PVRL1", "DLL3", "CD276", "CD274", "PDCD1LG2", "PDCD1", "LAG3", "TIGIT", "TNFRSF4", "TNFRSF9", "CTLA4", "AR", "KLK2", "KLK3", "TMPRSS2", "FKBP5", "IGF1R")
  vsd_sub <- vsd[luminal_basal_genes, ]
  
  # Prepare data for heatmap -----------------------------------------------------
  anno <- meta_data %>% select(-c("Dataset"))
  vsd_sub <- vsd_sub[, anno$`Sample Name`]
  rownames(anno) <- anno$`Sample Name`
  AR_Score <- AR_Score[anno$`Sample Name`]
  NE_Score <- NE_Score[anno$`Sample Name`]
  
  # Assign colors to subtypes of PCa, AR and NE Score in heatmap -----------------
  colors <- c("#5E4FA2", "#5E4FA2", "#3288BD", "#ABDDA4", "#E6F598", "#E6F598", "#E6F598", "#E6F598", "#FEE08B", "#FEE08B", "#F46D43", "#D53E4F", "#9E0142")
  names(colors) <- c("Healthy", "Healthy", "Benign", "Tumor", "Tumor", "Tumor", "Tumor", "ICC/IDC", "ICC/IDC", "CRPC", "mixed-NEPC", "NEPC", "Small Cell")
  colors <- c("#5E4FA2", "#3288BD", "#ABDDA4", "#E6F598", "#FEE08B", "#F46D43", "#9E0142")
  names(colors) <- c("Healthy", "Benign", "Adeno", "CRPC", "mixed-NEPC", "NEPC", "Small Cell")
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
  
  # Create heatmap column section ------------------------------------------------
  col_ha <- HeatmapAnnotation(Subtype = anno$Subtype,
                              "AR_Score" = AR_Score %>% as.numeric() %>% unname(),
                              "NE_Score" = NE_Score %>% as.numeric() %>% unname(),
                              col = anno_colors_col,
                              annotation_label = c("Subtype", "AR Score", "NE Score"),
                              annotation_name_side = "right",
                              border = c(Subtype = FALSE, AR_Score = FALSE, NE_Score = FALSE),
                              simple_anno_size = unit(5, "mm")
  )
  
  # Assign colors to rows of genes based on their category -----------------------
  anno_colors_row <- list(Gene = c("#FF0000", "#49FF00", "#0092FF"))
  names(anno_colors_row$Gene) <- c("ADC", "Immune", "AR reg.")
  
  # Create heatmap row section ---------------------------------------------------
  row_ha <- HeatmapAnnotation(Gene = rep(x = c("ADC", "Immune", "AR reg."), times = c(7, 9, 6)) %>% 
                                factor(levels = c("ADC", "Immune", "AR reg.")),
                              col = anno_colors_row,
                              which = "row",
                              show_annotation_name = FALSE,
                              simple_anno_size = unit(2, "mm")
                              
  )
  
  # Create heatmap ---------------------------------------------------------------
  ht_plot <- Heatmap(matrix = vsd_sub %>% as.matrix(),
                     name = "log(BCC+1)",
                     col = colorpanel(n = 500,
                                      low = "#440154FF",
                                      mid = "#21908CFF",
                                      high = "#FDE725FF"),
                     row_split = rep(x = c("ADC", "Immune", "AR reg."), times = c(7, 9, 6)) %>% factor(levels = c("ADC", "Immune", "AR reg.")),
                     cluster_columns = FALSE,
                     cluster_rows = FALSE,
                     row_gap = unit(x = c(1, 2), units = "mm"),
                     border = TRUE,
                     left_annotation = row_ha,
                     top_annotation = col_ha,
                     row_labels = c("STEAP1", "STEAP2", "PSCA", "PSMA", "TROP-2", "NECTIN1", "DLL3", "B7-H3", "PD-1", "PD-L2", "PD-L1", "LAG3", "TIGIT", "OX40", "4-1BB", "CTLA4", "AR", "KLK2", "KLK3", "TMPRSS2", "FKBP5", "IGF1R"),
                     show_column_names = FALSE,
                     row_title = NULL,
                     row_names_side = "left",
                     width = (ncol(vsd_sub) * unit(2, "mm")) * 1.2,
                     height = (nrow(vsd_sub) + 3) * unit(3.5, "mm") * 1.2,
                     row_names_gp = gpar(fontsize = 10)
  )
  
  # Save heatmap into pdf --------------------------------------------------------
  pdf(file = file.path(project_dir, "figures", 'Fig2b.pdf'), 
      width = convertUnit(x = (ncol(vsd_sub) * unit(2, "mm")) * 1.7, 
                          unitTo = "inches"), 
      height = convertUnit(x = (nrow(vsd_sub) + 3) * unit(3.5, "mm") * 1.7, 
                           unitTo = "inches"))
  ComplexHeatmap::draw(ht_plot)
  dev.off()
  dev.off()
  
}




library(RColorBrewer)

obj <- readRDS(file = "/media/singhlab/B684-19A61/Single cell datasets/Prostate/integrated_dataset/Integrated_clustered_annotated_data_trim.RData")
obj_lum <- subset(obj, idents = c("Luminal", "Club", "Hillock", "Neuroendocrine"))

meta_data <- obj_lum@meta.data[, c("Subtype", "nCount_RNA","Patient_ID", "Sample_ID")]
meta_data <- meta_data[!(meta_data$Patient_ID %in% c("GSE172357_BPH")), ]
# meta_data$Patient_ID <- NULL
order <- c(which(meta_data$Subtype %in% c("Healthy")), 
           which(meta_data$Subtype %in% c("Benign")), 
           which(meta_data$Subtype %in% c("Adeno")), 
           which(meta_data$Subtype %in% c("ICC/IDC")), 
           which(meta_data$Subtype %in% "CRPC"), 
           which(meta_data$Subtype %in% "NEPC"), 
           which(meta_data$Subtype %in% "Small Cell")
           ) %>% rev()
meta_data <- meta_data[order, ]
genes <- c("STEAP1", "STEAP2", "PSCA", "FOLH1", "TACSTD2", "PVRL1", "DLL3", "CD276", "CD274", "PDCD1LG2", "PDCD1", "LAG3", "TIGIT", "TNFRSF4", "TNFRSF9", "CTLA4", "AR", "KLK2", "KLK3", "TMPRSS2", "FKBP5", "IGF1R")

sample_ID <- meta_data$Sample_ID %>% unique() %>% as.matrix() %>% as.vector()
H <- vector(mode = "numeric", length = length(sample_ID))
names(H) <- sample_ID

df <- data.frame(matrix(nrow = length(genes), ncol = length(sample_ID)), row.names = genes)
colnames(df) <- sample_ID
pathology <- meta_data[, c("Sample_ID", "Subtype")] %>% unique() %>% remove_rownames()
rownames(pathology) <- pathology$Sample_ID
pathology <- pathology[rownames(pathology) %in% sample_ID, ]
cell_count <- vector(mode = "integer", length = length(sample_ID))
names(cell_count) <- sample_ID
med_reads_per_cell <- vector(mode = "numeric", length = length(sample_ID))
names(med_reads_per_cell) <- sample_ID

for (sample in sample_ID) {
  
  print(sample)
  
  obj_lum_sub <- subset(obj_lum, Sample_ID == sample)
  obj_lum_sub <- NormalizeData(obj_lum_sub, normalization.method = "RC")
  data <- FetchData(obj_lum_sub, vars = genes, layer = "data", assay = "RNA") %>% as.matrix()
  # data_mat <- GetAssayData(obj_lum_sub, layer = "counts", assay = "RNA")
  # lib_size <- colSums(data_mat)
  # cts <- FetchData(obj_lum_sub, vars = "CD276", layer = "counts") %>% as.matrix() %>% as.vector()
  # data <- (cts + 1) / lib_size
  # H_patient <- vector(mode = "numeric", length = 10000)
  # for (i in c(1:10000)) {
  #   
  #   data_mc <- sample(x = data, size = 0.2 * length(data), replace = FALSE)
  #   H_patient[i] <- ((- log(data_mc)) * data_mc) %>% sum()
  #   
  # }
  # H[patient] <- median(H_patient) / length(data) * 10000
  
  for (gene in genes) {
    
    
    df[gene, sample] <- sd(data[, gene] %>% as.vector()) / mean(data[, gene] %>% as.vector()) * 100
    
  }
  
  data <- GetAssayData(obj_lum_sub, layer = "counts", assay = "RNA")
  med_reads_per_cell[sample] <- colSums(data) %>% as.matrix() %>% as.vector() %>% median()
  cell_count[sample] <- dim(obj_lum_sub)[2]
  
}
pathology$Subtype <- factor(pathology$Subtype, levels = c("Healthy", "Benign", "Adeno", "ICC/IDC", "CRPC", "NEPC", "Small Cell"))
anno_colors_col <- list(Subtype = c("#5E4FA2", "#3288BD", "#ABDDA4", "#E6F598", "#FEE08B", "#F46D43", "#9E0142"),
                        "Cell counts" = colorRamp2(breaks = seq(min(cell_count), max(cell_count), length.out = 100),
                                                   hcl_palette = "SunsetDark",
                                                   transparency = 0,
                                                   reverse = TRUE),
                        "Reads per cell" = colorRamp2(breaks = seq(min(med_reads_per_cell), max(med_reads_per_cell), length.out = 100),
                                                      hcl_palette = "Purples 3",
                                                      transparency = 0,
                                                      reverse = TRUE)
                        )
names(anno_colors_col$Subtype) <- levels(pathology$Subtype)
col_ha <- HeatmapAnnotation(Subtype = pathology$Subtype,
                            "Cell counts" = cell_count %>% as.numeric() %>% unname(),
                            "Reads per cell" = med_reads_per_cell %>% as.numeric() %>% unname(),
                            col = anno_colors_col,
                            annotation_label = c("Subtype", "Cell counts", "Reads per cell"),
                            annotation_name_side = "right",
                            border = c(Pathology = FALSE, AR_Score = TRUE, NE_Score = TRUE),
                            simple_anno_size = unit(5, "mm")
)
anno_colors_row <- list(Gene = c("#FF0000", "#49FF00", "#0092FF"))
names(anno_colors_row$Gene) <- c("ADC", "Immune", "AR reg.")
# poly_ht <- data.frame(x = c(1, 1, 2, 1, 2, 2),
#                       y = c(1, 2, 2, 1, 2, 1),
#                       id = factor(c(1, 1, 1, 2, 2, 2), levels = c(1, 2)),
#                       value = factor(c(1, 1, 1, 2, 2, 2), levels = c(1, 2)))
row_ha <- HeatmapAnnotation(Gene = rep(x = c("ADC", "Immune", "AR reg."), times = c(7, 9, 6)) %>%
                              factor(levels = c("ADC", "Immune", "AR reg.")),
                            col = anno_colors_row,
                            which = "row",
                            show_annotation_name = FALSE,
                            simple_anno_size = unit(2, "mm")
)
ht_plot <- Heatmap(matrix = df %>% as.matrix(),
                   name = "% CV",
                   col = colorpanel(n = 500,
                                    low = "#2A6A99",
                                    mid = "white",
                                    high = "#B26E39"),
                   row_split = rep(x = c("ADC", "Immune", "AR reg."), times = c(7, 9, 6)) %>% factor(levels = c("ADC", "Immune", "AR reg.")),
                   cluster_columns = FALSE,
                   cluster_rows = FALSE,
                   row_gap = unit(x = c(1, 2), units = "mm"),
                   border = TRUE,
                   left_annotation = row_ha,
                   top_annotation = col_ha,
                   row_labels = c("STEAP1", "STEAP2", "PSCA", "PSMA", "TROP-2", "NECTIN1", "DLL3", "B7-H3", "PD-1", "PD-L2", "PD-L1", "LAG3", "TIGIT", "OX40", "4-1BB", "CTLA4", "AR", "KLK2", "KLK3", "TMPRSS2", "FKBP5", "IGF1R"),
                   show_column_names = FALSE,
                   row_title = NULL,
                   row_names_side = "left",
                   width = (ncol(df) * unit(2, "mm")) * 1.2,
                   height = (nrow(df) + 3) * unit(3.5, "mm") * 1.2,
                   row_names_gp = gpar(fontsize = 10)
)
pdf(file = file.path(project_dir, "figures", 'FigS8.pdf'), width = convertUnit(x = (ncol(df) * unit(2, "mm")) * 1.7, unitTo = "inches"), height = convertUnit(x = (nrow(df) + 3) * unit(3.5, "mm") * 1.7, unitTo = "inches"))
ComplexHeatmap::draw(ht_plot)
dev.off()
dev.off()
