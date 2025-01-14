library(RColorBrewer)

obj <- readRDS(file = "/media/singhlab/B684-19A61/Single cell datasets/Prostate/integrated_dataset/Integrated_clustered_annotated_data_trim.RData")
obj_lum <- subset(obj, idents = c("Luminal", "MSMB+", "UK3"))

meta_data <- obj_lum@meta.data[, c("Pathology", "nCount_RNA", "Patient_ID")]
meta_data <- meta_data[!(meta_data$Patient_ID %in% c("GSE172357_BPH_P4", "GSE141445_T_P2", "GSE172357_BPH_P6")), ]
# meta_data$Patient_ID <- NULL
order <- c(which(meta_data$Pathology %in% c("NEPC")),
           which(meta_data$Pathology %in% c("CRPC")),
           which(meta_data$Pathology %in% c("ICC/IDC")),
           which(meta_data$Pathology %in% c("Tumor")),
           which(meta_data$Pathology %in% "BPH"),
           which(meta_data$Pathology %in% "Benign"),
           which(meta_data$Pathology %in% "Normal")
)
meta_data <- meta_data[order, ]
genes <- c("STEAP1", "STEAP2", "PSCA", "FOLH1", "TACSTD2", "PVRL1", "DLL3", "CD276", "CD274", "PDCD1LG2", "PDCD1", "LAG3", "TIGIT", "TNFRSF4", "TNFRSF9", "CTLA4", "AR", "KLK2", "KLK3", "TMPRSS2", "FKBP5", "IGF1R")

patient_ID <- meta_data$Patient_ID %>% unique() %>% as.matrix() %>% as.vector()
H <- vector(mode = "numeric", length = length(patient_ID))
names(H) <- patient_ID

df <- data.frame(matrix(nrow = length(genes), ncol = length(patient_ID)), row.names = genes)
colnames(df) <- patient_ID
pathology <- meta_data[, c("Patient_ID", "Pathology")] %>% unique() %>% remove_rownames()
rownames(pathology) <- pathology$Patient_ID
pathology <- pathology[rownames(pathology) %in% patient_ID, ]
cell_count <- vector(mode = "integer", length = length(patient_ID))
names(cell_count) <- patient_ID
med_reads_per_cell <- vector(mode = "numeric", length = length(patient_ID))
names(med_reads_per_cell) <- patient_ID

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

anno_colors_col <- list(Pathology = brewer.pal(n = length(pathology$Pathology %>% levels()), name = "Set2"),
                        "Cell counts" = colorRamp2(breaks = seq(min(cell_count), max(cell_count), length.out = 100),
                                                hcl_palette = "SunsetDark",
                                                transparency = 0,
                                                reverse = TRUE),
                        "Reads per cell" = colorRamp2(breaks = seq(min(med_reads_per_cell), max(med_reads_per_cell), length.out = 100),
                                                hcl_palette = "Purples 3",
                                                transparency = 0,
                                                reverse = TRUE)
                        )
names(anno_colors_col$Pathology) <- levels(pathology$Pathology)
col_ha <- HeatmapAnnotation(Pathology = pathology$Pathology,
                            "Cell counts" = cell_count %>% as.numeric() %>% unname(),
                            "Reads per cell" = med_reads_per_cell %>% as.numeric() %>% unname(),
                            col = anno_colors_col,
                            annotation_label = c("Pathology", "Cell counts", "Reads per cell"),
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
                   width = (ncol(df) * unit(3.5, "mm")) * 1.2,
                   height = (nrow(df) + 3) * unit(3.5, "mm") * 1.2,
                   row_names_gp = gpar(fontsize = 10)
)
pdf(file = file.path(working_dir, 'FigS8.pdf'), width = convertUnit(x = (ncol(df) * unit(3.5, "mm")) * 1.5, unitTo = "inches"), height = convertUnit(x = (nrow(df) + 3) * unit(3.5, "mm") * 1.7, unitTo = "inches"))
draw(ht_plot)
dev.off()
dev.off()
