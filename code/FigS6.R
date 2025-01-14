figS6 <- function(obj, project_dir) {
  
  features <- read_xlsx(path = file.path("/projects/eshenderov-hpc/Shivang/projects/scRNA-seq_atlas/Prostate/integrated_dataset", "Genesets_Emir.xlsx")) %>%
    as.list()
  features$`AR Score` <- features$`AR Score`[!(features$`AR Score` %in% NA)]
  features$`NE Score` <- features$`NE Score`[!(features$`NE Score` %in% NA)]
  
  obj@meta.data$AR_Score <- 
    FetchData(obj, 
              assay = "RNA", 
              layer = "data", 
              vars = c(features$`AR Score`) %>% unique()) %>% 
    as.matrix() %>% 
    rowSums2()
  obj@meta.data$NE_Score <- 
    FetchData(obj, 
              assay = "RNA", 
              layer = "data", 
              vars = c(features$`NE Score`) %>% unique()) %>% 
    as.matrix() %>% 
    rowSums2()
  
  
  cell_comp_score <- obj@meta.data[, c("CoarseAnnotation", "AR_Score", "NE_Score")]
  score_boxplot <- box_plot(df_expr = cell_comp_score, 
                            x_axis_name = "CoarseAnnotation", 
                            y_axis_name = "Score", 
                            group_col_name = "Geneset", 
                            col_fill = c("#B2182B", "#2166AC"), 
                            show_legend = TRUE
  ) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 16, color = "black"), 
          axis.text.y = element_text(size = 16, color = "black"), 
          axis.title.x = element_text(size = 16, face = "bold"), 
          axis.title.y = element_text(size = 16, face = "bold"), 
          legend.text = element_text(size = 16, color = "black"), 
          legend.title = element_text(size = 16, face = "bold")
    )
  score_boxplot
  ggsave(file.path(project_dir, "figures", 'FigS6.tiff'), device = "tiff", width = 360, height = 180, units = "mm", dpi = 250)
  
}