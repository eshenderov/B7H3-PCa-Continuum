tableS3 <- function (obj, project_dir) {
  
  meta_data <- obj@meta.data
  summarized_meta_data <- meta_data %>% 
    group_by(Dataset, Subtype) %>% 
    summarize("# Patients" = n_distinct(Patient_ID), 
              "Cell Count" = n(), 
              "UMI (Mean)" = mean(nCount_RNA), 
              "Genes (Mean)" = mean(nFeature_RNA), 
              "Mito RNA (Mean)" = mean(percent.mito)
              )
  
  summarized_meta_data %>% 
    write.table(x = ., 
                file = file.path(project_dir, 
                                 "TableS3.csv"), 
                row.names = FALSE, 
                col.names = TRUE
                )
  
}
