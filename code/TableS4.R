tableS4 <- function(project_dir) {
  
  markers <- prostate_celltype_markers(dataset = "whole_prostate")
  
  markers_lineage <- markers$Lineage %>% lapply(X = ., FUN = function(x) {
    
    paste(x, collapse = ", ")
    
  }) %>% as.data.frame()
  colnames(markers_lineage) <- names(markers$Lineage)
  rownames(markers_lineage) <- "Markers"
  markers_lineage <- df_transpose(markers_lineage)

  markers_coarse_annotation <- markers$Coarse%>% lapply(X = ., FUN = function(x) {
    
    paste(x, collapse = ", ")
    
  }) %>% as.data.frame()
  colnames(markers_coarse_annotation) <- names(markers$Coarse)
  rownames(markers_coarse_annotation) <- "Markers"
  markers_coarse_annotation <- df_transpose(markers_coarse_annotation)
  rbind(markers_lineage, markers_coarse_annotation) %>% 
    write.table(x = , 
                file = file.path(project_dir, 
                                 "tables", 
                                 "TableS4.xlsx"), 
                row.names = TRUE, 
                col.names = TRUE
                )
  
}