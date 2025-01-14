figS3 <- function(obj, project_dir) {
  
  features <- c("HLA-A", "HLA-B", "HLA-C", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1")
  features_titles <- c("HLA-A", "HLA-B", "HLA-C", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1")
  
  FP <- plot_density_over_UMAP(obj, 
                               features = features, 
                               feature_titles = features_titles, 
                               umap_slot_name = "umap_2D")
  
  design <- "
  123
  444
  567
  888
  9ab
"
  plt <- FP[[1]] + FP[[2]] + FP[[3]] +  plot_spacer() + FP[[4]] + FP[[5]] + FP[[6]] + plot_spacer() + FP[[7]] + FP[[8]] + FP[[9]] + plot_layout(design = design, guides = 'collect') + plot_layout(design = design, guides = 'collect', heights = c(10, -1.5, 10, -1.5, 10)) & theme(plot.margin = unit(c(0, 0, 0, 0), units = "mm"))
  plt
  ggsave(file.path(project_dir, "figures", 'FigS3.tiff'), device = "tiff", width = 180, height = 180, units = "mm", dpi = 250)
  
}