figS4 <- functin(obj, project_dir) {
  
  features <- c("AR", "KLK2", "KLK3", "TMPRSS2", "FKBP5", "IGF1R")
  features_titles <- c("AR", "KLK2", "KLK3", "TMPRSS2", "FKBP5", "IGF1R")
  
  FP <- plot_density_over_UMAP(obj, 
                               features = features, 
                               feature_titles = features_titles, 
                               umap_slot_name = "umap_2D")
  
  design <- "
  123
  444
  567
"
  plt <- FP[[1]] + FP[[2]] + FP[[3]] + plot_spacer() + FP[[4]] + FP[[5]] + FP[[6]] + plot_layout(design = design, guides = 'collect') + plot_layout(design = design, guides = 'collect', heights = c(10, -1.5, 10)) & theme(plot.margin = unit(c(0, 0, 0, 0), units = "mm"))
  plt
  ggsave(file.path(project_dir, "figures", 'FigS4.tiff'), device = "tiff", width = 180, height = 120, units = "mm", dpi = 250)
  
  # Fig Supp 3 -------------------------------------------------------------------
  # Pathological composition of celltypes
  cell_comp_path <- obj@meta.data %>% 
    group_by(CellAnnFine, Pathology) %>% 
    summarize(Counts = n())
  ggplot(cell_comp_path, aes(x=CellAnnFine, fill=Pathology, y = Counts)) + 
    geom_bar(position = "fill", stat = "identity") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 28, color = "black"), 
          axis.text.y = element_text(size = 28, color = "black"), 
          axis.title.x = element_text(size = 28, face = "bold"), 
          axis.title.y = element_text(size = 28, face = "bold"), 
          legend.text = element_text(size = 28, color = "black"), 
          legend.title = element_text(size = 28, face = "bold")
    ) + 
    scale_fill_manual(values = brewer.pal(n = 7, name = "Set2"))
  ggsave(file.path(project_dir, "figures", 'FigS4.tiff'), device = "tiff", width = 360, height = 180, units = "mm", dpi = 250)
  
}