figS4 <- function(obj, project_dir, genes, gene_names) {
  
  features <- genes
  features_titles <- gene_names
  
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
ddd
eff
"
  
  plt <- (
    FP[[1]] + FP[[2]] + FP[[3]] + plot_spacer() +
      FP[[4]] + FP[[5]] + FP[[6]] + plot_spacer() +
      FP[[7]] + FP[[8]] + FP[[9]] + plot_spacer() +
      FP[[10]] + plot_spacer()
  ) +
    plot_layout(
      design  = design,
      guides  = "collect",
      # example heights: one per row (adjust as you like)
      heights = c(10, 1, 10, 1, 10, 1, 10)
    ) &
    theme(plot.margin = unit(c(0, 0, 0, 0), units = "mm"))
  
  plt
  
  ggsave(file.path(project_dir, "figures", 'FigS4.tiff'), device = "tiff", width = 180, height = 240, units = "mm", dpi = 600)
  #ggsave(file.path(project_dir, "figures", 'FigS4.pdf'), device = "pdf", width = 180, height = 120, units = "mm", dpi = 600)
  
  # # Fig Supp 3 -------------------------------------------------------------------
  # # Pathological composition of celltypes
  # cell_comp_path <- obj@meta.data %>% 
  #   group_by(CellAnnFine, Pathology) %>% 
  #   summarize(Counts = n())
  # ggplot(cell_comp_path, aes(x=CellAnnFine, fill=Pathology, y = Counts)) + 
  #   geom_bar(position = "fill", stat = "identity") + 
  #   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 28, color = "black"), 
  #         axis.text.y = element_text(size = 28, color = "black"), 
  #         axis.title.x = element_text(size = 28, face = "bold"), 
  #         axis.title.y = element_text(size = 28, face = "bold"), 
  #         legend.text = element_text(size = 28, color = "black"), 
  #         legend.title = element_text(size = 28, face = "bold")
  #   ) + 
  #   scale_fill_manual(values = brewer.pal(n = 7, name = "Set2"))
  # ggsave(file.path(project_dir, "figures", 'FigS4.tiff'), device = "tiff", width = 360, height = 180, units = "mm", dpi = 1800)
  # 
}
