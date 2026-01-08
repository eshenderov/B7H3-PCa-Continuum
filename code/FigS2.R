figS2 <- function (obj, project_dir) {
  
  obj_epithelial <- subset(obj, 
                           subset = CoarseAnnotation %in% 
                             c("Luminal", "Club", "Hillock", "Neuroendocrine")
  )
  
  obj_epithelial[["log10(nCount_RNA)"]] <- log10(obj_epithelial$nCount_RNA)
  p1 <- VlnPlot(obj_epithelial, 
                features = c("log10(nCount_RNA)"), 
                pt.size = 0.001, 
                alpha = 0) + 
    theme(axis.title.x = element_blank()) + 
    NoLegend()
  p2 <- VlnPlot(obj_epithelial, 
                features = c("nFeature_RNA"), 
                pt.size = 0.001, 
                alpha = 0) + 
    theme(axis.title.x = element_blank()) + 
    NoLegend()
  p3 <- VlnPlot(obj_epithelial, 
                features = c("percent.mito"), 
                pt.size = 0.001, 
                alpha = 0) + 
    theme(axis.title.x = element_blank()) + 
    NoLegend()
  p4 <- VlnPlot(obj_epithelial, 
                features = c("percent.ribo"), 
                pt.size = 0.001, 
                alpha = 0) + 
    theme(axis.title.x = element_blank()) + 
    NoLegend()
  plt <- wrap_plots(p1, p2, p3, p4)
  ggsave(filename = "FigS2.pdf", path = file.path(project_dir, "figures"), plot = plt, device = "pdf", width = 150, height = 150, units = "mm", dpi = 600)
  ggsave(filename = "FigS2.tiff", path = file.path(project_dir, "figures"), plot = plt, device = "tiff", width = 150, height = 150, units = "mm", dpi = 600)
  
}
