figS9 <- function (obj, project_dir) {

  obj_sub <- subset(obj, subset = (CoarseAnnotation %in% c("Luminal", "Basal")) & (Subtype %in% c("Healthy", "Benign", "BPH")))
  df_CD276 <- data.frame(Patient_ID = obj_sub$Patient_ID, 
                         Cell_Type = factor(obj_sub$CoarseAnnotation %>% as.character(), levels = c("Luminal", "Basal"), labels = c("Luminal", "Basal")), 
                         Subtype = factor(obj_sub$Subtype, levels = c("Healthy", "Benign")), 
                         Expr = FetchData(obj_sub, layer = "data", vars = "CD276")
  )
  df_CD276 <- df_CD276 %>% 
    group_by(Patient_ID, Cell_Type) %>% 
    summarize(CD276 = mean(CD276)) %>% 
    as.data.frame()
  plt <- box_plot(df_expr = df_CD276[, c(-1)], 
                  x_axis_name = "Cell_Type", 
                  y_axis_name = "CD276", 
                  col_fill = c("red", "blue"), 
                  group_plot = FALSE, 
                  show_legend = TRUE)
  plt <- plt + 
    labs(x = "Cell types", y = "Normalized CD276 Expression") + 
    theme(axis.title.x = element_text(size = 16, face = "bold"), 
          axis.title.y = element_text(size = 16, face = "bold"), 
          axis.text.x = element_text(size = 16), 
          axis.text.y = element_text(size = 16), 
          legend.text = element_text(size = 16), 
          legend.title = element_text(size = 16, face = "bold")
    )
  plt
  ggsave(file.path(project_dir, "figures", "FigS9.pdf"),
         device = "pdf", 
         width = 200, 
         height = 300, 
         units = "mm")
  
  ggsave(file.path(project_dir, "figures", 'FigS9.pdf'), 
         device = "tiff", 
         width = 200, 
         height = 300, 
         units = "mm")
  
  pdf(file = file.path(project_dir, "figures", 'FigS9.pdf'), 
      width  = convertUnit(unit(200, "mm"), "inches", valueOnly = TRUE),
      height = convertUnit(unit(300, "mm"), "inches", valueOnly = TRUE)
      )
  draw(plt)
  dev.off()
  dev.off()
  
  # Welch Two Sample t-test ------------------------------------------------------
  x <- df_CD276$CD276[df_CD276$Cell_Type %in% "Luminal"]
  y <- df_CD276$CD276[df_CD276$Cell_Type %in% "Basal"]
  t.test(x = x, y = y)
  
  # Result -----------------------------------------------------------------------
  # data:  x and y
  # t = 2.6032, df = 62.747, p-value = 0.01151
  # alternative hypothesis: true difference in means is not equal to 0
  # 95 percent confidence interval:
  #   0.008618151 0.065580536
  # sample estimates:
  #   mean of x  mean of y 
  # 0.10179828 0.06469894 
  
}