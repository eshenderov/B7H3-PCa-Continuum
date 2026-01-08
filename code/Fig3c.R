obj_epithelial <- NormalizeData(obj_epithelial, 
                                normalization.method = "LogNormalize",
                                scale.factor = 100000)
df_exp_epi <- data.frame(FetchData(object = obj_epithelial, 
                                   vars = c("NECTIN1", "NECTIN4", "TACSTD2", "CD276", "PSCA", "DLL3", "CoarseAnnotation", "Dataset", "Sample_ID", "Subtype"), 
                                   assay = "RNA", 
                                   layer = "counts")
)
colnames(df_exp_epi) <- c("NECTIN1", "NECTIN4", "TROP2", "B7-H3", "PSCA", "DLL3", "Celltype", "Dataset", "Sample", "Subtype")

df_exp_epi <- df_exp_epi %>% filter(Celltype %in% c("Luminal", "Neuroendocrine"))

dt <- as.data.table(df_exp_epi)

df_exp_epi_sum <- dt[, .(
  NECTIN1 = sum(NECTIN1, na.rm = TRUE),
  NECTIN4 = sum(NECTIN4, na.rm = TRUE), 
  TROP2 = sum(TROP2, na.rm = TRUE), 
  `B7-H3` = sum(`B7-H3`, na.rm = TRUE), 
  PSCA = sum(PSCA, na.rm = TRUE), 
  DLL3 = sum(DLL3, na.rm = TRUE)
), by = .(Dataset, Sample, Subtype)]

# Extract the row as a numeric vector,
# preserve names = sample names (column names of vsd_sub)
nectin4_high <- vsd_sub["NECTIN4", ]
nectin4_vec <- vsd_sub["NECTIN4", , drop = TRUE]
ord_idx <- order(nectin4_vec %>% as.numeric(), decreasing = TRUE)
nectin4_high <- nectin4_high[, ord_idx]
nectin4_high <- nectin4_high[, colnames(nectin4_high)[c(1:floor(0.92*dim(nectin4_high)[2]))]]
colnames(nectin4_high)[length(colnames(nectin4_high))]
df_exp_epi_sum <- as.data.frame(df_exp_epi_sum)
rownames(df_exp_epi_sum) <- df_exp_epi_sum$Sample
df_exp_epi_sum["HP95_T_unsorted", "NECTIN4"]

rownames(df_exp_epi_sum) <- df_exp_epi_sum$Sample
pie_df <- df_exp_epi_sum[, c("NECTIN1", "NECTIN4", "TROP2", "B7-H3", "PSCA", "DLL3")] %>% as.data.frame()
rownames(pie_df) <- df_exp_epi_sum$Sample
colnames(pie_df) <- c("NECTIN1", "NECTIN4", "TACSTD2", "CD276", "PSCA", "DLL3")

th <- 5

pie_df$CD276_NECTIN4_CoExp <- dplyr::case_when(
  pie_df$CD276 > th & pie_df$NECTIN4 > th ~ "Both",
  pie_df$CD276 > th & pie_df$NECTIN4 <= th ~ "CD276",
  pie_df$CD276 <= th & pie_df$NECTIN4 > th ~ "NECTIN4",
  TRUE ~ "None"
)

pie_df$CD276_PSCA_CoExp <- dplyr::case_when(
  pie_df$CD276 > th & pie_df$PSCA > th ~ "Both",
  pie_df$CD276 > th & pie_df$PSCA <= th ~ "CD276",
  pie_df$CD276 <= th & pie_df$PSCA > th ~ "PSCA",
  TRUE ~ "None"
)

# Set factor levels in desired order
pie_df$CD276_NECTIN4_CoExp <- factor(
  pie_df$CD276_NECTIN4_CoExp,
  levels = c("Both", "CD276", "NECTIN4", "None")
)

# Set factor levels in desired order
pie_df$CD276_PSCA_CoExp <- factor(
  pie_df$CD276_PSCA_CoExp,
  levels = c("Both", "CD276", "PSCA", "None")
)

pie_df <- as.data.frame(pie_df)

# Count proportions
plot_CD276_NECTIN4_df <- pie_df %>%
  dplyr::count(CD276_NECTIN4_CoExp) %>%
  mutate(percent = n / sum(n) * 100)

plot_CD276_PSCA_df <- pie_df %>%
  dplyr::count(CD276_PSCA_CoExp) %>%
  mutate(percent = n / sum(n) * 100)

# Pie chart

# Pie Chart colors
pie_colors <- c(
  "Both"    = "#9E6A9D",  # red
  "CD276"   = "#3F7FBF",  # blue
  "NECTIN4" = "#2FA37B",  # green
  "None"    = "#B5B5B5"   # gray
)

# assume plot_CD276_NECTIN1_df already completed to include all 4 categories
df <- plot_CD276_NECTIN4_df %>%
  mutate(CD276_NECTIN4_CoExp = factor(CD276_NECTIN4_CoExp, levels = names(pie_colors))) %>%
  complete(CD276_NECTIN4_CoExp = names(pie_colors), fill = list(percent = 0)) %>%
  arrange(CD276_NECTIN4_CoExp) %>%
  mutate(
    ymax = cumsum(percent),
    ymin = lag(ymax, default = 0),
    ymid = (ymin + ymax) / 2,
    lab  = ifelse(percent > 0, paste0(CD276_NECTIN4_CoExp, " (", round(percent, 1), "%)"), "")
  )

p_CD276_NECTIN4 <- ggplot(df, aes(x = 1, y = percent, fill = CD276_NECTIN4_CoExp)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y", clip = "off") +   # allow labels outside
  # Less horizontal span -> bigger pie; still some room for labels
  xlim(c(0.5, 1.5)) +
  scale_fill_manual(
    name   = "Co-expression",
    values = pie_colors,
    breaks = names(pie_colors),
    drop   = FALSE
  ) +
  theme_void(base_size = 14) +
  theme(
    legend.position  = "right",
    legend.title     = element_text(size = 16, face = "bold"),
    legend.text      = element_text(size = 14),
    legend.key.size  = unit(1.2, "cm"),
    plot.margin      = margin(10, 40, 10, 10)  # extra right margin for longer labels
  )
p_CD276_NECTIN4

ggsave(file.path(file.path(project_dir, "figures"), paste0("PDX_PSCA_NECTIN1/FigPDX_h", ".tiff")),
       plot = p_CD276_NECTIN1, device = "tiff",
       width = 120, height = 120, units = "mm", dpi = 600)
ggsave(file.path(file.path(project_dir, "figures"), paste0("PDX_PSCA_NECTIN1/FigPDX_h", ".pdf")),
       plot = p_CD276_NECTIN1, device = "pdf",
       width = 120, height = 120, units = "mm", dpi = 600)

pie_colors <- c(
  "Both"    = "#9E6A9D",  # red
  "CD276"   = "#3F7FBF",  # blue
  "PSCA"    = "#2FA37B",  # green
  "None"    = "#B5B5B5"   # gray
)

# Pie chart
df <- plot_CD276_PSCA_df %>%
  mutate(CD276_PSCA_CoExp = factor(CD276_PSCA_CoExp, levels = names(pie_colors))) %>%
  complete(CD276_PSCA_CoExp = names(pie_colors), fill = list(percent = 0)) %>%
  arrange(CD276_PSCA_CoExp) %>%
  mutate(
    ymax = cumsum(percent),
    ymin = lag(ymax, default = 0),
    ymid = (ymin + ymax) / 2,
    lab  = ifelse(percent > 0, paste0(CD276_PSCA_CoExp, " (", round(percent, 1), "%)"), "")
  )

p_CD276_PSCA <- ggplot(df, aes(x = 1, y = percent, fill = CD276_PSCA_CoExp)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y", clip = "off") +   # allow labels outside
  # Less horizontal span -> bigger pie; still some room for labels
  xlim(c(0.5, 1.5)) +
  scale_fill_manual(
    name   = "Co-expression",
    values = pie_colors,
    breaks = names(pie_colors),
    drop   = FALSE
  ) +
  theme_void(base_size = 14) +
  theme(
    legend.position  = "right",
    legend.title     = element_text(size = 16, face = "bold"),
    legend.text      = element_text(size = 14),
    legend.key.size  = unit(1.2, "cm"),
    plot.margin      = margin(10, 40, 10, 10)  # extra right margin for longer labels
  )
p_CD276_PSCA

ggsave(file.path(file.path(project_dir, "figures"), paste0("PDX_PSCA_NECTIN1/FigPDX_i", ".tiff")),
       plot = p_CD276_PSCA, device = "tiff",
       width = 120, height = 120, units = "mm", dpi = 600)
ggsave(file.path(file.path(project_dir, "figures"), paste0("PDX_PSCA_NECTIN1/FigPDX_i", ".pdf")),
       plot = p_CD276_PSCA, device = "pdf",
       width = 120, height = 120, units = "mm", dpi = 600)


