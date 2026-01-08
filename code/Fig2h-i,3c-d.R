library(data.table)

# ---------- Params ----------
epithelial_types <- c("Luminal", "Neuroendocrine")
genes <- genes  # your character vector of features

# ---------- 1) Subset once & normalize once ----------
obj_epithelial <- subset(
  obj,
  subset = CoarseAnnotation %in% epithelial_types
)

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
  TROP2   = sum(TROP2,   na.rm = TRUE), 
  `B7-H3` = sum(`B7-H3`, na.rm = TRUE), 
  PSCA    = sum(PSCA,    na.rm = TRUE), 
  DLL3    = sum(DLL3,    na.rm = TRUE),
  N       = .N
), by = .(Dataset, Sample, Subtype)]

df_exp_epi_mean <- dt[, .(
  NECTIN1 = mean(NECTIN1, na.rm = TRUE),
  NECTIN4 = mean(NECTIN4, na.rm = TRUE), 
  TROP2 = mean(TROP2, na.rm = TRUE), 
  `B7-H3` = mean(`B7-H3`, na.rm = TRUE), 
  PSCA = mean(PSCA, na.rm = TRUE), 
  DLL3 = mean(DLL3, na.rm = TRUE)
), by = .(Dataset, Sample, Subtype)]


obj_epithelial2 <- NormalizeData(obj_epithelial, 
                                normalization.method = "RC",
                                scale.factor = 100000)
df_exp_epi2 <- data.frame(FetchData(object = obj_epithelial2, 
                               vars = c("NECTIN1", "NECTIN4", "CD276", "PSCA", "DLL3", "CoarseAnnotation", "Dataset", "Sample_ID", "Subtype"), 
                               assay = "RNA", 
                               layer = "data")
)
colnames(df_exp_epi2) <- c("NECTIN1", "NECTIN4", "B7-H3", "PSCA", "DLL3", "Celltype", "Dataset", "Sample", "Subtype")

df_exp_epi2 <- df_exp_epi2 %>% filter(Celltype %in% c("Luminal", "Neuroendocrine"))

dt2 <- as.data.table(df_exp_epi2)

df_exp_epi_summ2 <- dt[, .(
  NECTIN1 = mean(NECTIN1, na.rm = TRUE),
  NECTIN4 = mean(NECTIN4, na.rm = TRUE), 
  `B7-H3` = mean(`B7-H3`, na.rm = TRUE), 
  PSCA = mean(PSCA, na.rm = TRUE), 
  DLL3 = mean(DLL3, na.rm = TRUE)
), by = .(Dataset, Sample, Subtype)]



bin1 <- !(df_exp_epi == 0)
bin1[bin1] <- 1
bin1 <- as.data.frame(bin1)
bin2 <- !(df_exp_epi2 == 0)
bin2[bin2] <- 1
bin2 <- as.data.frame(bin2)

defected <- rownames(bin1)[!(bin1$NECTIN1 == bin2$NECTIN1)]
raw_ct <- FetchData(obj_epithelial, vars = "NECTIN1", assay = "RNA", layer = "counts")



genes_sub <- c("STEAP1", "STEAP2", "PSCA", "FOLH1", "KLK2", "TACSTD2", "NECTIN1", "NECTIN4", "DLL3", "CD276", "FKBP5", "IGF1R", "AR", "TMPRSS2")

mat_corr_sc <- vsd_sub[genes_sub, ] %>% t() %>% cor()
pheatmap(mat_corr_sc)

mat_corr_pdx <- expr[genes_sub, ] %>% t() %>% cor()
pheatmap(mat_corr_pdx)


mat_sc_scaled <- t(scale(t(vsd_sub)))
mat_pdx_scaled <- t(scale(t(expr[genes, ])))


df2 <- mat_pdx_scaled %>% 
  as.data.frame() %>% 
  mutate(gene = rownames(.)) %>% 
  select(gene, everything())

# Convert to long format: one row per gene×sample
df_long <- df2 %>% 
  pivot_longer(cols = -gene,
               names_to = "sample",
               values_to = "expr")

# Now plot CDFs for each gene
ggplot(df_long, aes(x = expr, colour = gene)) +
  stat_ecdf(geom = "step") +
  labs(x = "Expression (scaled)", y = "ECDF",
       title = "CDFs of gene expression across samples for each gene") +
  theme_minimal() +
  theme(legend.position = "right")




# Z: numeric matrix, genes in rows, samples in cols (scaled already)
# link: "probit" (pnorm) or "logit" (plogis)
# tau: softness of the mapping; larger = softer (less extreme)
# clamp: optional clipping to avoid extreme tails saturating the link
fuzzy_or_and_weighted <- function(Z,
                                  link = c("probit", "logit", "identity"),
                                  tau = 1,
                                  clamp = 4,
                                  scale_rows = FALSE,
                                  diag_to_na = TRUE,
                                  return_penalties = TRUE,
                                  cat = NULL,
                                  category_weight = NULL) {
  stopifnot(is.matrix(Z))
  link <- match.arg(link)
  n <- ncol(Z)
  
  if (isTRUE(scale_rows)) {
    Z <- t(scale(t(Z)))
  }
  
  # Construct μ – the membership/expression values
  if (link != "identity") {
    if (!is.null(clamp)) {
      Z[Z >  clamp] <-  clamp
      Z[Z < -clamp] <- -clamp
    }
    mu <- if (link == "probit") {
      pnorm(Z / tau)
    } else {
      plogis(Z / tau)
    }
  } else {
    # identity link: rescale rows to [0,1]
    row_min <- apply(Z, 1, min, na.rm = TRUE)
    row_max <- apply(Z, 1, max, na.rm = TRUE)
    mu <- (Z - row_min) / (row_max - row_min)
  }
  
  # Construct sample weights based on categories
  if (!is.null(cat)) {
    stopifnot(length(cat) == n)
    cats   <- unique(cat)
    K      <- length(cats)
    if (is.null(category_weight)) {
      w_cat <- rep(1 / K, K)
      names(w_cat) <- cats
    } else {
      stopifnot(length(category_weight) == K)
      names(category_weight) <- cats
      w_cat <- category_weight / sum(category_weight)
    }
    # sample weights = category weight divided equally among samples in category
    w_samp <- sapply(cat, function(cn) {
      w_cat[[as.character(cn)]] / sum(cat == cn)
    })
  } else {
    w_samp <- rep(1 / n, n)
  }
  
  # Weighted mean per gene
  m <- rowSums(mu * matrix(w_samp, nrow = nrow(mu), ncol = ncol(mu), byrow = TRUE), na.rm = TRUE)
  
  # Weighted co‐mean P_{g,h}
  Wsqrt <- sqrt(w_samp)
  mu_w  <- mu * matrix(Wsqrt, nrow = nrow(mu), ncol = ncol(mu), byrow = TRUE)
  P     <- tcrossprod(mu_w)  # gives sum_i (w_i * mu_{g,i} * mu_{h,i})
  
  AND <- P
  OR  <- outer(m, m, "+") - P
  
  if (diag_to_na) {
    diag(AND) <- NA_real_
    diag(OR)  <- NA_real_
  }
  
  out <- list(AND = AND, OR = OR)
  if (return_penalties) {
    out$AND_unfavorable <- 1 - AND
    out$OR_unfavorable  <- 1 - OR
  }
  return(out)
}


# --- Example usage ---
# res <- fuzzy_or_and(Z, link = "probit", tau = 1)
# AND <- res$AND
# OR  <- res$OR

and_or <- fuzzy_or_and_weighted(Z = mat_sc_scaled, 
                       link = "identity", 
                       cat = meta_data$Subtype, 
                       category_weight  = NULL)
plt <- pheatmap(
  and_or$AND[genes_sub, genes_sub],
  display_numbers   = TRUE,
  number_format     = "%.2f",
  number_color      = "black",
  fontsize_number   = 10,
  name = "AND score",
  heatmap_legend_param = list(
    title = "AND score",     # legend title
    title_gp = gpar(fontsize = 12, fontface = "bold")
  )
)
pdf(file.path(project_dir, "figures", "Fig_sc_bispecific_AND.pdf"),
    width  = 240/25.4,   # ~9.45 inches
    height = 240/25.4)   # ~9.45 inches
ComplexHeatmap::draw(plt)
dev.off()
plt <- pheatmap(
  and_or$OR[genes_sub, genes_sub],
  display_numbers   = TRUE,
  number_format     = "%.2f",
  number_color      = "black",
  fontsize_number   = 10,
  name = "AND score",
  heatmap_legend_param = list(
    title = "AND score",     # legend title
    title_gp = gpar(fontsize = 12, fontface = "bold")
  )
)
pdf(file.path(project_dir, "figures", "Fig_sc_bispecific_OR.pdf"),
    width  = 240/25.4,   # ~9.45 inches
    height = 240/25.4)   # ~9.45 inches
ComplexHeatmap::draw(plt)
dev.off()



and_or <- fuzzy_or_and_weighted(Z = mat_pdx_scaled, 
                       link = "identity", 
                       cat = meta_pdx$Subtype, 
                       category_weight  = NULL)
plt <- pheatmap(
  and_or$AND[genes_sub, genes_sub],
  display_numbers   = TRUE,
  number_format     = "%.2f",
  number_color      = "black",
  fontsize_number   = 10,
  name = "AND score",
  heatmap_legend_param = list(
    title = "AND score",     # legend title
    title_gp = gpar(fontsize = 12, fontface = "bold")
  )
)
pdf(file.path(project_dir, "figures", "Fig_pdx_bispecific_AND.pdf"),
    width  = 240/25.4,   # ~9.45 inches
    height = 240/25.4)   # ~9.45 inches
ComplexHeatmap::draw(plt)
dev.off()
plt <- pheatmap(
  and_or$OR[genes_sub, genes_sub],
  display_numbers   = TRUE,
  number_format     = "%.2f",
  number_color      = "black",
  fontsize_number   = 10,
  name = "OR score",
  heatmap_legend_param = list(
    title = "OR score",     # legend title
    title_gp = gpar(fontsize = 12, fontface = "bold")
  )
)
pdf(file.path(project_dir, "figures", "Fig_pdx_bispecific_OR.pdf"),
    width  = 240/25.4,   # ~9.45 inches
    height = 240/25.4)   # ~9.45 inches
ComplexHeatmap::draw(plt)
dev.off()


