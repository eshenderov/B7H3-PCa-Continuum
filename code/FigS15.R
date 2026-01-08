library(synergyfinder)
library(tidyverse)

syn <- readRDS("/projects/eshenderov-hpc/Shivang/projects/B7-H3 ADC Manuscript/data/SynergyFinder_r_object_2025-11-25.rds")

df <- data.frame(Enzalutamide = syn[["synergy_scores"]][["conc1"]], 
                 MGC018 = syn[["synergy_scores"]][["conc2"]], 
                 'ZIP Synergy' = syn[["synergy_scores"]][["ZIP_synergy"]],
                 'HSA Synergy' = syn[["synergy_scores"]][["HSA_synergy"]], 
                 'Loewess Synergy' = syn[["synergy_scores"]][["Loewe_synergy"]], 
                 'Bliss Synergy' = syn[["synergy_scores"]][["Bliss_synergy"]]
)



write.csv2(x = df, file = '/projects/eshenderov-hpc/Shivang/projects/B7-H3 ADC Manuscript/data/synergy_res.csv', col.names = TRUE, row.names = FALSE)
write.csv(x = df, file = '/projects/eshenderov-hpc/Shivang/projects/B7-H3 ADC Manuscript/data/synergy_res.csv', row.names = FALSE)




data("mathews_screening_data")
df_avg_conc <- read.csv(file = "/projects/eshenderov-hpc/Shivang/projects/B7-H3 ADC Manuscript/data/Matrix average.csv")
rownames(df_avg_conc) <- df_avg_conc$X
df_avg_conc$X <- NULL
colnames(df_avg_conc) <- c(0, 2500, 5000, 10000, 20000)
df_avg_conc$conc_r <- rownames(df_avg_conc)

df_long <- df_avg_conc %>%
  pivot_longer(
    cols = c("0", "2500", "5000", "10000", "20000"), # Selects columns starting with "Score_"
    names_to = "conc_c",           # Name of the new column for original column names
    values_to = "response"          # Name of the new column for values
  )
df_long$drug_row <- "MGC018"
df_long$drug_col <- "Enzalutamide"
df_long$conc_r_unit <- "nM"
df_long$conc_c_unit <- "nM"
df_long$block_id <- 1
df_long$conc_c <- as.numeric(df_long$conc_c)
df_long$conc_r <- as.numeric(df_long$conc_r)


res <- ReshapeData(
  data = df_long,
  data_type = "viability",
  impute = TRUE,
  impute_method = NULL,
  noise = TRUE,
  seed = 1)
res <- CalculateSynergy(
  data = res,
  method = c("ZIP", "HSA", "Bliss", "Loewe"),
  Emin = NA,
  Emax = NA,
  correct_baseline = "non")

res <- CalculateSensitivity(
  data = res,
  correct_baseline = "non"
)


sensitive_columns <- c(
  "block_id", "drug1", "drug2",
  "ic50_1", "ic50_2",
  "ri_1", "ri_2",
  "css1_ic502", "css2_ic501", "css")
res$drug_pairs[, sensitive_columns]


for (i in c(1, 2)){
  PlotDoseResponseCurve(
    data = res,
    plot_block = 1,
    drug_index = i,
    plot_new = FALSE,
    record_plot = FALSE
  )
}


Plot2DrugHeatmap(
  data = res,
  plot_block = 1,
  drugs = c(1, 2),
  plot_value = "response",
  dynamic = FALSE,
  summary_statistic = c("mean",  "median")
)

Plot2DrugHeatmap(
  data = res,
  plot_block = 1,
  drugs = c(1, 2),
  plot_value = "Bliss_synergy",
  dynamic = FALSE,
  summary_statistic = c( "quantile_25", "quantile_75")
)




Plot2DrugSurface(
  data = res,
  plot_block = 1,
  drugs = c(1, 2),
  plot_value = "response",
  dynamic = FALSE,
  summary_statistic = c("mean", "quantile_25", "median", "quantile_75")
)
Plot2DrugSurface(
  data = res,
  plot_block = 1,
  drugs = c(1, 2),
  plot_value = "Bliss_synergy",
  dynamic = TRUE,
  summary_statistic = c("mean", "quantile_25", "median", "quantile_75")
)
