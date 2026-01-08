library(readxl)
library(writexl)
library(stringr)
source("~/GitHub/Shivang/CosMx/R/data_transformation_functions.R")

calc_CI_spearman <- function(n_obs, r, alpha) {
  
  stderr <- 1.0 / sqrt(n_obs - 3)
  delta <- 1.96 * stderr
  lower <- tanh(atanh(r) - delta)
  upper <- tanh(atanh(r) + delta)
  c(lower, upper)
  
}

data <- read_xlsx(path = "/media/singhlab/B684-19A6/Nature PCAN/data/rt-pcr.xlsx", col_names = TRUE)
cor.test(x = data$`B7-H3`, y = data$AR, method = "spearman")
n_obs <- 33
r <- 0.6580882
ci_interval <- calc_CI_spearman(n_obs = n_obs, r = r, alpha = 0.95)

cor.test(x = data$`B7-H3`, y = data$`PD-L1`, method = "spearman")
n_obs <- 33
r <- 0.5855615
ci_interval <- calc_CI_spearman(n_obs = n_obs, r = r, alpha = 0.95)

cor.test(x = data$`PD-L1`, y = data$AR, method = "spearman")
n_obs <- 33
r <- 0.375
ci_interval <- calc_CI_spearman(n_obs = n_obs, r = r, alpha = 0.95)
plot(x = data$`AR`, y = data$`B7-H3`)


data <- read_xlsx(path = "/media/singhlab/B684-19A6/Nature PCAN/data/cell_line_prot_data.xlsx", sheet = "Full protein matrix")
data <- data %>% as.data.frame()
colnames_data <- data[1, ]
colnames(data) <- colnames_data
rownames_data <- data[, 1] %>% as.matrix()
rownames(data) <- rownames_data
data <- data[c(-1), ]
data <- data[, c(-1)]
data <- df_transpose(data)
write_xlsx(x = data, path = "/media/singhlab/B684-19A6/Nature PCAN/data/cell_line_prot_data_formatted.xlsx")

rownames(data) <- rownames(data) %>% 
  str_remove(string = ., pattern = ".*;") %>% 
  str_remove(string = ., pattern = "_HUMAN")
colnames(data) <- colnames(data) %>% 
  str_remove(string = ., pattern = ".*;")
write_xlsx(x = data, path = "/media/singhlab/B684-19A6/Nature PCAN/data/cell_line_prot_data_formatted.xlsx")

genes <- c("PRSS24", "STEAP1", "STEAP2", "PSCA", "FOLH1", "TROP", "PVRL1", "DLL3", "CD276", "CD274", "PDCD1LG2", "PDCD1", "LAG3", "TIGIT", "TNFRSF4", "TNFRSF9", "CTLA4", "AR", "KLK2", "KLK3", "TMPRSS2", "FKBP5", "IGF1R")
cell_lines_prostate <- c("VCaP", "PC-3", "PWR-1E", "LNCaP-Clone-FGC", "22RV1", "DU-145", "BPH-1", "NCI-H660")
data_subset <- data[genes, cell_lines_prostate]
data$Gene <- rownames(data)
