library(matrixStats)
library(rstudioapi)
library(Seurat)
library(tidyverse)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
source("~/GitHub/shivangsharma/SpaCeMInD/R/file_manipulation.R")
source("~/GitHub/shivangsharma/BulkRNAseq/R/data_transformation.R")
source("~/GitHub/shivangsharma/BulkRNAseq/R/data_visualization.R")
source("~/GitHub/shivangsharma/SpaCeMInD/R/statistical_functions.R")

# QC ---------------------------------------------------------------------------
p1 <- DimPlot(object = obj_dedoublet, group.by = "Dataset", alpha = 0.2, raster = FALSE, label = TRUE)
p2 <- DimPlot(object = obj, group.by = "CoarseAnnotation", alpha = 0.2, raster = FALSE, label = TRUE)
obj_dedoublet <- subset(obj, cells = colnames(obj)[obj$`Doublet Status` %in% c("Doublet")], invert = TRUE)

project_dir <- getActiveDocumentContext()$path %>% dirname() %>% cd..()
data_obj_loc <- read.csv(file = file.path(project_dir, "data", "data_loc.txt"), header = FALSE)[1, 1]
obj <- readRDS(data_obj_loc)

plot_density <- function(object, feature, feature_title) {
  
  x <- obj@reductions$umap.cca@cell.embeddings[, 1]
  y <- obj@reductions$umap.cca@cell.embeddings[, 2]
  w <- FetchData(obj, assay = "RNA", layer = "data", vars = feature)[, feature]
  
  dens <- ks::kde(x = obj@reductions$umap.cca@cell.embeddings[, c(1, 2)], w = w / sum(w) * length(w), bgridsize = rep(1000, 2))
  ix <- findInterval(obj@reductions$umap.cca@cell.embeddings[, c(1)], dens$eval.points[[1]])
  iy <- findInterval(obj@reductions$umap.cca@cell.embeddings[, c(2)], dens$eval.points[[2]])
  ii <- cbind(ix, iy)
  z <- dens$estimate[ii]
  names(z) <- colnames(obj)
  obj$z <- z
  plt <- FeaturePlot(obj, 
                     features = c("z"), 
                     raster = FALSE,
                     pt.size = 0.0005, 
                     label = FALSE, 
                     label.color = "black") + 
    ggtitle(label = feature_title) + 
    theme(plot.title = element_text(vjust = - 3.5), 
          plot.margin = unit(c(-5, -5, -5, -5), units = "mm")) + 
    NoLegend() + 
    NoAxes()
  # plt <- rasterise(input = plt, layers = "Point", dpi = 600)
  
}

# Fig Supp 1: HLA Expression plots ---------------------------------------------

# Fig Supp 2: AR Reg. genes ----------------------------------------------------

# Distribution of datasets across celltypes
axes <- list(x = -12, y = -15, x_len = 5, y_len = 5)
umap_plot <-  DimPlot(obj, 
                      reduction = "umap.cca", 
                      label = TRUE, 
                      repel = TRUE, 
                      raster = FALSE, 
                      pt.size = 0.0005
) + 
  ggtitle("UMAP Clusters") + 
  annotate("segment", 
           x = axes$x, xend = axes$x + c(axes$x_len, 0), 
           y = axes$y, yend = axes$y + c(0, axes$y_len), 
           arrow = arrow(type = "closed", length = unit(10, 'pt')), 
           linewidth = 1
  ) + 
  annotate(geom = "text", x = axes$x + axes$x_len / 2, y = axes$y - 1.5, label = "UMAP1", fontface = "bold", size = 5) + 
  annotate(geom = "text", x = axes$x - 1.5, y = axes$y + axes$y_len / 2, label = "UMAP2", angle = "90", fontface = "bold", size = 5) + 
  theme(plot.title = element_text(hjust = 0.5, vjust = - 3.5), 
        axis.line = element_line(arrow = arrow()), 
        plot.margin = unit(c(-5, -5, -5, -5), units = "mm")) + 
  NoAxes() + 
  guides(col = guide_legend(ncol = 2, override.aes = list(size=6)))

datasets_PCa <- obj$Dataset_orig %>% unique()
UMAP_plt <- vector(mode = "list", length = length(datasets_PCa))
for (i in c(1:length(datasets_PCa))) {
  
  obj_sub <- subset(obj, subset = Dataset_orig == datasets_PCa[i])
  UMAP_plt[[i]] <-  DimPlot(obj_sub, 
                            reduction = "umap.cca", 
                            label = FALSE, 
                            repel = TRUE, 
                            raster = FALSE, 
                            pt.size = 0.0005
  ) + 
    ggtitle(datasets_PCa[i]) + 
    theme(plot.title = element_text(hjust = 0.5, vjust = - 3.5, size = 16), 
          plot.margin = unit(c(-5, -5, -5, -5), units = "mm")) + 
    NoLegend + 
    NoAxes()
  
}
design <- "
  11234
  11555
  11678
  99999
  abcde
"
gc()
plt <- umap_plot + UMAP_plt[[1]] + UMAP_plt[[2]] + UMAP_plt[[3]] + plot_spacer() + UMAP_plt[[4]] + UMAP_plt[[5]] + UMAP_plt[[6]] + plot_spacer() + guide_area() + UMAP_plt[[7]] + UMAP_plt[[8]] + UMAP_plt[[9]] + UMAP_plt[[10]] + plot_layout(design = design, guides = 'collect', heights = c(10, -1.03, 10, -1.03, 10)) & theme(plot.margin = unit(c(0, 0, 0, 0), units = "mm"))
plt
ggsave(file.path(project_dir, "figures", "Fig_Supp_3b.tiff"), device = "tiff", width = 400, height = 240, units = "mm", dpi = 250)




# Fig Supp 4 -------------------------------------------------------------------
# Identify highest scoring celltypes for AR and NE scores ---------------------

# Fig Supp 5 -------------------------------------------------------------------
# CD276 expression in Luminal vs Basal






# Fig Supp 2 -------------------------------------------------------------------
# Alternate color schemes and trash code ---------------------------------------
data_rel <- sweep(x = data, MARGIN = 2, STAT = data["PDCD1LG2", ] %>% as.integer(), FUN = "-")
data_rel <- data_rel[c("AR", "CD276", "KLK2", "KLK3", "TMPRSS2", "FOLH1", "TACSTD2", "CD111", "FKBP5", "CD274", "IGF1R", "PDCD1LG2", "PDCD1", "LAG3", "TIGIT", "TNFRSF4", "TNFRSF9", "CTLA4", "HLA-A", "HLA-B", "HLA-C", "SLC2A1", "DLL3"), ]
vsd <- vst(data[, order] %>% as.matrix())
dm <- model.matrix(~Pathology + Dataset, meta_data)
treatment.design <- dm[,1:4]
vsd <- limma::removeBatchEffect(vsd, batch = meta_data$Dataset)
luminal_basal_genes <- c("AR", "CD276", "KLK2", "KLK3", "TMPRSS2", "FOLH1", "TACSTD2", "CD111", "FKBP5", "CD274", "IGF1R", "PDCD1LG2", "PDCD1", "LAG3", "TIGIT", "TNFRSF4", "TNFRSF9", "CTLA4", "HLA-A", "HLA-B", "HLA-C", "SLC2A1", "DLL3")
vsd <- as.data.frame(vsd)
vsd_sub <- vsd[luminal_basal_genes, ]
vsd_sub_scaled <- vsd_sub*data_rel
samples_rem <- which(colnames(vsd_sub_scaled) %in% c("GSE172357_BPH_P4", "GSE172357_BPH_P13"))
vsd_sub <- vsd_sub_scaled[, - samples_rem]
pheatmap(vsd_sub_scaled %>% log1p(), cluster_rows = FALSE, cluster_cols = FALSE)
gene_cv <- rowVars(vsd_sub %>% as.matrix())/rowMeans(vsd_sub %>% as.matrix())*100

obj_epi <- subset(obj, subset = CellAnnFine %in% c("MSMB+", "Luminal", "Hillock", "Club", "Basal", "UK1", "UK2", "UK3"))
obj_epi <- NormalizeData(obj_epi)
temp <- obj_epi@assays$RNA@layers$counts
rownames(temp) <- rownames(obj_epi)
expr <- temp["CD276", ]
df <- data.frame(Expr = expr, CellType = obj_epi$CellAnnFine %>% as.character() %>% as.factor())
barplot(Expr~CellType,
        data=df,
        main="CD276 expression",
        xlab="Cell type",
        ylab="log-normalized counts",
        col=c("red", "orange", "green", "purple", "blue"),
        border="black", 
        outline = FALSE
)

df_mean <- c(mean(df$Expr[df$CellType %in% "MSMB+"]), mean(df$Expr[df$CellType %in% "Luminal"]), mean(df$Expr[df$CellType %in% "Basal"]), mean(df$Expr[df$CellType %in% "Hillock"]), mean(df$Expr[df$CellType %in% "Club"]), mean(df$Expr[df$CellType %in% "UK1"]), mean(df$Expr[df$CellType %in% "UK2"]), mean(df$Expr[df$CellType %in% "UK3"]))
df_sd <- c(sd(df$Expr[df$CellType %in% "MSMB+"]), sd(df$Expr[df$CellType %in% "Luminal"]), sd(df$Expr[df$CellType %in% "Basal"]), sd(df$Expr[df$CellType %in% "Hillock"]), sd(df$Expr[df$CellType %in% "Club"]), sd(df$Expr[df$CellType %in% "UK1"]), sd(df$Expr[df$CellType %in% "UK2"]), sd(df$Expr[df$CellType %in% "UK3"]))
df_names

df_summ <- data.frame(Celltype = df_names, Mean = df_mean, SD = df_sd)
ggplot(df_summ) +
  geom_bar( aes(x=Celltype, y=Mean), stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=Celltype, ymin=Mean-SD, ymax=Mean+SD), width=0.4, colour="orange", alpha=0.9, size=1.3)
df_summ <- df %>% group_by(CellType) %>% dplyr::summarize(Mean = mean(x = .), SD = sd(x = .))
ggplot(df) + geom_bar(aes(x = CellType, y = Expr), stat="identity", fill="skyblue", alpha=0.7)

# Fig Supp 3 -------------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(data.table)
source("~/GitHub/shivangsharma/CosMx/R/statistical_functions.R")

working_dir <- "/media/singhlab/B684-19A6/CosMxJuly2022Experiment/Single cell datasets/Prostate"
pattern_dataset_names <- "^GSE.*"
pattern_RData_names <- ".*_lean.RData"
genes <- c("CD276")
order_dataset <- c("GSE176031_T", "GSE176031_B", "GSE185344_T", "GSE185344_B", "GSE193337_T", "GSE193337_B", "GSE141445_T", "GSE172357_BPH", "GSE120716_N")

stacked_barplot <- function(working_dir,
                            pattern_dataset_names,
                            pattern_RData_names,
                            genes,
                            order_celltype,
                            order_dataset) {
  
  datasets <- dir(path = working_dir,
                  pattern = pattern_dataset_names,
                  full.names = FALSE)
  datasets_full_path <- file.path(working_dir, datasets)
  RData_full_path <- list.files(path = datasets_full_path,
                                pattern = pattern_RData_names,
                                full.names = TRUE)
  n_datasets <- length(RData_full_path)
  plt <- list()
  
  for(gene in genes) {
    
    print(gene)
    df_stacked_barplot <- data.table(matrix(nrow = 0, ncol = 4)) %>%
      setNames(., c("Cell_ID", "Celltype", "Dataset", "Expression"))
    
    for(i in c(1:n_datasets)) {
      
      print(datasets[i])
      obj <- readRDS(file = RData_full_path[i])
      obj <- NormalizeData(obj, normalization.method = "LogNormalize")
      Idents(obj) <- "CellAnn"
      df_stacked_barplot <- data.table(Cell_ID = colnames(obj),
                                       Celltype = Idents(obj),
                                       Dataset = datasets[i] %>% as.factor(),
                                       Expression = obj@assays$RNA$data[gene, ]) %>%
        rbind(df_stacked_barplot, .)
      
    }
    
    df_stacked_barplot$Dataset <- factor(df_stacked_barplot$Dataset,
                                         levels = order_dataset)
    # df_stacked_barplot <- df_stacked_barplot[df_stacked_barplot$Expression > 0, ]
    
    summarized_data <- df_stacked_barplot %>%
      group_by(Celltype, Dataset) %>%
      summarize(Expression = mean(log1p(Expression)),
                Count = n(),
                Fraction = sum(Expression > 0)/n())
    
    plt[[gene]] <- ggplot(data = summarized_data,
                          mapping = aes(x = Dataset,
                                        y = Expression,
                                        fill = Celltype)) +
      geom_bar(stat = "identity")
    plt[[gene]]
    
  }
  
}