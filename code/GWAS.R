geno <- read.csv(file = "~/Downloads/GCST90274713.tsv", sep = "\t")
geno$gene <- NA
gene_loc <- read.csv(file = "~/Downloads/mart_export.txt", sep = "\t")
colnames(gene_loc) <- c("start", "end", "chr", "gene")
rownames(gene_loc) <- gene_loc$gene
for(gene in gene_loc$gene) {
  print(gene)
  geno_loc <- which((geno$base_pair_location >= gene_loc[gene, "start"]) & (geno$base_pair_location <= gene_loc[gene, "end"]) & (geno$chromosome == gene_loc$chr))
  geno[geno_loc, "gene"] <- gene

}
geno_sub <- geno[!is.na(geno$gene), ]
geno_sub_sub <- geno_sub[geno_sub$p_value <= 5e-2, ]
geno_sub_sub$gene <- geno_sub_sub$gene %>% factor(levels = gene_loc$gene)
geno_summ <- geno_sub_sub %>% group_by(gene) %>% summarize(beta = mean(abs(beta)))
geno_summ$length <- gene_loc$end - gene_loc$start
geno_summ$norm_counts <- geno_summ$counts/geno_summ$length
