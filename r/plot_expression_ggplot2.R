#!/usr/bin/env Rscript
#Horizontal Gene Transfer from Fungi to Myxozoa: Genomic Evidence and Evolutionary Implications		
#Amr Galal Abdelraheem Ibrahim; Edson Aparecido Adriano	

	
#Creator: Ibrahim A. G. A. (amrgalalibrahim@gmail.com)

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(reshape2)
})

option_list <- list(
  make_option(c("--norm_matrix"), type="character", help="Normalized counts TSV from DESeq2 (gene_id + samples)"),
  make_option(c("--samples"), type="character", help="Sample metadata TSV with columns: sample, condition"),
  make_option(c("--out_prefix"), type="character", default="expr", help="Output prefix"),
  make_option(c("--top_n"), type="integer", default=50, help="Top N variable genes for heatmap (default 50)")
)
opt <- parse_args(OptionParser(option_list=option_list))

nm <- read.delim(opt$norm_matrix, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
gene_id <- nm[[1]]
mat <- as.matrix(nm[,-1])
rownames(mat) <- gene_id

meta <- read.delim(opt$samples, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
meta$condition <- as.factor(meta$condition)
meta <- meta[match(colnames(mat), meta$sample), , drop=FALSE]

# log transform
log_mat <- log2(mat + 1)

# pick top variable genes
vars <- apply(log_mat, 1, var, na.rm=TRUE)
top_genes <- names(sort(vars, decreasing=TRUE))[seq_len(min(opt$top_n, length(vars)))]
hm <- log_mat[top_genes, , drop=FALSE]

# Convert to long
df <- data.frame(gene_id=rownames(hm), hm, check.names=FALSE)
df_long <- reshape2::melt(df, id.vars="gene_id", variable.name="sample", value.name="log2_norm_count")
df_long <- merge(df_long, meta, by.x="sample", by.y="sample", all.x=TRUE)

# Heatmap
p_hm <- ggplot(df_long, aes(x=sample, y=gene_id, fill=log2_norm_count)) +
  geom_tile() +
  theme_minimal(base_size=10) +
  theme(axis.text.y=element_text(size=6)) +
  labs(x="Sample", y="Gene", fill="log2(norm+1)", title="Expression heatmap (top variable genes)")

ggsave(paste0(opt$out_prefix, "_heatmap.pdf"), p_hm, width=10, height=8)

# Mean expression per condition (line plot)
df_mean <- aggregate(log2_norm_count ~ gene_id + condition, data=df_long, FUN=mean)
p_line <- ggplot(df_mean, aes(x=condition, y=log2_norm_count, group=gene_id)) +
  geom_line(alpha=0.3) +
  theme_minimal(base_size=10) +
  labs(x="Condition", y="Mean log2(norm+1)", title="Gene expression profiles by condition")

ggsave(paste0(opt$out_prefix, "_profiles.pdf"), p_line, width=10, height=6)

# Session info
sink(paste0(opt$out_prefix, "_plot_sessionInfo.txt"))
sessionInfo()
sink()

message("Wrote plots: ", opt$out_prefix, "_heatmap.pdf and _profiles.pdf")
