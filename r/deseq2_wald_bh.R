#!/usr/bin/env Rscript
#Horizontal Gene Transfer from Fungi to Myxozoa: Genomic Evidence and Evolutionary Implications		
#Amr Galal Abdelraheem Ibrahim; Edson Aparecido Adriano	

	
#Creator: Ibrahim A. G. A. (amrgalalibrahim@gmail.com)

suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
})

option_list <- list(
  make_option(c("--counts"), type="character", help="Counts matrix TSV (genes x samples) with header, first col gene_id"),
  make_option(c("--samples"), type="character", help="Sample metadata TSV with columns: sample, condition (and others)"),
  make_option(c("--design"), type="character", default="~ condition", help="Design formula (default: ~ condition)"),
  make_option(c("--contrast"), type="character", default=NULL,
              help="Contrast string 'factor,levelA,levelB' meaning levelA vs levelB (e.g., condition,infection,spore)"),
  make_option(c("--out_prefix"), type="character", default="deseq2", help="Output prefix"),
  make_option(c("--min_count"), type="integer", default=10, help="Filter genes with rowSums < min_count (default 10)")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$counts) || is.null(opt$samples)) {
  stop("Provide --counts and --samples", call.=FALSE)
}

counts <- read.delim(opt$counts, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
if (ncol(counts) < 2) stop("Counts matrix must have gene_id + >=1 sample columns")

gene_id <- counts[[1]]
count_mat <- as.matrix(counts[,-1])
rownames(count_mat) <- gene_id

meta <- read.delim(opt$samples, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
if (!all(c("sample","condition") %in% colnames(meta))) {
  stop("samples.tsv must include columns: sample, condition", call.=FALSE)
}
meta$condition <- as.factor(meta$condition)

# ensure sample order matches counts columns
if (!all(colnames(count_mat) %in% meta$sample)) {
  missing <- setdiff(colnames(count_mat), meta$sample)
  stop(paste("These samples in counts.tsv are missing from samples.tsv:", paste(missing, collapse=", ")), call.=FALSE)
}
meta <- meta[match(colnames(count_mat), meta$sample), , drop=FALSE]

# Filter low counts
keep <- rowSums(count_mat) >= opt$min_count
count_mat <- count_mat[keep, , drop=FALSE]

dds <- DESeqDataSetFromMatrix(
  countData = round(count_mat),
  colData = meta,
  design = as.formula(opt$design)
)

dds <- DESeq(dds, test="Wald")  # Wald test as described

# Results
if (!is.null(opt$contrast)) {
  parts <- strsplit(opt$contrast, ",")[[1]]
  if (length(parts) != 3) stop("--contrast must be 'factor,levelA,levelB'", call.=FALSE)
  res <- results(dds, contrast=parts)
} else {
  res <- results(dds)
}

# BH adjustment is default in DESeq2 (p.adjust.method="BH")
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[, c("gene_id", setdiff(colnames(res_df), "gene_id"))]

out_res <- paste0(opt$out_prefix, "_results.tsv")
write.table(res_df, file=out_res, sep="\t", quote=FALSE, row.names=FALSE)

# Normalized counts for plotting
norm_counts <- counts(dds, normalized=TRUE)
norm_df <- data.frame(gene_id=rownames(norm_counts), norm_counts, check.names=FALSE)
out_norm <- paste0(opt$out_prefix, "_normalized_counts.tsv")
write.table(norm_df, file=out_norm, sep="\t", quote=FALSE, row.names=FALSE)

# Session info for reproducibility
out_sess <- paste0(opt$out_prefix, "_sessionInfo.txt")
sink(out_sess)
cat("DESeq2 pipeline (Wald test + BH)\\n")
sessionInfo()
sink()

message("Wrote: ", out_res)
message("Wrote: ", out_norm)
message("Wrote: ", out_sess)
