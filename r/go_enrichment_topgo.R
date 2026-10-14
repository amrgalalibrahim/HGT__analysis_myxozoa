#!/usr/bin/env Rscript
#Horizontal Gene Transfer from Fungi to Myxozoa: Genomic Evidence and Evolutionary Implications		
#Amr Galal Abdelraheem Ibrahim; Edson Aparecido Adriano	

	
#Creator: Ibrahim A. G. A. (amrgalalibrahim@gmail.com)

suppressPackageStartupMessages({
  library(optparse)
  library(topGO)
})

option_list <- list(
  make_option(c("--gene_universe"), type="character", help="Text file: all genes, one per line"),
  make_option(c("--sig_genes"), type="character", help="Text file: significant genes, one per line"),
  make_option(c("--gene2go"), type="character", help="TSV: gene_id<TAB>GO:0000000;GO:0000001"),
  make_option(c("--ontology"), type="character", default="BP", help="GO ontology: BP/MF/CC (default BP)"),
  make_option(c("--out_prefix"), type="character", default="topgo", help="Output prefix"),
  make_option(c("--algorithm"), type="character", default="weight01", help="topGO algorithm (default weight01)"),
  make_option(c("--statistic"), type="character", default="fisher", help="Test statistic (default fisher)")
)
opt <- parse_args(OptionParser(option_list=option_list))

universe <- readLines(opt$gene_universe)
sig <- readLines(opt$sig_genes)
sig_set <- setNames(as.integer(universe %in% sig), universe)

# gene2go mapping
g2g <- read.delim(opt$gene2go, header=FALSE, stringsAsFactors=FALSE, sep="\\t")
colnames(g2g) <- c("gene_id","go_terms")

geneID2GO <- list()
for (i in seq_len(nrow(g2g))) {
  gid <- g2g$gene_id[i]
  terms <- unlist(strsplit(g2g$go_terms[i], ";"))
  terms <- terms[terms != ""]
  geneID2GO[[gid]] <- terms
}

GOdata <- new("topGOdata",
              ontology=opt$ontology,
              allGenes=sig_set,
              annot=annFUN.gene2GO,
              gene2GO=geneID2GO)

result <- runTest(GOdata, algorithm=opt$algorithm, statistic=opt$statistic)

tab <- GenTable(GOdata,
                Fisher=result,
                orderBy="Fisher",
                topNodes=50)

out_tab <- paste0(opt$out_prefix, "_", opt$ontology, "_topGO.tsv")
write.table(tab, file=out_tab, sep="\\t", quote=FALSE, row.names=FALSE)

sink(paste0(opt$out_prefix, "_topgo_sessionInfo.txt"))
sessionInfo()
sink()

message("Wrote: ", out_tab)
