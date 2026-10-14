#!/usr/bin/env Rscript
#Horizontal Gene Transfer from Fungi to Myxozoa: Genomic Evidence and Evolutionary Implications		
#Amr Galal Abdelraheem Ibrahim; Edson Aparecido Adriano	

	
#Creator: Ibrahim A. G. A. (amrgalalibrahim@gmail.com)

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("--install"), action="store_true", default=FALSE,
              help="If set, attempts to install missing packages (ape, phytools, ggtree)"),
  make_option(c("--out"), type="character", default="phylo_sessionInfo.txt", help="Output session info file")
)
opt <- parse_args(OptionParser(option_list=option_list))

pkgs_cran <- c("ape", "phytools")
pkgs_bioc <- c("ggtree")

if (opt$install) {
  for (p in pkgs_cran) {
    if (!requireNamespace(p, quietly=TRUE)) install.packages(p, repos="https://cloud.r-project.org")
  }
  if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")
  for (p in pkgs_bioc) {
    if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, ask=FALSE, update=FALSE)
  }
}

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(ggtree)
})

sink(opt$out)
cat("Phylogenetic/statistical R environment\\n")
sessionInfo()
sink()

message("Wrote: ", opt$out)
