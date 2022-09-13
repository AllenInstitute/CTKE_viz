# Read in files provided in NeMO archives for plotting

suppressPackageStartupMessages({
  library(plyr)
  library(tidyverse)
  library(Matrix)
})

# Load rownames (features) and colnames (barcodes)
features <- read.table("features.tsv", sep="\t", header = FALSE)
barcodes <- read.table("barcode.tsv",sep = ",", header = TRUE)
colnames(features) <- c("ens", "gene","dimname")
colnames(barcodes) <- c("num", "sample_name")

# Load .mtx file and set dimnames
count_mat <- readMM("matrix.mtx")
count_mat@Dimnames[[1]] <- features$gene
count_mat@Dimnames[[2]] <- barcodes$sample_name