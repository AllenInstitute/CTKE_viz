# Example of how to make all dot plots

library(plyr)
library(tidyverse)
library(Matrix)
library(feather)
library(scrattch.hicat)
library(scrattch.vis)
library(ggpubr)
library(ggfittext)

# Set dir
mydir <- "//allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/human"
setwd(mydir)

# Load necessary input files
taxonomy_file <- read_csv("HumanTaxonomy_cldf.csv")

tenxcounts <- readRDS("Human_M1_10xV3_Matrix.RDS")
tenxcounts <- as(tenxcounts, "dgTMatrix")
tenxmeta <- read_feather("Human_M1_10xV3_Metadata.feather")   
colnames(tenxmeta)[1] <- "sample_name"

allgenes_tidy <- read_csv("AllHumanGenes_Tidy.csv")

# Load plotting functions
source("make_dotplot.R")
source("make_all_dots.R")


# Generate all dot plots
make_all_dots(counts = tenxcounts,
              metadata = tenxmeta,
              taxonomy_file = taxonomy_file,
              genes = allgenes_tidy,
              outdir = mydir,
              datasetname = "Human10X")