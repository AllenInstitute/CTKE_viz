# This is an example of how all cluster metrics plots for a taxonomy can be generated in a few lines of code

# Load scripts
source("make_cluster_metrics.R")

# Load metadata file
allmet <- read_csv("scRNA_10xv2_A Metadata.csv")

# Load taxonomy file
clupdate <- read_csv("cl_updated.df.csv")

# Set output directory
outdir <- "~/automated_cards"


# Make all cluster metrics plots
make_cluster_metrics(metadata = allmet,
                     taxonomy_file = clupdate,
                     outdir = outdir,
                     datasetname = "scRNA10xV2")
