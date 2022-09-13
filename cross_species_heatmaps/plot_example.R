# Example of using "build_clickable_heatmaps.R" to plot some human cell types

source("build_clickable_heatmaps.R")

# Load taxonomy info
humancldf <- read_csv("HumanTaxonomy_cldf.csv")
marmosetcldf <- read_csv("MarmosetTaxonomy_cldf.csv")
mousecldf <- read_csv("cl_updated.df.csv")


## Load in cell type homology tables
human_to_mouse <- read_csv("human_to_mouse_clusters.csv")
human_to_marmoset <- read_csv("human_to_marmoset_clusters.csv")



# Load confusion matrix
human_confusion <- read_csv("Human_CrossSpecies_Similarity.csv")

# Make plots

for (i in 1:length(humancldf$cluster_label)) {
  plot_cross_species_heatmaps(main_confusion = human_confusion,
                              main_cldf = humancldf,
                              node = humancldf$cluster_label[i],
                              spec1 = "marmoset",
                              spec1_tbl = human_to_marmoset,
                              spec1_cldf = marmosetcldf,
                              spec2 = "mouse",
                              spec2_tbl = human_to_mouse,
                              spec2_cldf = mousecldf)
  
}
