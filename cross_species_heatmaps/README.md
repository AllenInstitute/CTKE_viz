## Cross-species heatmaps

Authors: Ray Sanchez, Nik Jorstad

This folder contains scripts and files necessary to create the clickable cross-species heatmaps shown at the cluster level on cell type cards. This information is also
provided to the ontology team to support their development work on linkages between species taxonomimes. The subfolder "input files" contains the following files/file types:

* **species_CrossSpecies_Similarity.csv** - this file contains pairwise Euclidean distance metrics comparing every cluster in "species" to every cluster in the other two species. Given for all 3 species
* **species1_to_species2_clusters.csv** - contains information about homologous cell types between two species (for example, human_to_mouse_clusters). Manually built from supplementary tables in Bakken et al., 2021 and the **species_cross-species.csv** files also present in this folder
* **nik_script.R** - Nik's script to calculate Euclidean distances. Used to build species_CrossSpecies_Similarity files

The main directory contains the following R scripts:
* **build_clickable_heatmaps.R** - the main function for plotting heatmaps. Requires files from "input files" and necessary input parameters are detailed in the script
* **plot_example.R** - example script showing how all heatmaps for human M1 cell types can be generated with a few lines of code
