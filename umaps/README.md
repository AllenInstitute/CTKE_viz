## UMAP Neighborhood plots

Authors: Cindy van Velthoven, Ray Sanchez

This folder contains script, "make_all_umap.R", that was adapted from Cindy's original UMAP plotting code (in "Cindy_originalcode" subfolder). The script requires
the following input files:

* **standard metadata file** as described in the cluster metrics subfolder
* **umap coordinates** - data frame with 3 columns - sample_name, UMAP coord 1, and UMAP coord 2. sample names must match those in the metadata file
* **taxonomy file** - standard annotated cluster table providing id, color, and label for every node in the taxonomy. Used as an extra "sanity check" to complement the metadata file, and is especially important if plotting a dataset that may be missing a node present in the taxonomy
