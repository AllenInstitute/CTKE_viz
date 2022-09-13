## **Cluster Metrics**
**Authors:** Ray Sanchez


This subfolder contains the following files:

* **Make_cluster_metrics.R** - primary function used to generate all cluster metrics plots for a taxonomy. Requires functions in "CM_helper_functions.R"
* **CM_helper_functions.R** - contains 2 functions to make individual violin plots, with and without highlighting the type featured on the cell type card
* **plot_all_CM.R** - example script showing how all plots can be generated in a few lines of code using the above functions
* **cl_updated.df.csv** - standard annotated cluster table providing id, color, and label for every node in the taxonomy. Used as an extra "sanity check" to complement the metadata file, and is especially important if plotting a dataset that may be missing a node present in the taxonomy 
* **scRNA10_xv2_A Metadata.csv** - standard metadata file used for plotting. Same format as what's used in shiny and requires the following columns:
    * sample_name
    * _id, _color, and _label columns for each class, subclass and type in the taxonomy
    * Column with number of genes detected, named as "genes_detected_label"
    * Column with number of UMIs, named as "unique.counts"

Note: The main functions could be easily adapted to handle any metadata available (donor age, confusion, etc.) but currently only accept genes detected and UMIs. 
