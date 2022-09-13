## Dot plots for Transcriptomics

Authors: Ray Sanchez

This folder contains a script with functions used to make all of the dot plots for cell type cards. The main plotting function draws heavily from the [scrattch.vis](https://github.com/AllenInstitute/scrattch.vis)
package for RNA-seq data visualization. The script "example.R" shows how all dot plots for a taxonomy can be generated with one function call. It requires the following inputs:

* **count matrix** - standard cellxgene matrix with sample names as column names and genes as row names. Can be read into R from any format (.mtx, .RDS, anndata, etc.) but function accepts a sparse matrix
* **metadata file** - standard annotations file described in the [cluster metrics folder](https://github.com/AllenInstitute/celltype_cards_contenthub/tree/main/all_code/cluster_metrics)
* **genes to plot** - example given here is "AllHumanGenes_Tidy.csv", which provides the NS-Forest markers for the Human M1 10X data. This file represents the standard file format necessary for "make_all_dots.R" to work. 
* **taxonomy file** - standard taxonomy file described elsewhere

**Note:** the current plotting method assumes a predetermined list of marker genes is available. For internal versions of cell type cards, this function will need to be updated such that a set of marker genes can be calculated on the fly. 
