## Spatial transcriptomics panel

Authors: Brian Long (with some modifications by Ray Sanchez)

This folder contains a Python notebook to make example location plots and cortical depth histograms from the MERFISH data collected by Xiaowei Zhuang's lab. 
It requires a counts file in .h5ad format and metadata file with spatial coordinates for each cell (both files are too large to include here). Brian has written a nice description
of the code and its requirements [here](https://github.com/AllenInstitute/celltype_cards_spatial/tree/85a99a5973fde9a9567ffd2e7b16c9c61eeac0a1).

Running this code is a bit more
manual than some of the others, but TBD whether it's worth improving, as this will not be the way spatial transcriptomics data are displayed on the whole mouse brain explorer
(and likely also the internal version of cell type cards). 
