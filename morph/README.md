## Morphology

Authors: Alice Mukora, Ray Sanchez

This folder contains 2 sets of scripts. The first is a set of Python notebooks written by Alice Mukora that generate average cortical depth histograms and example morphology
reconstructions from the Tolias lab Patch-seq data. The second is two R scripts used to arrange the reconstructions into a grid for display on cell type cards. The most important
of these 2 scripts is "arrange_morph_SVGs.R", which contains a function to generate grids for all cell type cards given the following inputs:

* **ctc_summary.csv** - this file is output by Alice's code and provides filenames, brief text descriptions, input .swc files, and input cell IDs for all reconstructions
* **Panel_Image Info.csv** - for morphology panel. This is described in detail in the [content placement CSV scripts](https://github.com/AllenInstitute/celltype_cards_contenthub/tree/main/all_code/createCSVs_productteam) but is used here to determine which files should be included in each grid. (Note: this could probably be improved/generalized to not rely on Panel_Image Info, but was a late add to cell type cards in response to user feedback, and this was the easiest way to accomplish it).
* **morphology metadata CSV** - this can be used to order placement of reconstructions by different features (i.e. cortical depth) but for now I use taxonomy order to do this
* **taxonomy file** - standard cluster annotation table referenced in all other scripts. Used to determine order of images (subclass or cluster ID)


"arrange_reconstructions.R" is an older version of "arrange_morph_SVGs.R" that is specific to this dataset and has not been functionized. 
