# CTKE_viz
Code used to generate data visualizations and text for the Cell Type Knowledge Explorer.

# Overview
This folder contains subfolders with scripts used to create all the content in the Cell Type Knowledge Explorer. This repo does NOT provide an out-of-the-box package for plotting the data modalities described below, but the provided scripts can be used to reproduce any content displayed on the current version of the CTKE. Many of these scripts may be generalizable to
other products/analyses, and those that are currently not in this state may be improved upon with some minimal focused effort (see individual subfolders). 

# Subfolders

**cluster_metrics**
 - Scripts used to create all the cluster metrics violin plots in the cards

**createCSVs_productteam**
 - Scripts used to create the CSVs with necessary metadata for the product team to populate cell type cards

**cross_species_heatmaps**
 - Scripts and input files used to build clickable heatmap SVGs on cell type cards

**ephys**
 - Contains an R script used to plot summary electrophysiology metrics for all taxons in Patch-seq dataset
 - Contains a series of iPython notebooks used to read in .NWB files and plot example action potential traces from Patch-seq

**load_metadata**
 - Contains helper scripts that load in files from NeMO and update transcriptomic dataset annotations. Specific to this effort and probably not necessary in the future.

**morph**
 - contains a series of iPython notebooks used to generate morphological reconstructions from .SWC files and plot histograms
   of cortical depth of reconstructions, and a script to arrange reconstructions saved in SVG format into a labeled grid

**rnaseq_dotplots**
 - Contains scripts used to make dot plots for all of the cell type cards

**spatial_loc**
 - Contains one iPython notebook used to generate locations of cell types and histograms of cortical depth from MERFISH data

**text_summaries**
 - Contains R scripts used to generate short text summaries shown at the top of cell type cards. 

**umaps**
 - Contains scripts used to generate UMAPs highlighted by cell types

# Level of Support
We are not currently supporting this code, but simply releasing it to the community AS IS but are not able to provide any guarantees of support. The community is welcome to submit issues, but you should not expect an active response.

# License
The license for this repo is available on Github at: https://github.com/AllenInstitute/CTKE_viz/blob/main/LICENSE.txt
