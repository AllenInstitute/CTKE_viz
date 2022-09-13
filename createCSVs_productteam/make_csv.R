# Auto-generate CSVs with critical metadata needed for the product team to populate cell type cards

# This script provides a function to create the main Panel_Image Info CSV requested by the product team for each panel. 
# The columns of the final output CSV include: Panel_Image (a text description of the image), Panel_Title 
# (in this case, "Gene Expression"), Taxon_Accession_ID, Dataset_ID (in this case the Broad dataset), Image_File_Name, 
# and Order (the order each image should appear, from left to right or top to bottom if it's for one of the dropdown menus).

# This approach was used for all datasets for all panels:
# - 1. Load Accession IDs and metadata
# - 2. Clean names as necessary to match files
# - 3. Read in files and create vectors of the easy stuff (Dataset_ID for example)
# - 4. Match Accession IDs and cbind everything to a dataframe
# - 5. Write dataframe to CSV
# - 6. Rinse & repat


#' @param taxonomy_file - standard taxonomy file as a data frame
#' @param nomenclature_file - data frame of standard nomenclature file, that includes an added column titled "Cell Type Card?" with Yes or No
#'                        This is added to filter out the stuff not going on cards (all the aliases that also have accession IDs)
#' @param panel_title - Could be Gene Expression, Cluster Metrics, etc.
#' @param dataset_name - Could be "Human10X", "MarmosetSSv4", etc.
#' @param image_order - order that image should appear on the page 
#' @param filedir - directory with image files to build CSV from
#' @param outdir - directory to save CSV to
#' @param desc_str - a string describing the images, for example: "Dot plot of marker gene expression, Human SSv4"
#' @param outname - name for output file

# CSVs can then be combined using a simple script (load and cbind or use something like vroom() )

make_csv <- function(taxonomy_file,
                     nomenclature_file,
                     panel_title,
                     dataset_name,
                     image_order,
                     filedir,
                     outdir,
                     desc_str,
                     outname){
  suppressPackageStartupMessages({
    library(tidyverse)
    library(tidyr)
    library(dplyr)
    library(plyr)
  })

  clupdate <- taxonomy_file
  nomenclature <- nomenclature_file
  nomenclature <- nomenclature %>% select(`cell set preferred alias`,`cell set accession`,`Cell Type Card?`)
  nomenclature <- nomenclature[nomenclature$`Cell Type Card?` != "No",]
  colnames(nomenclature) <- c("alias","accession","card")

  nomtypes <- nomenclature %>% filter(str_detect(card, "Cell Type"))
  nomother <- nomenclature[nomenclature$alias != nomtypes$alias,]

  dir <- filedir
  filenames <- list.files(path = dir, pattern = ".svg")
  filenames_cv <- filenames[grep("_child_view",filenames)]
  filenames_default <- setdiff(filenames,filenames_cv)


  paneltitles <- rep(panel_title, length(filenames_default))
  dataset <- rep(dataset_name, length(filenames_default))
  order <- rep(image_order,length(filenames_default))


  taxons <- sub(paste0(".*",dataset_name,"_"),"",filenames)
  taxons <- gsub(".svg","",taxons)
  taxons[taxons == "L2-3 IT"] <- "L2/3 IT"
  taxons[taxons == "L2-3 IT_1"] <- "L2/3 IT_1"
  taxons[taxons == "L2-3 IT_2"] <- "L2/3 IT_2"
  taxons[taxons == "L2-3 IT_3"] <- "L2/3 IT_3"
  taxons[taxons == "L4-5 IT"] <- "L4/5 IT"
  taxons[taxons == "L4-5 IT_1"] <- "L4/5 IT_1"
  taxons[taxons == "L4-5 IT_2"] <- "L4/5 IT_2"
  taxons[taxons == "L5-6 NP"] <- "L5/6 NP"
  taxons[taxons == "L5-6 NP_1"] <- "L5/6 NP_1"
  taxons[taxons == "L5-6 NP_2"] <- "L5/6 NP_2"
  taxons[taxons == "L5-6 NP_3"] <- "L5/6 NP_3"
  taxons[taxons == "L5-6 NP CT"] <- "L5/6 NP CT"
  taxons <- unique(taxons)

  overlapnom <- nomenclature[nomenclature$alias %in% taxons,]
  overlapnom <- overlapnom[order(overlapnom$alias),]


  descriptions <- paste0(desc_str,"_", taxons)

  CSV <- cbind.data.frame(descriptions, paneltitles, overlapnom$accession, dataset, filenames, order)
  colnames(CSV) <- c("Panel_Image","Panel_Title","Taxon_Accession_ID","Dataset_ID","Image_File_Name","Order")

  write.table(CSV,paste0(outdir,"/",outname,".csv"),sep=",",row.names = F)
}
