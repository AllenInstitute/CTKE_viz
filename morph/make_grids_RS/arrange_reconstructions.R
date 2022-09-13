# Trying to combine multiple SVGs into a nice grid

library(svglite)
library(gridExtra)
library(grImport2)
library(rsvg)
library(gridSVG)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(ggpubr)

# Read in file information
file_summary <- read_csv("//allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/morph/delivery/ctc_summary.csv")
file_summary <- file_summary[1:134,]

# Read in morphology summary data
morph_meta <- read_csv("//allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/morph/Tolias_m1_patchseq_morph_features.csv")
morph_meta <- morph_meta[morph_meta$`cell id` %in% file_summary$`Input cell/sample id`,]

# Read in panel_image info for morph
image_placement <- read_csv("//allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/morph/Panel_Image Info_copy.csv")
image_placement <- image_placement[grepl("Example Reconstructions",image_placement$Panel_Title),]


##  Now split PI csv into Siblings and All Reconstruction panels and loop through accession IDs
label_helper <- function(x,label,size){
  new <- pictureGrob(x,ext = "gridSVG", delayContent = FALSE)
  text_lab <- text_grob(label = label, size = size, family = "Arial", face = "bold")
  newg <- grid.arrange(top=text_lab,new)
  newg
}

IP_allrecons <- image_placement[image_placement$Panel_Title == "Example Reconstructions - All Reconstructions",]
IP_allrecons$taxon_name <- sub(" exemplar reconstruction.*", "", IP_allrecons$Panel_Image)
mouse_cldf <- read_csv("//allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr/cl_updated.df.csv")
setwd("/allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/morph/delivery/image_grids")

alltaxons <- unique(IP_allrecons$taxon_name)
for (i in 1:length(alltaxons)) {

# Subset based on taxon name 
  taxon <- alltaxons[i]
  test_PI <- IP_allrecons[IP_allrecons$taxon_name == taxon,]
  test_summary <- file_summary[match(test_PI$Image_File_Name,file_summary$`Image Filename`),]
  tempmeta <- morph_meta[morph_meta$`cell id` %in% gsub(".svg","",test_PI$Image_File_Name),]
  cellids <- gsub(".svg","",test_PI$Image_File_Name)
  tempmeta <- tempmeta[order(match(tempmeta$`cell id`,cellids)),]
  test_PI$depth <- tempmeta$`normalized depth`
  test_PI$cluster_label <- test_summary$`Image info`[test_summary$`Image Filename` %in% test_PI$Image_File_Name]
  test_PI$cluster_label <- gsub(" exemplar reconstruction","",test_PI$cluster_label)
  test_PI$cluster_id <- mouse_cldf$cluster_id[match(test_PI$cluster_label,mouse_cldf$cluster_label)]
  test_PI <- test_PI[with(test_PI, order(test_PI$cluster_id)),]
  test_PI$Order <- seq(1,nrow(test_PI),1)
  

# Read SVG files in order
  setwd("/allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/morph/delivery/all_images")
  all_SVG <- list()
  for (i in 1:nrow(test_PI)) {
    rsvg_svg(test_PI$Image_File_Name[i], paste0(test_PI$Image_File_Name[i],"-cairo.svg"))
    temp <- readPicture(paste0(test_PI$Image_File_Name[i],"-cairo.svg"))
    all_SVG <- c(all_SVG, temp)
  }

# Add labels to reconstructions and plot in taxonomy order
  setwd("/allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/morph/delivery/image_grids3")
  
  e <- new.env()
  testlabs <- test_PI$cluster_label
  all_SVG <- setNames(all_SVG,testlabs)
  listinds <- which(duplicated(names(all_SVG)))
  names(all_SVG)[listinds] <- paste0(names(all_SVG)[listinds]," ")
  list2env(all_SVG, envir = e)
  
  
  labinds <- which(duplicated(testlabs))
  testlabs[labinds] <- paste0(testlabs[listinds]," ")
  new_e <- list()
  
  if (length(e) > 6){
    
  for (i in 1:length(e)) {
    cellid <- names(e)[order(match(names(e), testlabs))][i]
    temp <- label_helper(e[[cellid]], label = str_wrap(names(e)[order(match(names(e), testlabs))][i], width=10),size = 16)
    new_e[[i]] <- temp
    }
  }
  
  
  if (length(e) <= 6){
    
    for (i in 1:length(e)) {
      cellid <- names(e)[order(match(names(e), testlabs))][i]
      temp <- label_helper(e[[cellid]], label = names(e)[order(match(names(e), testlabs))][i],size = 22)
      new_e[[i]] <- temp
      }  
    }
 
  if (length(new_e) <= 8){
    g <- do.call("arrangeGrob", c(new_e, ncol = length(new_e)))
  } else {
    g <- do.call("arrangeGrob", c(new_e, ncol = 8))
  }
  
  
  taxon <- gsub("/","-",taxon)
  
  
  #svg(paste0(taxon,"_All_Reconstructions_Grid.svg"),width=16, height=10)
  #plot_grid(g)
  #dev.off()
  ggsave(paste0(taxon,"_All_Reconstructions_Grid.svg"),g,width = 16, height = 10)

}



#####

## Repeat with sibling reconstructions
IP_sibs <- image_placement[image_placement$Panel_Title == "Example Reconstructions - Compare to Sibling Types",]
IP_sibs$taxon_name <- sub(" exemplar reconstruction.*", "", IP_sibs$Panel_Image)
IP_sibs$taxon_name <- sub(" sibling", "", IP_sibs$taxon_name)

alltaxons <- unique(IP_sibs$taxon_name)
for (i in 1:length(alltaxons)) {
  
  # Subset based on taxon name 
  taxon <- alltaxons[i]
  test_PI <- IP_sibs[IP_sibs$taxon_name == taxon,]
  test_summary <- file_summary[match(test_PI$Image_File_Name,file_summary$`Image Filename`),]
  tempmeta <- morph_meta[morph_meta$`cell id` %in% gsub(".svg","",test_PI$Image_File_Name),]
  cellids <- gsub(".svg","",test_PI$Image_File_Name)
  tempmeta <- tempmeta[order(match(tempmeta$`cell id`,cellids)),]
  test_PI$depth <- tempmeta$`normalized depth`
  test_PI$cluster_label <- test_summary$`Image info`[test_summary$`Image Filename` %in% test_PI$Image_File_Name]
  test_PI$cluster_label <- gsub(" exemplar reconstruction","",test_PI$cluster_label)
  test_PI$cluster_id <- mouse_cldf$cluster_id[match(test_PI$cluster_label,mouse_cldf$cluster_label)]
  test_PI <- test_PI[with(test_PI, order(test_PI$cluster_id)),]
  test_PI$Order <- seq(1,nrow(test_PI),1)
  
  
  # Read SVG files in order
  setwd("/allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/morph/delivery/all_images")
  all_SVG <- list()
  for (i in 1:nrow(test_PI)) {
    rsvg_svg(test_PI$Image_File_Name[i], paste0(test_PI$Image_File_Name[i],"-cairo.svg"))
    temp <- readPicture(paste0(test_PI$Image_File_Name[i],"-cairo.svg"))
    all_SVG <- c(all_SVG, temp)
  }
  
  
  # Add labels to reconstructions and plot in taxonomy order
  setwd("/allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/morph/delivery/image_grids3")
  
  e <- new.env()
  testlabs <- test_PI$cluster_label
  all_SVG <- setNames(all_SVG,testlabs)
  listinds <- which(duplicated(names(all_SVG)))
  names(all_SVG)[listinds] <- paste0(names(all_SVG)[listinds]," ")
  list2env(all_SVG, envir = e)
  
  new_e <- list()
  
  labinds <- which(duplicated(testlabs))
  testlabs[labinds] <- paste0(testlabs[listinds]," ")
  if (length(e) > 6){
    
    for (i in 1:length(e)) {
      cellid <- names(e)[order(match(names(e), testlabs))][i]
      temp <- label_helper(e[[cellid]], label = str_wrap(names(e)[order(match(names(e), testlabs))][i], width=10),size = 16)
      new_e[[i]] <- temp
    }
  }
  
  
  if (length(e) <= 6){
    
    for (i in 1:length(e)) {
      cellid <- names(e)[order(match(names(e), testlabs))][i]
      temp <- label_helper(e[[cellid]], label = names(e)[order(match(names(e), testlabs))][i],size = 22)
      new_e[[i]] <- temp
    }  
  }
  
  if (length(new_e) <= 8){
    g <- do.call("arrangeGrob", c(new_e, ncol = length(new_e)))
  } else {
    g <- do.call("arrangeGrob", c(new_e, ncol = 8))
  }
  taxon <- gsub("/","-",taxon)
  ggsave(paste0(taxon,"_Sibling_Reconstructions_Grid.svg"),g,width = 16, height = 10)
  
}




