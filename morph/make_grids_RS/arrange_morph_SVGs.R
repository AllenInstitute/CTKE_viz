#' Generalized function for arranging morphology reconstructions. Currently it orders them by taxon ID, and so requires the 
#' taxonomy file to be included. Although the function is written to arrange the morph SVGs specifically, it could be adapted
#' to order any image by taxon ID, OR to order by any metric included in the "morph_meta" file. 
#' 
#' @param file_summary - data frame describing each SVG file in the directory to be combined. This is provided by Alice and
#'                       includes the following columns "Image Info", "Image Filename", "Input file" (.swc file used to create
#'                       the reconstruction), and "Input cell/sample id"
#'
#' @param img_placement - data frame that is a subset of standard "Panel_Image Info.csv". Provides placement schema for files to 
#'                        be included in a reconstruction grid
#'                        
#' @param morph_meta - metadata file containing column titled "cell id". Can contain any metdata
#' 
#' @param taxonomy_file - data frame of standard taxonomy file
#' 
#' @param indir - character string of directory containing all SVGs to be looped through
#' 
#' @param outdir - character string of output directory

arrange_morph_SVGs <- function(file_summary,
                               img_placement,
                               moprh_meta,
                               taxonomy_file,
                               indir,
                               outdir) {
  
  suppressPackageStartupMessages({
    library(svglite)
    library(gridExtra)
    library(grImport2)
    library(rsvg)
    library(gridSVG)
    library(ggplot2)
    library(cowplot)
    library(tidyverse)
    library(ggpubr)}) 

  taxonomy_file <- cldf
  
  # Subset morph_meta to include only cells that are included in file_summary
  morph_meta <- morph_meta[morph_meta$`cell id` %in% file_summary$`Input cell/sample id`,]

  # Clean image placement file
  img_placement <- img_placement[grepl("Example Reconstructions",img_placement$Panel_Title),]


  # Helper function to add cell type name to reconstruction
  label_helper <- function(x,label,size){
    new <- pictureGrob(x,ext = "gridSVG", delayContent = FALSE)
    text_lab <- text_grob(label = label, size = size, family = "Arial", face = "bold")
    newg <- grid.arrange(top=text_lab,new)
    newg
  }

  
  # split PI csv into Siblings and All Reconstruction panels and loop through accession IDs
  IP_allrecons <- img_placement[img_placement$Panel_Title == "Example Reconstructions - All Reconstructions",]
  IP_allrecons$taxon_name <- sub(" exemplar reconstruction.*", "", IP_allrecons$Panel_Image)
  

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
    test_PI$cluster_id <- cldf$cluster_id[match(test_PI$cluster_label,cldf$cluster_label)]
    test_PI <- test_PI[with(test_PI, order(test_PI$cluster_id)),]
    test_PI$Order <- seq(1,nrow(test_PI),1)
  
  
    # Read SVG files in order
    setwd(indir)
    all_SVG <- list()
    for (i in 1:nrow(test_PI)) {
      rsvg_svg(test_PI$Image_File_Name[i], paste0(test_PI$Image_File_Name[i],"-cairo.svg"))
      temp <- readPicture(paste0(test_PI$Image_File_Name[i],"-cairo.svg"))
      all_SVG <- c(all_SVG, temp)
    }
  
    # Add labels to reconstructions and plot in taxonomy order
    setwd(outdir)
  
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
  
  
    ggsave(paste0(taxon,"_All_Reconstructions_Grid.svg"),g,width = 16, height = 10)
  
  }



  #####

  ## Repeat with sibling reconstructions
  IP_sibs <- img_placement[img_placement$Panel_Title == "Example Reconstructions - Compare to Sibling Types",]
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
    test_PI$cluster_id <- cldf$cluster_id[match(test_PI$cluster_label,cldf$cluster_label)]
    test_PI <- test_PI[with(test_PI, order(test_PI$cluster_id)),]
    test_PI$Order <- seq(1,nrow(test_PI),1)
  
  
    # Read SVG files in order
    setwd(indir)
    all_SVG <- list()
    for (i in 1:nrow(test_PI)) {
      rsvg_svg(test_PI$Image_File_Name[i], paste0(test_PI$Image_File_Name[i],"-cairo.svg"))
      temp <- readPicture(paste0(test_PI$Image_File_Name[i],"-cairo.svg"))
      all_SVG <- c(all_SVG, temp)
    }
  
  
    # Add labels to reconstructions and plot in taxonomy order
    setwd(outdir)
  
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

}



