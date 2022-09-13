#' Helper functions required for main cluster metrics plotting function "make_cluster_metrics.R". Contains functions to 
#' make a single violin plot, with options for highlighted types and no highlighting. I almost never run these in isolation but
#' documenting in case we need to for some reason. 
#' 
#' @param data - a data frame containing taxonomy metadata that is standard for use in Shiny, containing "genes_detected_label" 
#'               columns and "unique.counts" columns
#'               
#' @param filter_level - level of taxonomy to filter by (e.g. if I'm plotting types, I want to filter by subclass_label)
#' 
#' @param X - level of taxonomy to plot individual taxons at (e.g. if I'm plotting subclasses, x = subclass_label)
#' 
#' @param Y - cluster metric. Either "genes_detected_label" or "unique.counts"
#' 
#' @param level_id - data frame containing label, id and color of items to be plotted (this is constructed from the taxonomy file
#'                   provided as part of the main function)
#'                   
#' @param outname - character string specifying name of output file
#' 
#' @param nodename - name of taxon to be highlighted, e.g. "L5 ET"
#' 
#' 
#' 
#' NOTE: These last 3 params are redundant given level_id, need to come back and fix this but it's fine for now
#' 
#' @param setnames - character vector specifying names of taxons to be plotted (this is constructed from the taxonomy file
#'                   provided as part of the main function)
#'                   
#' @param setids - character vector specifying ids of taxons to be plotted (this is constructed from the taxonomy file
#'                   provided as part of the main function)
#'                   
#' @param setcolors - character vector specifying colors of taxons to be plotted (this is constructed from the taxonomy file
#'                   provided as part of the main function)



plot_cluster_metrics_nohighlight <- function( data, 
                                              filter_level, 
                                              filterby,
                                              X,
                                              Y,
                                              level_id,
                                              outname,
                                              nodename,
                                              setnames,
                                              setids,
                                              setcolors) {
  library(svglite)
  library(feather)
  library(tidyverse)
  library(ggplot2)
  library(forcats)
  library(ggnewscale)
  
  # Pull relevant arguments to feed into ggplot
  args <- as.list(sys.call())
  labtest <- as.character(args[[length(args)-6]])
  
  
  # Subset data per user input
  filter_level <- enquo(filter_level)
  X <- enquo(X)
  Y <- enquo(Y)
  level_id <- enquo(level_id)
  sub.dat <- data %>%
    filter(!!filter_level %in% filterby)
  
  # Pull in all labels and ids from taxonomy to ensure nothing is missing in the plot
  sub.dat <- sub.dat %>%
    mutate(X = as_factor(!!X)) 
  
  sub.dat <- sub.dat %>%
    arrange(!!level_id)
  
  sub.dat <- sub.dat %>%
    mutate(level_id = as.numeric(!!level_id)) %>%
    mutate(level_id = cut(!!level_id, breaks = length(setnames)))
  
  
  levels(sub.dat$level_id) <- setids
  catorder <- setnames
  sub.dat$X <- factor(sub.dat$X, catorder)
  
  # Set title of y-axis
  if (labtest == "genes_detected_label") {
    ytitle <- "Number of Genes Detected"
  }
  else {
    ytitle <- "Number of UMIs"
  }
  
  # set which taxon to highlight in the plot
  #lineloc <- which(levels(as.factor(sub.dat$X)) == nodename)
  
  
  
  plot <- sub.dat %>%
    ggplot(aes(x = X,
               y = !!Y,
               fill = X)) +
    
    geom_violin(scale = "width", width = 1) +
    scale_fill_manual(values = setcolors, drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    #geom_vline(aes(xintercept = lineloc-0.5), alpha = 0.2) +
    #geom_vline(aes(xintercept = lineloc+0.5), alpha = 0.2) +
    
    ylab(ytitle) + xlab("") +
    
    theme_light() +
    
    theme(
      aspect.ratio = .75,
      legend.position = "none",
      axis.text.x = element_text(angle =  45, hjust = 1, size = 14),
      axis.text.y = element_text(size = 14),
      strip.background.x = element_blank(),
      strip.background.y = element_rect(fill = "#2c3e50"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  print(plot)
  
  # Write SVG and CSV files to local directory
  nodename <- gsub("/","-",nodename)
  ggsave(filename = paste0(outname,nodename,".svg"))
  #write.table(sub.dat, paste0(outname, nodename, "_data.csv"), sep=",", 
  #            row.names = F)
  
}




plot_cluster_metrics <- function( data, 
                                  filter_level, 
                                  filterby,
                                  X,
                                  Y,
                                  level_id,
                                  outname,
                                  nodename,
                                  setnames,
                                  setids,
                                  setcolors) {
  library(svglite)
  library(feather)
  library(tidyverse)
  library(ggplot2)
  library(forcats)
  library(ggnewscale)
  
  # Pull relevant arguments to feed into ggplot
  args <- as.list(sys.call())
  labtest <- as.character(args[[length(args)-6]])
  
  
  # Subset data per user input
  filter_level <- enquo(filter_level)
  X <- enquo(X)
  Y <- enquo(Y)
  level_id <- enquo(level_id)
  sub.dat <- data %>%
    filter(!!filter_level %in% filterby)
  
  # Pull in all labels and ids from taxonomy to ensure nothing is missing in the plot
  sub.dat <- sub.dat %>%
    mutate(X = as_factor(!!X)) 
  
  sub.dat <- sub.dat %>%
    arrange(!!level_id)
  
  sub.dat <- sub.dat %>%
    mutate(level_id = as.numeric(!!level_id)) %>%
    mutate(level_id = cut(!!level_id, breaks = length(setnames)))
  
  
  levels(sub.dat$level_id) <- setids
  catorder <- setnames
  sub.dat$X <- factor(sub.dat$X, catorder)
  
  # Set title of y-axis
  if (labtest == "genes_detected_label") {
    ytitle <- "Number of Genes Detected"
  }
  else {
    ytitle <- "Number of UMIs"
  }
  
  # set which taxon to highlight in the plot
  lineloc <- which(levels(as.factor(sub.dat$X)) == nodename)
  
  
  
  plot <- sub.dat %>%
    ggplot(aes(x = X,
               y = !!Y,
               fill = X)) +
    geom_rect(data = sub.dat,
              xmin = lineloc-0.5,
              xmax = lineloc+0.5,
              ymin = 0,
              ymax = Inf,
              fill = "grey91",
              alpha=0.18
    ) +
    
    geom_violin(scale = "width", width = 1) +
    scale_fill_manual(values = setcolors, drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    ylab(ytitle) + xlab("") +
    
    theme_light() +
    
    theme(
      aspect.ratio = .75,
      legend.position = "none",
      axis.text.x = element_text(angle =  45, hjust = 1, size = 14),
      axis.text.y = element_text(size = 14),
      strip.background.x = element_blank(),
      strip.background.y = element_rect(fill = "#2c3e50"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  print(plot)
  
  # Write SVG and CSV files to local directory
  nodename <- gsub("/","-",nodename)
  ggsave(filename = paste0(outname,nodename,".svg"))
  #write.table(sub.dat, paste0(outname, nodename, "_data.csv"), sep=",", 
  #            row.names = F)
  
}