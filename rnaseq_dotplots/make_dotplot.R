##
#' This is a modified version of the "group_dot_plot" function from scrattch.vis that I used to make dotplots
#' for cell type cards. It takes a few more arguments, but creates placeholders for taxons that may be missing
#' from the dataset to keep things visually consistent, and also adds light grey bars around a given taxon
#' to highlight it for the card it's on. The same modifications could easily be made to the group_violin_plot
#' function from scrattch.vis to get similar plots. 
#' 
#' There are two functions included: 
#'  - modified_dot3, which does the actual plotting
#'  - plot_all_dots, which is a wrapper for modified_dot3. Allows one to select markers to highlight for a given type
#' 
#' 
#' 
#' @param data - count matrix with sample names
#' @param anno - metadata data frame
#' @param genes - list of genes to plot
#' @param grouping - level of taxonomy to group at (e.g. "subclass")
#' @param parentlevel - parent level of taxon to be plotted (e.g. "class")
#' @param setnames - character vector of names to compare featured taxon against
#' @param setids - vector of IDs to compare featured taxon against
#' @param nodename - featured taxon (e.g. "L5 ET" if making a plot for the L5 ET cell type card)
#' 
#' ## plot_all_dots args
#' @param filter_by - "parentlevel" argument for modified_dot3 (should probably change the name of this param)
#' @param dataname - character string specifiying name of dataset, to be included in output file name
#' @param markers - character vector of markers to be highlighted directly on the plot. Must be a subset of "genes" from
#'                  modified_dot3
#' 

modified_dot3 <- function(data,
                          anno,
                          genes,
                          grouping,
                          parentlevel,
                          setnames,
                          setids,
                          nodename,
                          group_order = NULL,
                          keep.seq = FALSE,
                          fill_stat = "tmean",
                          size_stat = "prop_gt0",
                          max_size = 48,
                          markers = NULL,
                          log_scale = TRUE,
                          normalize_rows = FALSE,
                          colorset = c("#fee8c8","#fdbb84","#f58d71","#e34a33","orangered","#dd1f13"),
                          font_size = 15, 
                          label_height = 2,
                          show_counts = TRUE, 
                          rotate_counts = FALSE,
                          max_width = 14,
                          return_type = "plot") {
  require(plyr)
  require(tidyverse)
  require(scrattch.hicat) # might not be required
  require(scrattch.vis)
  require(ggfittext)
  
  # Reverse so that genes go from top to bottom
  # instead of bottom to top.
  genes <- rev(genes)
  gene_labels <- genes
  gene_colors <- rep("black",length(gene_labels))
  
  
  if(!is.null(markers)) {
    gene_inds <- which(gene_labels %in% markers)
    gene_colors[gene_inds] <- "#219653"
    gene_labels[gene_inds] <- paste0("\u25b6"," ",gene_labels[gene_inds])
  }
  
  
  group_cols <- group_columns(grouping)
  
  # Filter data to genes and samples in anno
  gene_data <- filter_gene_data(data, 
                                genes, 
                                anno, 
                                group_cols,
                                group_order, 
                                "sample_name")
  
  # Filter annotations if group_order is provided
  if(!is.null(group_order)) {
    anno <- anno[anno[[group_cols$id]] %in% group_order,]
  }
  
  gene_fill_stats <- group_stats(gene_data,
                                 value_cols = genes,
                                 anno = anno,
                                 grouping = group_cols$label,
                                 stat = fill_stat)
  
  # Get maximum values for each gene before rescaling to plot space.
  max_vals_unscaled <- max_gene_vals(gene_fill_stats, genes)
  
  gene_size_stats <- group_stats(gene_data,
                                 value_cols = genes,
                                 anno = anno,
                                 grouping = group_cols$label,
                                 stat = size_stat)
  
  if(log_scale) {
    gene_fill_stats <- scale_gene_data(gene_fill_stats, genes, scale_type = "log10")
  }
  
  # Convert the data values to heatmap colors
  gene_fill_data <- data_df_to_colors(gene_fill_stats,
                                      value_cols = genes,
                                      per_col = normalize_rows,
                                      colorset = colorset)
  
  names(gene_fill_data)[match(genes, names(gene_fill_data))] <- paste0(genes, "_fill")
  names(gene_size_stats)[match(genes, names(gene_size_stats))] <- paste0(genes, "_size")
  
 
  
  # Left-join data to anno. This will ensure that data is filtered for the cells provided in anno
  plot_anno <- anno %>%
    select(one_of(group_cols$id, group_cols$label, group_cols$color)) %>%
    unique()
  
  group_counts <- anno %>%
    group_by_(group_cols$id) %>%
    summarise(group_n = n())
  
  plot_data <- plot_anno %>%
    left_join(gene_fill_data, by = group_cols$label) %>%
    left_join(gene_size_stats, by = group_cols$label) %>%
    left_join(group_counts, by = group_cols$id)
  
  # Create a placeholder for cell types that might be missing from dataset but are still present in the taxonomy
  missing <- setnames[which(!(setnames %in% unique(plot_data[[2]])))]
  if (length(missing) > 0) {
    for (i in 1:length(missing)) {
      missing_id <- setids[[2]][setids[[1]] %in% missing[i]]
      missing_label <- missing[i]
      missing_color <- setids[[3]][setids[[1]] %in% missing[i]]
      plot_data <- rbind(plot_data, c(missing_id,as.character(missing_label),as.character(missing_color), rep("#FFFFFF",length(genes)), rep(0,length(genes) + 1), missing_id))
    }
  }
  
  # Add x-positions for each group
  plot_data <- add_group_xpos(plot_data,
                              group_cols = group_cols,
                              group_order = group_order)
  
  
  n_stats <- get_n_stats(plot_data, group_cols, genes)
  
  
  firstcols <- colnames(plot_data)[1:3]
  fillcols <- colnames(plot_data)[grepl("fill",colnames(plot_data))]
  charcols <- c(firstcols, fillcols)
  
  plot_data <- plot_data %>% mutate_at(colnames(plot_data)[!(colnames(plot_data) %in% charcols)], as.numeric)
  plot_data <- plot_data[which(unique(plot_data[[2]]) %in% setnames),]
  plot_data <- plot_data[order(plot_data[[1]]),]
  plot_data$xpos <- as.numeric(plot_data[[1]])
  
  if (keep.seq == TRUE) {
    plot_data$xpos <- seq(1,nrow(plot_data),1)
  }
  
  n_stats <- get_n_stats(plot_data, group_cols, genes)
  header_labels <-build_header_labels(data = plot_data, 
                                      grouping = grouping,
                                      group_order = group_order,
                                      ymin = n_stats$genes + 1, 
                                      label_height = label_height, 
                                      label_type = "simple")
  
  label_y_size <- max(header_labels$ymax) - min(header_labels$ymin)
  
  group_data <- plot_data %>%
    select(xpos, group_n) %>%
    mutate(label_y = n_stats$genes + label_y_size * 0.05,
           group_n_y = max(header_labels$ymax) - 0.1 * label_y_size)
  
  colnames(plot_data) <- make.names(colnames(plot_data))
  
  # Scale max dot size based on number of genes to be plotted (for legibility)
  if (length(genes) <= 25) {
    max_size <- 48
  }
  
  if (length(genes) > 25) {
    max_size <- 28
  }
  
  # Plot setup. Place a light grey box where the type of interest is (specified by nodename)
  p <- ggplot() +
    scale_fill_identity(name="Median expression level") +
    scale_size_area(max_size = pt2mm(max_size)) +
    scale_y_continuous("", 
                       breaks = 1:length(gene_labels) + 0.45, 
                       labels = gene_labels, 
                       expand = c(0, 1)) +
    scale_x_discrete(name="Proportion of cells expressing gene", 
                     expand = c(0, 0)) +
    theme_light(font_size) +
    theme(axis.text = element_text(size = rel(1), face = "italic",family = "Helvetica",colour = gene_colors),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_text(vjust=-3.75),
          panel.border = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          text = element_text(family = "Helvetica"),
          legend.position = "bottom", legend.title = element_blank(),
          plot.background= element_rect(color= "transparent")) +
    geom_hline(aes(yintercept = 1:(n_stats$genes)), size = 0.2,alpha = 0.001) 
  lineloc <- plot_data$xpos[plot_data[[2]] == nodename]
  p <- p + geom_rect(data = plot_data,
                     xmin = lineloc-0.5,
                     xmax = lineloc+0.5,
                     ymin = 0,
                     ymax = Inf,
                     fill = "grey91",
                     alpha=0.18
  ) 
  
  # Create dots
  plot_data <- plot_data[order(plot_data$xpos),]
  for(i in 1:length(genes)) {
    gene <- make.names(genes[[i]])
    gene_fill <- paste0(gene, "_fill")
    gene_size <- paste0(gene, "_size")
    p <- p + 
      geom_point(data = plot_data,
                 aes_string(x = "xpos",
                            y = i + 0.5, 
                            fill = gene_fill,
                            size = gene_size,
                            shape = 16,
                            stroke = 0),
                 pch = 21, 
                 show.legend = TRUE) 
    
    
  }
  
  # Create heatmap legend
  colors <- colorRampPalette(colorset)(101)
  legend_data <- data.frame(xmin = 1:101,
                            xmax = 1:101 + 1,
                            ymin = 0,
                            ymax = 0.1,
                            fill = colors)
  
  
  min_val <- round(min(gene_fill_stats[,c(2:ncol(gene_fill_stats))]),2)
  max_val <- round(max(gene_fill_stats[,c(2:ncol(gene_fill_stats))]),2)
  scale_name <- "Log-normalized counts"
  
  legend_plot <- ggplot(legend_data) +
    geom_rect(aes(xmin = xmin, xmax = xmax, 
                  ymin = ymin, ymax = ymax, 
                  fill = fill)) +
    geom_segment(aes(x = min(xmin), xend = max(xmax), 
                     y = 0, yend = 0)) +
    scale_fill_identity() +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(scale_name, 
                       breaks = c(0, 25, 50, 75, 100),
                       labels = round(seq(min_val, max_val, by = (max_val - min_val) / 4), 2)) +
    theme_classic() +
    coord_cartesian(clip="off") +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.line.x = element_blank())
  
  
  ggplot_header_labels2 <- function (p, header_labels, header_polygons = NULL) 
  {
    p <- p + geom_rect(data = header_labels, aes(xmin = xmin+0.01, 
                                                 xmax = xmax-0.01, ymin = ymin, ymax = ymax, fill = color)) 
    
    if (!is.null(header_polygons)) {
      p <- p + geom_polygon(data = header_polygons, aes(x = poly.x, 
                                                        y = poly.y, fill = color, group = id))
    }
    p
  }
  
  p <- ggplot_header_labels2(p,
                             header_labels = header_labels,
                             header_polygons = NULL)
  
  p <- p + theme(legend.position = "bottom", legend.title = element_blank())
  
  
  p <- p + theme(legend.position = "bottom") 
  
  
  
   
  # Add labels at the top of the plot  
  y_label <- max(header_labels$ymax) - 0.1 * label_y_size
  group_data <- group_data[order(group_data$xpos),]
  group_data$xmin <- header_labels$xmin
  group_data$xmax <- header_labels$xmax
  # Cluster counts
  if (show_counts) {
    if (rotate_counts) {
      p <- p + geom_text(data = group_data,
                         aes(x = (xmin + xmax)/2,
                             y = group_n_y+0.4, 
                             label = group_n),
                         angle = 0, 
                         size = pt2mm(font_size-1),
                         color = "black") +
        geom_fit_text(data = header_labels, 
                  aes(x = (xmin + xmax)/2,
                      y = y_label+1.1, label = label, 
                      family = "bold"), 
                  hjust = 0, vjust = 0.35,reflow = TRUE)
    } 
    else {
      p <- p + geom_text(data = group_data,
                         aes(x = (xmin + xmax)/2,
                             y = group_n_y+0.4, 
                             label = group_n,
                             color = "black"),
                         size = pt2mm(font_size-1)) +
        
        geom_fit_text(data = header_labels, 
                  aes(x = (xmin + xmax)/2,
                      y = y_label+1.1, label = label, 
                      family = "bold"), 
                  hjust = 0, vjust = 0.35,reflow = TRUE)
    } 
  }
  
  
  
  # Combine plot and legend and return
  p <- ggarrange(p, legend_plot, widths = c(4,0.2), heights = c(6,0.25),ncol=1,nrow = 2)
  
  
  
}



# Wrapper for modified_dot3 used in scripts to generate all plots for cell type cards
plot_all_dots <- function(data,
                          anno,
                          genes,
                          grouping,
                          colorset = c("#fee8c8","#fdbb84","#f58d71","#e34a33","orangered","#dd1f13"),
                          filter_by,
                          nodename,
                          setids,
                          dataname,
                          setnames,
                          markers,
                          keep.seq = FALSE,
                          font_size = 15) {
  
  if(keep.seq == TRUE) {
    plt <- modified_dot3(data = data, 
                         anno = anno, 
                         genes = genes, 
                         grouping = grouping, 
                         parentlevel = filter_by,
                         setnames = setnames,
                         setids = setids,
                         nodename = nodename,
                         colorset = colorset,
                         markers = markers,
                         log_scale = TRUE,
                         font_size = font_size,
                         max_width = 14,
                         rotate_counts = TRUE,
                         keep.seq = TRUE)
  }
  
  else {

    plt <- modified_dot3(data = data, 
                         anno = anno, 
                         genes = genes, 
                         grouping = grouping, 
                         parentlevel = filter_by,
                         setnames = setnames,
                         setids = setids,
                         nodename = nodename,
                         colorset = colorset,
                         markers = markers,
                         log_scale = TRUE,
                         font_size = font_size,
                         max_width = 14,
                         rotate_counts = TRUE)
  }
  
 
  
  # Write SVG and CSV files to local directory
  nodename <- gsub("/","-",nodename)
  svg(paste0("gene expression_",dataname,nodename,".svg"),width=12.8, height=19)
  print(plt)
  dev.off()
}

