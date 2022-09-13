##
# This is a modified version of the "group_dot_plot" function from scrattch.vis that I used to make dotplots
# for cell type cards. It takes a few more arguments, but creates placeholders for taxons that may be missing
# from the dataset to keep things visually consistent, and also adds light grey bars around a given taxon
# to highlight it for the card it's on. The same modifications could easily be made to the group_violin_plot
# function from scrattch.vis to get similar plots. 

# The wrapper "plot_all_dots" was written to streamline plotting multiple taxons at a time,
# but honestly could probably be consolidated into the modified_dot3 function. 

# Arguments for modified_dot3
# data - a count matrix with rownames as genes and colnames as sample ids
# anno - an annotation matrix
# genes - list of ordered genes to plot
# grouping - cluster, subclass, or class
# parentlevel - next step up in hierarchy
# setnames - names of taxonomy siblings
# setids - ids of taxonomy siblings
# nodename - the desired file name for your output image (I usually just make it the taxon I'm plotting for a card), and which 
#            taxon to add the highlighted grey bar to
# group_order - for custom ordering. Leave this NULL
# keep.seq - used to get around an issue I had with the marmoset taxonomy ordering siblings incorrectly. Leaving as FALSE is fine
# fill_stat - what the dot fill color represents
# size_stat - same as above for size
# max_size - max size of dots
# log_scale - log scale the data
# normalize_rows - false by default
# colorset - specify color gradients of dot. I use white -> red by default
# font_size
# label_height - made short by default via suggestion from Brian S.
# show_counts - show number in cluster at top of plot
# rotate_counts = rotate the number in cluster
# max_width = not sure if I use this for anything
# return_type = "plot"


library(scrattch.vis)


# Main dot plot function
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
                          fill_stat = "median",
                          size_stat = "prop_gt0",
                          max_size = 22,
                          log_scale = TRUE,
                          normalize_rows = FALSE,
                          colorset = c("#fee8c8","#fdbb84","#e34a33","orangered"),
                          font_size = 8, 
                          label_height = 2,
                          show_counts = TRUE, 
                          rotate_counts = FALSE,
                          max_width = 14,
                          return_type = "plot") {
  # Reverse so that genes go from top to bottom
  # instead of bottom to top.
  genes <- rev(genes)
  gene_labels <- genes
  
  
  
  
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
  
  # Add x-positions for each group
  plot_data <- add_group_xpos(plot_data,
                              group_cols = group_cols,
                              group_order = group_order)
  
  
  n_stats <- get_n_stats(plot_data, group_cols, genes)
  
  missing <- setnames[which(!(setnames %in% unique(plot_data[[2]])))]
  
  if (length(missing) > 0) {
    for (i in 1:length(missing)) {
      missing_id <- setids[[2]][setids[[1]] %in% missing[i]]
      missing_label <- missing[i]
      missing_color <- setids[[3]][setids[[1]] %in% missing[i]]
      plot_data <- rbind(plot_data, c(missing_id,as.character(missing_label),as.character(missing_color), rep("#FFFFFF",length(genes)), rep(0,length(genes) + 1), missing_id))
    }
  }
  
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
  
  # Plot setup
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
    theme(axis.text = element_text(size = rel(1), face = "italic",family = "Helvetica"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_text(vjust=-3.75),
          panel.border = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          text = element_text(family = "Helvetica"),
          legend.position = "bottom", legend.title = element_blank()) +
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
  
  p <- p + scale_color_manual(name = "Expression Level", 
                              breaks = c(colorset[1], colorset[4]), 
                              values = colorset,
                              labels = c("Low", "High"))
  
  
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
  
  
  
  
  y_label <- max(header_labels$ymax) - 0.1 * label_y_size
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
                         color = "#696969") +
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
                             color = "#696969"),
                         size = pt2mm(font_size-1)) +
        
        geom_fit_text(data = header_labels, 
                  aes(x = (xmin + xmax)/2,
                      y = y_label+1.1, label = label, 
                      family = "bold"), 
                  hjust = 0, vjust = 0.35,reflow = TRUE)
    } 
  }
  
  
}







#######

# Wrapper for dot plot function to plot taxonomy nodes that may or may not actually be useful. Basically
# allows for filtering (filter_level) of the data and to add the name of a dataset to the output file when you 
# are plotting multiple taxonomy nodes at once.
plot_all_dots <- function(data,
                          anno,
                          genes,
                          grouping,
                          colorset,
                          filter_by,
                          nodename,
                          setids,
                          dataname,
                          setnames,
                          keep.seq,
                          font_size = 12) {
  
  
  if (keep.seq == TRUE) {
    plt <- modified_dot3(data = data, 
                         anno = anno, 
                         genes = genes, 
                         grouping = grouping, 
                         parentlevel = filter_by,
                         setnames = setnames,
                         setids = setids,
                         nodename = nodename,
                         colorset = colorset,
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
                         log_scale = TRUE,
                         font_size = font_size,
                         max_width = 14,
                         rotate_counts = TRUE,
                         keep.seq = FALSE)
  }
  
  print(plt)
  
  # Write SVG and CSV files to local directory
  nodename <- gsub("/","-",nodename)
  ggsave(filename = paste0("gene expression_",dataname,nodename,".svg"),width=12.3, height=9.8, units="in")
  #alldat <- cbind.data.frame(anno,data[data$sample_name %in% anno$sample_name,])
  #alldat <- alldat[!duplicated(as.list(alldat))]
  #write.table(alldat, paste0("gene expression_",dataname, nodename, "_data.csv"), sep=",", 
  #            row.names = F)
  
}


######################-- EXAMPLE --######################
### Run "plot_all_dots" in a loop to make dot plots for all glutamatergic subclasses in the
### marmoset taxonomy. Data and real script at: 
### \\allen\programs\celltypes\workgroups\humancelltypes\RayS\CTCard_Content\marmoset

for (i in 1:length(glutsubnames)) {
  
  plot_all_dots(data = plotcounts,
                anno = tenxmeta,
                genes = glutsubgenes,
                grouping = "subclass",
                filter_by = "Glutamatergic",
                nodename = glutsubnames[i],
                dataname = "marmoset10X_",
                setnames = glutsubnames,
                font_size = 12,
                header_text = 12)
  
  
}

