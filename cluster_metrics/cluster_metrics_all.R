# This script is for displaying the 2 different cluster metrics for cell type cards as violin plots, 
# in the same style as Nik Jorstad's original plots in the Google Drive. 
# It will output the plot in SVG format and a CSV file of the metadata selection necessary to reproduce the plot.
# This script plots each dataset one at a time (so is a series )

# - Ray S.
#
# The function below will take the following arguments to produce a plot:
#
#
# data -->         metadata dataframe, probably loaded from .feather, .RDA or .mtx files

# filter_level --> The level of the tree you want to filter at (class_label, subclass_label or cluster_label).
#                  For example, entering subclass_label will allow you to filter out any subclass or subclasses

# filterby -->     Character vector of desired nodes (class, subclasses or clusters) to filter out of
#                  entire metadata set. This should match filter_level. For example, if you enter class_label,
#                  you could enter "Non-Neuronal", "GABAergic", or "Glutamatergic".

# X -->            The level of the tree you want to visualize at (class_label, subclass_label or cluster_label). For
#                  example, if you filtered Glutamatergic class you could choose to show subclass level data only or
#                  all glutamatergic clusters 

# Y -->            Metric of interest (genes_detected_label or total_umi_label)
#
#

# level_id -->     dataframe containing labels, colors and IDs at the level you want to plot
#
#

# outname -->      What you want your filenames to be called
#
#

# nodename -->     What taxon the plot is for. This is used to "highlight" the taxon without changing the plot order
#
#

# setnames, setids, setcolors --> Information for the grouping you want to plot at. This is important for creating
#                                 placeholders in the plot for datasets where certain taxons may be missing



library(scrattch.vis)
library(tidyverse)
library(mgsub)
library(cowplot)
library(Matrix)
library(feather)
options(stringsAsFactors = F)

# Load cleaned metadata
source("load_metadata.R")
######

# Plotting function 
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
    
    geom_violin(scale = "width", width = 1) +
    scale_fill_manual(values = setcolors, drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    geom_vline(aes(xintercept = lineloc-0.5), alpha = 0.2) +
    geom_vline(aes(xintercept = lineloc+0.5), alpha = 0.2) +
    
    ylab(ytitle) + xlab("") +
    
    theme_light() +
    
    theme(
      aspect.ratio = .75,
      legend.position = "none",
      axis.text.x = element_text(angle =  45, hjust = 1),
      strip.background.x = element_blank(),
      strip.background.y = element_rect(fill = "#2c3e50"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  print(plot)
  
  # Write SVG and CSV files to local directory
  nodename <- gsub("/","-",nodename)
  ggsave(filename = paste0(outname,nodename,".svg"))
  write.table(sub.dat, paste0(outname, nodename, "_data.csv"), sep=",", 
              row.names = F)
  
}




# Make plots for all datasets
#######
# snRNA SMART
##########





setwd("/allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr/snRNA SMART/cluster metrics")
# Plot all class node plots 

for (i in 1:length(classnames)) {
  plot_cluster_metrics(data = smartnucleimet, 
                       filter_level = class_label, 
                       filterby = c("GABAergic","Glutamatergic","Non-Neuronal","Non-Neural"),
                       X = class_label, 
                       Y = genes_detected_label, 
                       level_id = class_id,
                       outname = "genes_detected_snRNASMART_",
                       nodename = classnames[i],
                       setnames = classnames,
                       setids = classids$class_id,
                       setcolors = classids$class_color)
}


for (i in 1:length(classnames)) {
  plot_cluster_metrics(data = smartnucleimet, 
                       filter_level = class_label, 
                       filterby = c("GABAergic","Glutamatergic","Non-Neuronal","Non-Neural"),
                       X = class_label, 
                       Y = unique.counts, 
                       level_id = class_id,
                       outname = "UMIs_snRNASMART_",
                       nodename = classnames[i],
                       setnames = classnames,
                       setids = classids$class_id,
                       setcolors = classids$class_color)
}

# Plot all glut subclass node plots

for (i in 1:length(glutsubnames)) {
  plot_cluster_metrics(data = smartnucleimet, 
                       filter_level = class_label, 
                       filterby = "Glutamatergic",
                       X = subclass_label, 
                       Y = genes_detected_label, 
                       level_id = subclass_id,
                       outname = "genes_detected_snRNASMART_",
                       nodename = glutsubnames[i],
                       setnames = glutsubnames,
                       setids = glutsubids$subclass_id,
                       setcolors = glutsubids$subclass_color)
}

for (i in 1:length(glutsubnames)) {
  plot_cluster_metrics(data = smartnucleimet, 
                       filter_level = class_label, 
                       filterby = "Glutamatergic",
                       X = subclass_label, 
                       Y = unique.counts, 
                       level_id = subclass_id,
                       outname = "UMIs_snRNASMART_",
                       nodename = glutsubnames[i],
                       setnames = glutsubnames,
                       setids = glutsubids$subclass_id,
                       setcolors = glutsubids$subclass_color)
}



# Plot all GABA subclass node plots

for (i in 1:length(gabasubnames)) {
  plot_cluster_metrics(data = smartnucleimet, 
                       filter_level = class_label, 
                       filterby = "GABAergic",
                       X = subclass_label, 
                       Y = unique.counts, 
                       level_id = subclass_id,
                       outname = "UMIs_snRNASMART_",
                       nodename = gabasubnames[i],
                       setnames = gabasubnames,
                       setids = gabasubids$subclass_id,
                       setcolors = gabasubids$subclass_color)
}

for (i in 1:length(gabasubnames)) {
  plot_cluster_metrics(data = smartnucleimet, 
                       filter_level = class_label, 
                       filterby = "GABAergic",
                       X = subclass_label, 
                       Y = genes_detected_label, 
                       level_id = subclass_id,
                       outname = "genes_deteected_snRNASMART_",
                       nodename = gabasubnames[i],
                       setnames = gabasubnames,
                       setids = gabasubids$subclass_id,
                       setcolors = gabasubids$subclass_color)
}



# Plot all Non-neural subclass node plots
for (i in 1:length(nonneuralsubnames)) {
  plot_cluster_metrics(data = smartnucleimet, 
                       filter_level = class_label, 
                       filterby = "Non-Neural",
                       X = subclass_label, 
                       Y = unique.counts, 
                       level_id = subclass_id,
                       outname = "UMIs_snRNASMART_",
                       nodename = nonneuralsubnames[i],
                       setnames = nonneuralsubnames,
                       setids = nonneuralsubids$subclass_id,
                       setcolors = nonneuralsubids$subclass_color)
}

for (i in 1:length(nonneuralsubnames)) {
  plot_cluster_metrics(data = smartnucleimet, 
                       filter_level = class_label, 
                       filterby = "Non-Neural",
                       X = subclass_label, 
                       Y = genes_detected_label, 
                       level_id = subclass_id,
                       outname = "genes_detected_label_snRNASMART_",
                       nodename = nonneuralsubnames[i],
                       setnames = nonneuralsubnames,
                       setids = nonneuralsubids$subclass_id,
                       setcolors = nonneuralsubids$subclass_color)
}


# Plot all Non-neuronal subclass node plots

for (i in 1:length(nonneuronalsubnames)) {
  plot_cluster_metrics(data = smartnucleimet, 
                       filter_level = class_label, 
                       filterby = "Non-Neuronal",
                       X = subclass_label, 
                       Y = genes_detected_label, 
                       level_id = subclass_id,
                       outname = "genes_detected_label_snRNASMART_",
                       nodename = nonneuronalsubnames[i],
                       setnames = nonneuronalsubnames,
                       setids = nonneuronalsubids$subclass_id,
                       setcolors = nonneuronalsubids$subclass_color)
}

for (i in 1:length(nonneuronalsubnames)) {
  plot_cluster_metrics(data = smartnucleimet, 
                       filter_level = class_label, 
                       filterby = "Non-Neuronal",
                       X = subclass_label, 
                       Y = unique.counts, 
                       level_id = subclass_id,
                       outname = "UMIs_snRNASMART_",
                       nodename = nonneuronalsubnames[i],
                       setnames = nonneuronalsubnames,
                       setids = nonneuronalsubids$subclass_id,
                       setcolors = nonneuronalsubids$subclass_color)
}





# Plot types yikes
setwd("/allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr/snRNA SMART/cluster metrics/test")


missingsubs <- setdiff(as.character(subclassnames), unique(as.character(smartnucleimet$subclass_label)))
missingnodes <- setdiff(as.character(typenames), unique(as.character(smartnucleimet$cluster_label)))

snSMARTsub <- setdiff(as.character(subclassnames),missingsubs)

# UMIs
for (i in 1:length(snSMARTsub)) {
  temp <- smartnucleimet %>% filter(subclass_label %in% snSMARTsub[i])
  temp <- temp[order(temp$cluster_id),]
  checktypes <- setdiff(clupdate$cluster_label[clupdate$subclass_label == snSMARTsub[i]], unique(as.character(temp$cluster_label)))
  subsize <- length(unique(temp$cluster_label[temp$subclass_label == snSMARTsub[i]]))
  
  if(length(checktypes) == 0 && subsize > 1) {
    temptypes <- as.character(unique(temp$cluster_label))
    for (j in 1:length(temptypes)) {
      ttnames <- temptypes
      ttids <- clupdate %>% filter(cluster_label %in% temptypes)
      plot_cluster_metrics(data = smartnucleimet, 
                           filter_level = subclass_label, 
                           filterby = snSMARTsub[i],
                           X = cluster_label, 
                           Y = unique.counts, 
                           level_id = cluster_id,
                           outname = "UMIs_snRNASMART_",
                           nodename = ttnames[j],
                           setnames = ttnames,
                           setids = ttids$cluster_id,
                           setcolors = ttids$cluster_color)
      
    }
  }
  
  if (subsize == 1 & length(checktypes) == 0) {
    temptypes <- as.character(unique(temp$cluster_label))
    whichparent <- unique(clupdate$class_label[clupdate$subclass_label == snSMARTsub[i]])
    whichsubs <- unique(clupdate$subclass_label[clupdate$class_label == whichparent])
    whichids <- unique(clupdate$subclass_id[clupdate$class_label == whichparent])
    whichcols <- unique(clupdate$subclass_color[clupdate$class_label == whichparent])
    plot_cluster_metrics(data = smartnucleimet, 
                         filter_level = class_label, 
                         filterby = whichparent,
                         X = subclass_label, 
                         Y = unique.counts, 
                         level_id = subclass_id,
                         outname = "UMIs_snRNASMART_TYPE_",
                         nodename = snSMARTsub[i],
                         setnames = whichsubs,
                         setids = whichids,
                         setcolors = whichcols)
  }
  
  
  
  if (length(checktypes > 0)) {
    temptypes <- c(as.character(unique(temp$cluster_label)), checktypes)
    whichids <- clupdate$cluster_id[clupdate$cluster_label %in% temptypes]
    whichtypes <- clupdate$cluster_label[clupdate$cluster_label %in% temptypes]
    whichcols <- clupdate$cluster_color[clupdate$cluster_label %in% temptypes]
    
    for (j in 1:length(temptypes)) {
      ttnames <- whichtypes
      ttids <- clupdate %>% filter(cluster_label %in% temptypes)
      plot_cluster_metrics(data = smartnucleimet, 
                           filter_level = subclass_label, 
                           filterby = snSMARTsub[i],
                           X = cluster_label, 
                           Y = unique.counts, 
                           level_id = cluster_id,
                           outname = "UMIs_snRNASMART_",
                           nodename = ttnames[j],
                           setnames = ttnames,
                           setids = whichids,
                           setcolors = whichcols)
      
    }
  }
  
  
}

# Genes detected
for (i in 1:length(snSMARTsub)) {
  temp <- smartnucleimet %>% filter(subclass_label %in% snSMARTsub[i])
  temp <- temp[order(temp$cluster_id),]
  checktypes <- setdiff(clupdate$cluster_label[clupdate$subclass_label == snSMARTsub[i]], unique(as.character(temp$cluster_label)))
  subsize <- length(unique(temp$cluster_label[temp$subclass_label == snSMARTsub[i]]))
  
  if(length(checktypes) == 0 && subsize > 1) {
    temptypes <- as.character(unique(temp$cluster_label))
    for (j in 1:length(temptypes)) {
      ttnames <- temptypes
      ttids <- clupdate %>% filter(cluster_label %in% temptypes)
      plot_cluster_metrics(data = smartnucleimet, 
                           filter_level = subclass_label, 
                           filterby = snSMARTsub[i],
                           X = cluster_label, 
                           Y = genes_detected_label, 
                           level_id = cluster_id,
                           outname = "genes_detected_snRNASMART_",
                           nodename = ttnames[j],
                           setnames = ttnames,
                           setids = ttids$cluster_id,
                           setcolors = ttids$cluster_color)
      
    }
  }
  
  if (subsize == 1 & length(checktypes) == 0) {
    temptypes <- as.character(unique(temp$cluster_label))
    whichparent <- unique(clupdate$class_label[clupdate$subclass_label == snSMARTsub[i]])
    whichsubs <- unique(clupdate$subclass_label[clupdate$class_label == whichparent])
    whichids <- unique(clupdate$subclass_id[clupdate$class_label == whichparent])
    whichcols <- unique(clupdate$subclass_color[clupdate$class_label == whichparent])
    plot_cluster_metrics(data = smartnucleimet, 
                         filter_level = class_label, 
                         filterby = whichparent,
                         X = subclass_label, 
                         Y = genes_detected_label, 
                         level_id = subclass_id,
                         outname = "genes_detected_snRNASMART_TYPE_",
                         nodename = snSMARTsub[i],
                         setnames = whichsubs,
                         setids = whichids,
                         setcolors = whichcols)
  }
  
  
  
  if (length(checktypes > 0)) {
    temptypes <- c(as.character(unique(temp$cluster_label)), checktypes)
    whichids <- clupdate$cluster_id[clupdate$cluster_label %in% temptypes]
    whichtypes <- clupdate$cluster_label[clupdate$cluster_label %in% temptypes]
    whichcols <- clupdate$cluster_color[clupdate$cluster_label %in% temptypes]
    
    for (j in 1:length(temptypes)) {
      ttnames <- whichtypes
      ttids <- clupdate %>% filter(cluster_label %in% temptypes)
      plot_cluster_metrics(data = smartnucleimet, 
                           filter_level = subclass_label, 
                           filterby = snSMARTsub[i],
                           X = cluster_label, 
                           Y = genes_detected_label, 
                           level_id = cluster_id,
                           outname = "genes_detected_snRNASMART_",
                           nodename = ttnames[j],
                           setnames = ttnames,
                           setids = whichids,
                           setcolors = whichcols)
      
    }
  }
  
  
}


write.table(c(missingsubs,missingnodes),"snRNA SMART Missing Taxons.txt")


######
# Broad data
#########


setwd("/allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr/snRNA 10xV3 Broad/cluster metrics")
# Plot all class node plots 

for (i in 1:length(classnames)) {
  plot_cluster_metrics(snRNA_10xv3_B_met, class_label, c("GABAergic","Glutamatergic","Non-Neuronal","Non-Neural"),
                       class_label, genes_detected_label,"genes_detected_snRNA10xv3B_",classnames[i])
}
for (i in 1:length(classnames)) {
  plot_cluster_metrics(snRNA_10xv3_B_met, class_label, c("GABAergic","Glutamatergic","Non-Neuronal","Non-Neural"),
                       class_label, unique.counts,"UMIs_snRNA10xv3B_",classnames[i])
}

# Plot all glut subclass node plots

for (i in 1:length(glutsubnames)) {
  plot_cluster_metrics(snRNA_10xv3_B_met, class_label, c("Glutamatergic"),
                       subclass_label, unique.counts,"UMIs_snRNA10xv3B_",glutsubnames[i])
}

for (i in 1:length(glutsubnames)) {
  plot_cluster_metrics(snRNA_10xv3_B_met, class_label, c("Glutamatergic"),
                       subclass_label, genes_detected_label,"genes_detected_snRNA10xv3B_",glutsubnames[i])
}



# Plot all GABA subclass node plots

for (i in 1:length(gabasubnames)) {
  plot_cluster_metrics(snRNA_10xv3_B_met, class_label, c("GABAergic"),
                       subclass_label, unique.counts,"UMIs_snRNA10xv3B_",gabasubnames[i])
}

for (i in 1:length(gabasubnames)) {
  plot_cluster_metrics(snRNA_10xv3_B_met, class_label, c("GABAergic"),
                       subclass_label, genes_detected_label,"genes_detected_snRNA10xv3B_",gabasubnames[i])
}



# Plot all Non-neural subclass node plots

for (i in 1:length(nonneuronalsubnames)) {
  plot_cluster_metrics(snRNA_10xv3_B_met, class_label, c("Non-Neural"),
                       subclass_label, unique.counts,"UMIs_snRNA10xv3B_",nonneuralsubnames[i])
}

for (i in 1:length(nonneuralsubnames)) {
  plot_cluster_metrics(snRNA_10xv3_B_met, class_label, c("Non-Neural"),
                       subclass_label, genes_detected_label,"genes_detected_snRNA10xv3B_",nonneuralsubnames[i])
}


# Plot all Non-neuronal subclass node plots

for (i in 1:length(nonneuronalsubnames)) {
  plot_cluster_metrics(snRNA_10xv3_B_met, class_label, c("Non-Neuronal"),
                       subclass_label, unique.counts,"UMIs_snRNA10xv3B_",nonneuronalsubnames[i])
}

for (i in 1:length(nonneuronalsubnames)) {
  plot_cluster_metrics(snRNA_10xv3_B_met, class_label, c("Non-Neuronal"),
                       subclass_label, genes_detected_label,"genes_detected_snRNA10xv3B_",nonneuronalsubnames[i])
}









# Plot types yikes
setwd("/allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr/snRNA 10xV3 Broad/cluster metrics/test")


for (i in 1:length(subclassnames)) {
  temp <- snRNA_10xv3_B_met %>% filter(subclass_label %in% subclassnames[i])
  temptypes <- unique(temp$cluster_label)
  for (j in 1:length(temptypes)) {
    plot_cluster_metrics(snRNA_10xv3_B_met, subclass_label, subclassnames[i],
                         cluster_label, unique.counts,"UMIs_snRNA10xv3B_",temptypes[j])
  }
}


for (i in 1:length(subclassnames)) {
  temp <- snRNA_10xv3_B_met %>% filter(subclass_label %in% subclassnames[i])
  temptypes <- unique(temp$cluster_label)
  for (j in 1:length(temptypes)) {
    plot_cluster_metrics(snRNA_10xv3_B_met, subclass_label, subclassnames[i],
                         cluster_label, genes_detected_label,"genes_detected_snRNA10xv3B_",temptypes[j])
  }
}

#######

#######
# snRNA 10X v3 A
########

setwd("/allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr/snRNA 10xV3 AIBS/cluster metrics")
# Plot all class node plots 

for (i in 1:length(classnames)) {
  plot_cluster_metrics(snRNA_10xv3_A_met, class_label, c("GABAergic","Glutamatergic","Non-Neuronal","Non-Neural"),
                       class_label, genes_detected_label,"genes_detected_snRNA10xv3A_",classnames[i])
}
for (i in 1:length(classnames)) {
  plot_cluster_metrics(snRNA_10xv3_A_met, class_label, c("GABAergic","Glutamatergic","Non-Neuronal","Non-Neural"),
                       class_label, unique.counts,"UMIs_snRNA10xv3A_",classnames[i])
}

# Plot all glut subclass node plots

for (i in 1:length(glutsubnames)) {
  plot_cluster_metrics(snRNA_10xv3_A_met, class_label, c("Glutamatergic"),
                       subclass_label, unique.counts,"UMIs_snRNA10xv3A_",glutsubnames[i])
}

for (i in 1:length(glutsubnames)) {
  plot_cluster_metrics(snRNA_10xv3_A_met, class_label, c("Glutamatergic"),
                       subclass_label, genes_detected_label,"genes_detected_snRNA10xv3A_",glutsubnames[i])
}



# Plot all GABA subclass node plots

for (i in 1:length(gabasubnames)) {
  plot_cluster_metrics(snRNA_10xv3_A_met, class_label, c("GABAergic"),
                       subclass_label, unique.counts,"UMIs_snRNA10xv3A_",gabasubnames[i])
}

for (i in 1:length(gabasubnames)) {
  plot_cluster_metrics(snRNA_10xv3_A_met, class_label, c("GABAergic"),
                       subclass_label, genes_detected_label,"genes_detected_snRNA10xv3A_",gabasubnames[i])
}



# Plot all Non-neural subclass node plots

for (i in 1:length(nonneuralsubnames)) {
  plot_cluster_metrics(snRNA_10xv3_A_met, class_label, c("Non-Neural"),
                       subclass_label, unique.counts,"UMIs_snRNA10xv3A_",nonneuralsubnames[i])
}

for (i in 1:length(nonneuralsubnames)) {
  plot_cluster_metrics(snRNA_10xv3_A_met, class_label, c("Non-Neural"),
                       subclass_label, genes_detected_label,"genes_detected_snRNA10xv3A_",nonneuralsubnames[i])
}


# Plot all Non-neuronal subclass node plots

for (i in 1:length(nonneuronalsubnames)) {
  plot_cluster_metrics(snRNA_10xv3_A_met, class_label, c("Non-Neuronal"),
                       subclass_label, unique.counts,"UMIs_snRNA10xv3A_",nonneuronalsubnames[i])
}

for (i in 1:length(nonneuronalsubnames)) {
  plot_cluster_metrics(snRNA_10xv3_A_met, class_label, c("Non-Neuronal"),
                       subclass_label, genes_detected_label,"genes_detected_snRNA10xv3A_",nonneuronalsubnames[i])
}









# Plot types yikes
setwd("/allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr/snRNA 10xV3 AIBS/cluster metrics/test")


for (i in 1:length(subclassnames)) {
  temp <- snRNA_10xv3_A_met %>% filter(subclass_label %in% subclassnames[i])
  temptypes <- unique(temp$cluster_label)
  for (j in 1:length(temptypes)) {
    plot_cluster_metrics(snRNA_10xv3_A_met, subclass_label, subclassnames[i],
                         cluster_label, unique.counts,"UMIs_snRNA10xv3A_",temptypes[j])
  }
}


for (i in 1:length(subclassnames)) {
  temp <- snRNA_10xv3_A_met %>% filter(subclass_label %in% subclassnames[i])
  temptypes <- unique(temp$cluster_label)
  for (j in 1:length(temptypes)) {
    plot_cluster_metrics(snRNA_10xv3_A_met, subclass_label, subclassnames[i],
                         cluster_label, genes_detected_label,"genes_detected_snRNA10xv3A_",temptypes[j])
  }
}

########
# snRNA 10X v2 A
########

setwd("/allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr/snRNA 10xV2 AIBS/cluster metrics")
# Plot all class node plots 

for (i in 1:length(classnames)) {
  plot_cluster_metrics(snRNA_10xv2_A_met, class_label, c("GABAergic","Glutamatergic","Non-Neuronal","Non-Neural"),
                       class_label, genes_detected_label,"genes_detected_snRNA10xv2A_",classnames[i])
}
for (i in 1:length(classnames)) {
  plot_cluster_metrics(snRNA_10xv2_A_met, class_label, c("GABAergic","Glutamatergic","Non-Neuronal","Non-Neural"),
                       class_label, unique.counts,"UMIs_snRNA10xv2A_",classnames[i])
}

# Plot all glut subclass node plots

for (i in 1:length(glutsubnames)) {
  plot_cluster_metrics(snRNA_10xv2_A_met, class_label, c("Glutamatergic"),
                       subclass_label, unique.counts,"UMIs_snRNA10xv2A_",glutsubnames[i])
}

for (i in 1:length(glutsubnames)) {
  plot_cluster_metrics(snRNA_10xv2_A_met, class_label, c("Glutamatergic"),
                       subclass_label, genes_detected_label,"genes_detected_snRNA10xv2A_",glutsubnames[i])
}



# Plot all GABA subclass node plots

for (i in 1:length(gabasubnames)) {
  plot_cluster_metrics(snRNA_10xv2_A_met, class_label, c("GABAergic"),
                       subclass_label, unique.counts,"UMIs_snRNA10xv2A_",gabasubnames[i])
}

for (i in 1:length(gabasubnames)) {
  plot_cluster_metrics(snRNA_10xv2_A_met, class_label, c("GABAergic"),
                       subclass_label, genes_detected_label,"genes_detected_snRNA10xv2A_",gabasubnames[i])
}



# Plot all Non-neural subclass node plots

for (i in 1:length(nonneuralsubnames)) {
  plot_cluster_metrics(snRNA_10xv2_A_met, class_label, c("Non-Neural"),
                       subclass_label, unique.counts,"UMIs_snRNA10xv2A_",nonneuralsubnames[i])
}

for (i in 1:length(nonneuralsubnames)) {
  plot_cluster_metrics(snRNA_10xv2_A_met, class_label, c("Non-Neural"),
                       subclass_label, genes_detected_label,"genes_detected_snRNA10xv2A_",nonneuralsubnames[i])
}


# Plot all Non-neuronal subclass node plots

for (i in 1:length(nonneuronalsubnames)) {
  plot_cluster_metrics(snRNA_10xv2_A_met, class_label, c("Non-Neuronal"),
                       subclass_label, unique.counts,"UMIs_snRNA10xv2A_",nonneuronalsubnames[i])
}

for (i in 1:length(nonneuronalsubnames)) {
  plot_cluster_metrics(snRNA_10xv2_A_met, class_label, c("Non-Neuronal"),
                       subclass_label, genes_detected_label,"genes_detected_snRNA10xv2A_",nonneuronalsubnames[i])
}









# Plot types yikes
setwd("/allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr/snRNA 10xV2 AIBS/cluster metrics/test")


for (i in 1:length(subclassnames)) {
  temp <- snRNA_10xv2_A_met %>% filter(subclass_label %in% subclassnames[i])
  temptypes <- unique(temp$cluster_label)
  for (j in 1:length(temptypes)) {
    plot_cluster_metrics(snRNA_10xv2_A_met, subclass_label, subclassnames[i],
                         cluster_label, unique.counts,"UMIs_snRNA10xv2A_",temptypes[j])
  }
}


for (i in 1:length(subclassnames)) {
  temp <- snRNA_10xv2_A_met %>% filter(subclass_label %in% subclassnames[i])
  temptypes <- unique(temp$cluster_label)
  for (j in 1:length(temptypes)) {
    plot_cluster_metrics(snRNA_10xv2_A_met, subclass_label, subclassnames[i],
                         cluster_label, genes_detected_label,"genes_detected_snRNA10xv2A_",temptypes[j])
  }
}


#########
# scRNA 10X v3 A
########
setwd("/allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr/scRNA 10xV3 AIBS/cluster metrics")
# Plot all class node plots 

for (i in 1:length(classnames)) {
  plot_cluster_metrics(scRNA_10xv3_A_met, class_label, c("GABAergic","Glutamatergic","Non-Neuronal","Non-Neural"),
                       class_label, genes_detected_label,"genes_detected_scRNA10xv3A_",classnames[i])
}
for (i in 1:length(classnames)) {
  plot_cluster_metrics(scRNA_10xv3_A_met, class_label, c("GABAergic","Glutamatergic","Non-Neuronal","Non-Neural"),
                       class_label, unique.counts,"UMIs_scRNA10xv3A_",classnames[i])
}

# Plot all glut subclass node plots

for (i in 1:length(glutsubnames)) {
  plot_cluster_metrics(scRNA_10xv3_A_met, class_label, c("Glutamatergic"),
                       subclass_label, unique.counts,"UMIs_scRNA10xv3A_",glutsubnames[i])
}

for (i in 1:length(glutsubnames)) {
  plot_cluster_metrics(scRNA_10xv3_A_met, class_label, c("Glutamatergic"),
                       subclass_label, genes_detected_label,"genes_detected_scRNA10xv3A_",glutsubnames[i])
}



# Plot all GABA subclass node plots

for (i in 1:length(gabasubnames)) {
  plot_cluster_metrics(scRNA_10xv3_A_met, class_label, c("GABAergic"),
                       subclass_label, unique.counts,"UMIs_scRNA10xv3A_",gabasubnames[i])
}

for (i in 1:length(gabasubnames)) {
  plot_cluster_metrics(scRNA_10xv3_A_met, class_label, c("GABAergic"),
                       subclass_label, genes_detected_label,"genes_detected_scRNA10xv3A_",gabasubnames[i])
}



# Plot all Non-neural subclass node plots

for (i in 1:length(nonneuralsubnames)) {
  plot_cluster_metrics(scRNA_10xv3_A_met, class_label, c("Non-Neural"),
                       subclass_label, unique.counts,"UMIs_scRNA10xv3A_",nonneuralsubnames[i])
}

for (i in 1:length(nonneuralsubnames)) {
  plot_cluster_metrics(scRNA_10xv3_A_met, class_label, c("Non-Neural"),
                       subclass_label, genes_detected_label,"genes_detected_scRNA10xv3A_",nonneuralsubnames[i])
}


# Plot all Non-neuronal subclass node plots

for (i in 1:length(nonneuronalsubnames)) {
  plot_cluster_metrics(scRNA_10xv3_A_met, class_label, c("Non-Neuronal"),
                       subclass_label, unique.counts,"UMIs_scRNA10xv3A_",nonneuronalsubnames[i])
}

for (i in 1:length(nonneuronalsubnames)) {
  plot_cluster_metrics(scRNA_10xv3_A_met, class_label, c("Non-Neuronal"),
                       subclass_label, genes_detected_label,"genes_detected_scRNA10xv3A_",nonneuronalsubnames[i])
}









# Plot types yikes
setwd("/allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr/scRNA 10xV3 AIBS/cluster metrics/test")


for (i in 1:length(subclassnames)) {
  temp <- scRNA_10xv3_A_met %>% filter(subclass_label %in% subclassnames[i])
  temptypes <- unique(temp$cluster_label)
  for (j in 1:length(temptypes)) {
    plot_cluster_metrics(scRNA_10xv3_A_met, subclass_label, subclassnames[i],
                         cluster_label, unique.counts,"UMIs_scRNA10xv3A_",temptypes[j])
  }
}


for (i in 1:length(subclassnames)) {
  temp <- scRNA_10xv3_A_met %>% filter(subclass_label %in% subclassnames[i])
  temptypes <- unique(temp$cluster_label)
  for (j in 1:length(temptypes)) {
    plot_cluster_metrics(scRNA_10xv3_A_met, subclass_label, subclassnames[i],
                         cluster_label, genes_detected_label,"genes_detected_scRNA10xv3A_",temptypes[j])
  }
}


#######
# scRNA 10X v2 A
######

setwd("/allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr/scRNA 10xV2 AIBS/cluster metrics")
# Plot all class node plots 

for (i in 1:length(classnames)) {
  plot_cluster_metrics(scRNA_10xv2_A_met, class_label, c("GABAergic","Glutamatergic","Non-Neuronal","Non-Neural"),
                       class_label, genes_detected_label,"genes_detected_scRNA10xv2A_",classnames[i])
}
for (i in 1:length(classnames)) {
  plot_cluster_metrics(scRNA_10xv2_A_met, class_label, c("GABAergic","Glutamatergic","Non-Neuronal","Non-Neural"),
                       class_label, unique.counts,"UMIs_scRNA10xv2A_",classnames[i])
}

# Plot all glut subclass node plots

for (i in 1:length(glutsubnames)) {
  plot_cluster_metrics(scRNA_10xv2_A_met, class_label, c("Glutamatergic"),
                       subclass_label, unique.counts,"UMIs_scRNA10xv2A_",glutsubnames[i])
}

for (i in 1:length(glutsubnames)) {
  plot_cluster_metrics(scRNA_10xv2_A_met, class_label, c("Glutamatergic"),
                       subclass_label, genes_detected_label,"genes_detected_scRNA10xv2A_",glutsubnames[i])
}



# Plot all GABA subclass node plots

for (i in 1:length(gabasubnames)) {
  plot_cluster_metrics(scRNA_10xv2_A_met, class_label, c("GABAergic"),
                       subclass_label, unique.counts,"UMIs_scRNA10xv2A_",gabasubnames[i])
}

for (i in 1:length(gabasubnames)) {
  plot_cluster_metrics(scRNA_10xv2_A_met, class_label, c("GABAergic"),
                       subclass_label, genes_detected_label,"genes_detected_scRNA10xv2A_",gabasubnames[i])
}



# Plot all Non-neural subclass node plots

for (i in 1:length(nonneuralsubnames)) {
  plot_cluster_metrics(scRNA_10xv2_A_met, class_label, c("Non-Neural"),
                       subclass_label, unique.counts,"UMIs_scRNA10xv2A_",nonneuralsubnames[i])
}

for (i in 1:length(nonneuralsubnames)) {
  plot_cluster_metrics(scRNA_10xv2_A_met, class_label, c("Non-Neural"),
                       subclass_label, genes_detected_label,"genes_detected_scRNA10xv2A_",nonneuralsubnames[i])
}


# Plot all Non-neuronal subclass node plots

for (i in 1:length(nonneuronalsubnames)) {
  plot_cluster_metrics(scRNA_10xv2_A_met, class_label, c("Non-Neuronal"),
                       subclass_label, unique.counts,"UMIs_scRNA10xv2A_",nonneuronalsubnames[i])
}

for (i in 1:length(nonneuronalsubnames)) {
  plot_cluster_metrics(scRNA_10xv2_A_met, class_label, c("Non-Neuronal"),
                       subclass_label, genes_detected_label,"genes_detected_scRNA10xv2A_",nonneuronalsubnames[i])
}









# Plot types yikes
setwd("/allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr/scRNA 10xV2 AIBS/cluster metrics/test")


for (i in 1:length(subclassnames)) {
  temp <- scRNA_10xv2_A_met %>% filter(subclass_label %in% subclassnames[i])
  temptypes <- unique(temp$cluster_label)
  for (j in 1:length(temptypes)) {
    plot_cluster_metrics(scRNA_10xv2_A_met, subclass_label, subclassnames[i],
                         cluster_label, unique.counts,"UMIs_scRNA10xv2A_",temptypes[j])
  }
}


for (i in 1:length(subclassnames)) {
  temp <- scRNA_10xv2_A_met %>% filter(subclass_label %in% subclassnames[i])
  temptypes <- unique(temp$cluster_label)
  for (j in 1:length(temptypes)) {
    plot_cluster_metrics(scRNA_10xv2_A_met, subclass_label, subclassnames[i],
                         cluster_label, genes_detected_label,"genes_detected_scRNA10xv2A_",temptypes[j])
  }
}





########
# scRNA SMART
######

##########





setwd("/allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr/scRNA SMART/cluster metrics")
# Plot all class node plots 

for (i in 1:length(classnames)) {
  plot_cluster_metrics(smartcellmet, class_label, c("GABAergic","Glutamatergic","Non-Neuronal","Non-Neural"),
                       class_label, genes_detected_label,"genes_detected_scRNASMART_",classnames[i])
}
for (i in 1:length(classnames)) {
  plot_cluster_metrics(smartcellmet, class_label, c("GABAergic","Glutamatergic","Non-Neuronal","Non-Neural"),
                       class_label, unique.counts,"UMIs_scRNASMART_",classnames[i])
}

# Plot all glut subclass node plots

for (i in 1:length(glutsubnames)) {
  plot_cluster_metrics(smartcellmet, class_label, c("Glutamatergic"),
                       subclass_label, unique.counts,"UMIs_scRNASMART_",glutsubnames[i])
}

for (i in 1:length(glutsubnames)) {
  plot_cluster_metrics(smartcellmet, class_label, c("Glutamatergic"),
                       subclass_label, genes_detected_label,"genes_detected_scRNASMART_",glutsubnames[i])
}



# Plot all GABA subclass node plots

for (i in 1:length(gabasubnames)) {
  plot_cluster_metrics(smartcellmet, class_label, c("GABAergic"),
                       subclass_label, unique.counts,"UMIs_scRNASMART_",gabasubnames[i])
}

for (i in 1:length(gabasubnames)) {
  plot_cluster_metrics(smartcellmet, class_label, c("GABAergic"),
                       subclass_label, genes_detected_label,"genes_detected_scRNASMART_",gabasubnames[i])
}



# Plot all Non-neural subclass node plots

for (i in 1:length(nonneuralsubnames)) {
  plot_cluster_metrics(smartcellmet, class_label, c("Non-Neural"),
                       subclass_label, unique.counts,"UMIs_scRNASMART_",nonneuralsubnames[i])
}

for (i in 1:length(nonneuralsubnames)) {
  plot_cluster_metrics(smartcellmet, class_label, c("Non-Neural"),
                       subclass_label, genes_detected_label,"genes_detected_scRNASMART_",nonneuralsubnames[i])
}


# Plot all Non-neuronal subclass node plots

for (i in 1:length(nonneuronalsubnames)) {
  plot_cluster_metrics(smartcellmet, class_label, c("Non-Neuronal"),
                       subclass_label, unique.counts,"UMIs_scRNASMART_",nonneuronalsubnames[i])
}

for (i in 1:length(nonneuronalsubnames)) {
  plot_cluster_metrics(smartcellmet, class_label, c("Non-Neuronal"),
                       subclass_label, genes_detected_label,"genes_detected_scRNASMART_",nonneuronalsubnames[i])
}





# Plot types yikes
setwd("/allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr/scRNA SMART/cluster metrics/test")


for (i in 1:length(subclassnames)) {
  temp <- smartcellmet %>% filter(subclass_label %in% subclassnames[i])
  temptypes <- unique(temp$cluster_label)
  for (j in 1:length(temptypes)) {
    plot_cluster_metrics(smartcellmet, subclass_label, subclassnames[i],
                         cluster_label, unique.counts,"UMIs_scRNASMART_",temptypes[j])
  }
}

for (i in 1:length(subclassnames)) {
  temp <- smartcellmet %>% filter(subclass_label %in% subclassnames[i])
  temptypes <- unique(temp$cluster_label)
  for (j in 1:length(temptypes)) {
    plot_cluster_metrics(smartcellmet, subclass_label, subclassnames[i],
                         cluster_label, genes_detected_label,"genes_detected_scRNASMART_",temptypes[j])
  }
}
