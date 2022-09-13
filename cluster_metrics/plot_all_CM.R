# Takin a stab at functionizing this whole thing


# Part for loading the data
options(stringsAsFactors = F)
setwd("//allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr/scRNA 10xV2 AIBS")

allmet <- read_csv("scRNA_10xv2_A Metadata.csv")
setwd("//allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr/")
clupdate <- read_csv("cl_updated.df.csv")

outdir <- "//allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr/automated_cards"

make_cluster_metrics <- function(metadata,nomenclature,outdir,datasetname) {
  
  suppressPackageStartupMessages({
    library(scrattch.vis)
    library(tidyverse)
    library(mgsub)
    library(cowplot)
    library(feather)
    library(ggnewscale)
  })
  
setwd("//allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr/clustermetrics_new")
source("new_CM_functions.R")

setwd(outdir)
allmet <- metadata
clupdate <- nomenclature


# Part for extracting relevant things from the cl_update file. Need to remove all traces of 
# being class or subclass specific

typenames <- unique(clupdate$cluster_label)
typeids <- cbind.data.frame(typenames, unique(clupdate$cluster_id), unique(clupdate$cluster_color))

colnames(typeids) <- c("cluster_label","cluster_id","cluster_color")
typeids <- typeids[order(typeids$cluster_id),]
typenames <- typeids$cluster_label

subclassnames <- unique(clupdate$subclass_label)
subids <- cbind.data.frame(subclassnames, unique(clupdate$subclass_id),unique(clupdate$subclass_color))
colnames(subids) <- c("subclass_label","subclass_id","subclass_color")
subids <- subids[order(subids$subclass_id),]
subclassnames <-subids$subclass_label


classnames <- unique(clupdate$class_label)
classids <- cbind.data.frame(classnames, unique(clupdate$class_id),unique(clupdate$class_color))
colnames(classids) <- c("class_label","class_id","class_color")



setwd("//allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr/automated_cards")

# Plot class level with taxons highlighted

for (i in 1:length(classnames)) {
  plot_cluster_metrics(data = allmet, 
                       filter_level = class_label, 
                       filterby = classnames,
                       X = class_label, 
                       Y = genes_detected_label, 
                       level_id = class_id,
                       outname = paste0("genes_detected_",datname,"_"),
                       nodename = classnames[i],
                       setnames = classnames,
                       setids = classids$class_id,
                       setcolors = classids$class_color)
  
  plot_cluster_metrics(data = allmet, 
                       filter_level = class_label, 
                       filterby = classnames,
                       X = class_label, 
                       Y = unique.counts, 
                       level_id = class_id,
                       outname = paste0("UMIs_",datname,"_"),
                       nodename = classnames[i],
                       setnames = classnames,
                       setids = classids$class_id,
                       setcolors = classids$class_color)
}

# Plot class level no highlights
plot_cluster_metrics_nohighlight(data = allmet, 
                     filter_level = class_label, 
                     filterby = classnames,
                     X = class_label, 
                     Y = genes_detected_label, 
                     level_id = class_id,
                     outname = paste0("genes_detected_",datname,"_"),
                     nodename = "AllClasses_child_view",
                     setnames = classnames,
                     setids = classids$class_id,
                     setcolors = classids$class_color)



plot_cluster_metrics_nohighlight(data = allmet, 
                     filter_level = class_label, 
                     filterby = classnames,
                     X = class_label, 
                     Y = unique.counts, 
                     level_id = class_id,
                     outname = paste0("UMIs_",datname,"_"),
                     nodename = "AllClasses_child_view",
                     setnames = classnames,
                     setids = classids$class_id,
                     setcolors = classids$class_color)


############

# Plot all subclasses in one nested for loop

for (i in 1:length(classnames)) {
  tempsubs <- unique(clupdate$subclass_label[clupdate$class_label == classnames[i]])
  tempdat <- unique(clupdate[clupdate$class_label == classnames[i],])
  tempdat <- unique(tempdat %>% select(subclass_id,subclass_label,subclass_color)) 
  
  for (j in 1:length(tempsubs)) {
    plot_cluster_metrics(data = allmet, 
                        filter_level = class_label, 
                        filterby = classnames[i],
                        X = subclass_label, 
                        Y = unique.counts, 
                        level_id = subclass_id,
                        outname = paste0("UMIs_",datname,"_"),
                        nodename = tempsubs[j],
                        setnames = tempsubs,
                        setids = tempdat$subclass_id,
                        setcolors = tempdat$subclass_color)
    
    plot_cluster_metrics(data = allmet, 
                         filter_level = class_label, 
                         filterby = classnames[i],
                         X = subclass_label, 
                         Y = genes_detected_label, 
                         level_id = subclass_id,
                         outname = paste0("genes_detected_",datname,"_"),
                         nodename = tempsubs[j],
                         setnames = tempsubs,
                         setids = tempdat$subclass_id,
                         setcolors = tempdat$subclass_color)
    
  }
  
  plot_cluster_metrics_nohighlight(data = allmet, 
                       filter_level = class_label, 
                       filterby = classnames[i],
                       X = subclass_label, 
                       Y = unique.counts, 
                       level_id = subclass_id,
                       outname = paste0("UMIs_",datname,"_"),
                       nodename = paste0(classnames[i],"_child_view"),
                       setnames = tempsubs,
                       setids = tempdat$subclass_id,
                       setcolors = tempdat$subclass_color)
  
  plot_cluster_metrics_nohighlight(data = allmet, 
                                   filter_level = class_label, 
                                   filterby = classnames[i],
                                   X = subclass_label, 
                                   Y = genes_detected_label, 
                                   level_id = subclass_id,
                                   outname = paste0("genes_detected_",datname,"_"),
                                   nodename = paste0(classnames[i],"_child_view"),
                                   setnames = tempsubs,
                                   setids = tempdat$subclass_id,
                                   setcolors = tempdat$subclass_color)
  
}





########

# Plot types 
missingsubs <- setdiff(as.character(subclassnames), unique(as.character(allmet$subclass_label)))
missingnodes <- setdiff(as.character(typenames), unique(as.character(allmet$cluster_label)))

goodsub <- setdiff(as.character(subclassnames),missingsubs)

# UMIs
for (i in 1:length(goodsub)) {
  temp <- allmet %>% filter(subclass_label %in% goodsub[i])
  temp <- temp[order(temp$cluster_id),]
  checktypes <- setdiff(clupdate$cluster_label[clupdate$subclass_label == goodsub[i]], unique(as.character(temp$cluster_label)))
  subsize <- length(unique(temp$cluster_label[temp$subclass_label == goodsub[i]]))
  
  if(length(checktypes) == 0 && subsize > 1) {
    temptypes <- as.character(unique(temp$cluster_label))
    for (j in 1:length(temptypes)) {
      ttnames <- temptypes
      ttids <- clupdate %>% filter(cluster_label %in% temptypes)
      plot_cluster_metrics(data = allmet, 
                           filter_level = subclass_label, 
                           filterby = goodsub[i],
                           X = cluster_label, 
                           Y = unique.counts, 
                           level_id = cluster_id,
                           outname = paste0("UMIs_",datasetname,"_"),
                           nodename = ttnames[j],
                           setnames = ttnames,
                           setids = ttids$cluster_id,
                           setcolors = ttids$cluster_color)
      
    }
  }
  
  if (subsize == 1 & length(checktypes) == 0) {
    temptypes <- as.character(unique(temp$cluster_label))
    whichparent <- unique(clupdate$class_label[clupdate$subclass_label == goodsub[i]])
    whichsubs <- unique(clupdate$subclass_label[clupdate$class_label == whichparent])
    whichids <- unique(clupdate$subclass_id[clupdate$class_label == whichparent])
    whichcols <- unique(clupdate$subclass_color[clupdate$class_label == whichparent])
    plot_cluster_metrics(data = allmet, 
                         filter_level = class_label, 
                         filterby = whichparent,
                         X = subclass_label, 
                         Y = unique.counts, 
                         level_id = subclass_id,
                         outname = paste0("UMIs_",datasetname,"TYPE_"),
                         nodename = goodsub[i],
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
      plot_cluster_metrics(data = allmet, 
                           filter_level = subclass_label, 
                           filterby = goodsub[i],
                           X = cluster_label, 
                           Y = unique.counts, 
                           level_id = cluster_id,
                           outname = paste0("UMIs_",datasetname,"_"),
                           nodename = ttnames[j],
                           setnames = ttnames,
                           setids = whichids,
                           setcolors = whichcols)
      
    }
  }
  
  
}

# No highlight

for (i in 1:length(goodsub)) {
  temp <- allmet %>% filter(subclass_label %in% goodsub[i])
  temp <- temp[order(temp$cluster_id),]
  checktypes <- setdiff(clupdate$cluster_label[clupdate$subclass_label == goodsub[i]], unique(as.character(temp$cluster_label)))
  subsize <- length(unique(temp$cluster_label[temp$subclass_label == goodsub[i]]))
  
  if(length(checktypes) == 0 && subsize > 1) {
    temptypes <- as.character(unique(temp$cluster_label))
    ttnames <- temptypes
    ttids <- clupdate %>% filter(cluster_label %in% temptypes)
    plot_cluster_metrics_nohighlight(data = allmet, 
                                     filter_level = subclass_label, 
                                     filterby = goodsub[i],
                                     X = cluster_label, 
                                     Y = unique.counts, 
                                     level_id = cluster_id,
                                     outname = paste0("UMIs_",datasetname,"_"),
                                     nodename = paste0(goodsub[i],"_child_view"),
                                     setnames = ttnames,
                                     setids = ttids$cluster_id,
                                     setcolors = ttids$cluster_color)
    
    
  }
  
  if (subsize == 1 & length(checktypes) == 0) {
    temptypes <- as.character(unique(temp$cluster_label))
    whichparent <- unique(clupdate$class_label[clupdate$subclass_label == goodsub[i]])
    whichsubs <- unique(clupdate$subclass_label[clupdate$class_label == whichparent])
    whichids <- unique(clupdate$subclass_id[clupdate$class_label == whichparent])
    whichcols <- unique(clupdate$subclass_color[clupdate$class_label == whichparent])
    plot_cluster_metrics_nohighlight(data = allmet, 
                                     filter_level = class_label, 
                                     filterby = whichparent,
                                     X = subclass_label, 
                                     Y = unique.counts, 
                                     level_id = subclass_id,
                                     outname = paste0("UMIs_",datasetname,"TYPE_"),
                                     nodename = paste0(goodsub[i],"_child_view"),
                                     setnames = whichsubs,
                                     setids = whichids,
                                     setcolors = whichcols)
  }
  
  
  
  if (length(checktypes > 0)) {
    temptypes <- c(as.character(unique(temp$cluster_label)), checktypes)
    whichids <- clupdate$cluster_id[clupdate$cluster_label %in% temptypes]
    whichtypes <- clupdate$cluster_label[clupdate$cluster_label %in% temptypes]
    whichcols <- clupdate$cluster_color[clupdate$cluster_label %in% temptypes]
    ttnames <- whichtypes
    ttids <- clupdate %>% filter(cluster_label %in% temptypes)
    plot_cluster_metrics_nohighlight(data = allmet, 
                                     filter_level = subclass_label, 
                                     filterby = goodsub[i],
                                     X = cluster_label, 
                                     Y = unique.counts, 
                                     level_id = cluster_id,
                                     outname = paste0("UMIs_",datasetname,"_"),
                                     nodename = paste0(goodsub[i],"_child_view"),
                                     setnames = ttnames,
                                     setids = whichids,
                                     setcolors = whichcols)
    
    
  }
  
}

# Genes detected
for (i in 1:length(goodsub)) {
  temp <- allmet %>% filter(subclass_label %in% goodsub[i])
  temp <- temp[order(temp$cluster_id),]
  checktypes <- setdiff(clupdate$cluster_label[clupdate$subclass_label == goodsub[i]], unique(as.character(temp$cluster_label)))
  subsize <- length(unique(temp$cluster_label[temp$subclass_label == goodsub[i]]))
  
  if(length(checktypes) == 0 && subsize > 1) {
    temptypes <- as.character(unique(temp$cluster_label))
    for (j in 1:length(temptypes)) {
      ttnames <- temptypes
      ttids <- clupdate %>% filter(cluster_label %in% temptypes)
      plot_cluster_metrics(data = allmet, 
                           filter_level = subclass_label, 
                           filterby = goodsub[i],
                           X = cluster_label, 
                           Y = genes_detected_label, 
                           level_id = cluster_id,
                           outname = paste0("genes_detected_",datasetname,"_"),
                           nodename = ttnames[j],
                           setnames = ttnames,
                           setids = ttids$cluster_id,
                           setcolors = ttids$cluster_color)
      
    }
  }
  
  if (subsize == 1 & length(checktypes) == 0) {
    temptypes <- as.character(unique(temp$cluster_label))
    whichparent <- unique(clupdate$class_label[clupdate$subclass_label == goodsub[i]])
    whichsubs <- unique(clupdate$subclass_label[clupdate$class_label == whichparent])
    whichids <- unique(clupdate$subclass_id[clupdate$class_label == whichparent])
    whichcols <- unique(clupdate$subclass_color[clupdate$class_label == whichparent])
    plot_cluster_metrics(data = allmet, 
                         filter_level = class_label, 
                         filterby = whichparent,
                         X = subclass_label, 
                         Y = genes_detected_label, 
                         level_id = subclass_id,
                         outname = paste0("genes_detected_",datasetname,"TYPE_"),
                         nodename = goodsub[i],
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
      plot_cluster_metrics(data = allmet, 
                           filter_level = subclass_label, 
                           filterby = goodsub[i],
                           X = cluster_label, 
                           Y = genes_detected_label, 
                           level_id = cluster_id,
                           outname = paste0("genes_detected_",datasetname,"_"),
                           nodename = ttnames[j],
                           setnames = ttnames,
                           setids = whichids,
                           setcolors = whichcols)
      
    }
  }
  
  
}

# No highlight
for (i in 1:length(goodsub)) {
  temp <- allmet %>% filter(subclass_label %in% goodsub[i])
  temp <- temp[order(temp$cluster_id),]
  checktypes <- setdiff(clupdate$cluster_label[clupdate$subclass_label == goodsub[i]], unique(as.character(temp$cluster_label)))
  subsize <- length(unique(temp$cluster_label[temp$subclass_label == goodsub[i]]))
  
  if(length(checktypes) == 0 && subsize > 1) {
    temptypes <- as.character(unique(temp$cluster_label))
    ttnames <- temptypes
    ttids <- clupdate %>% filter(cluster_label %in% temptypes)
    plot_cluster_metrics_nohighlight(data = allmet, 
                                     filter_level = subclass_label, 
                                     filterby = goodsub[i],
                                     X = cluster_label, 
                                     Y = genes_detected_label, 
                                     level_id = cluster_id,
                                     outname = paste0("genes_detected_",datasetname,"_"),
                                     nodename = paste0(goodsub[i],"_child_view"),
                                     setnames = ttnames,
                                     setids = ttids$cluster_id,
                                     setcolors = ttids$cluster_color)
    
    
  }
  
  if (subsize == 1 & length(checktypes) == 0) {
    temptypes <- as.character(unique(temp$cluster_label))
    whichparent <- unique(clupdate$class_label[clupdate$subclass_label == goodsub[i]])
    whichsubs <- unique(clupdate$subclass_label[clupdate$class_label == whichparent])
    whichids <- unique(clupdate$subclass_id[clupdate$class_label == whichparent])
    whichcols <- unique(clupdate$subclass_color[clupdate$class_label == whichparent])
    plot_cluster_metrics_nohighlight(data = allmet, 
                                     filter_level = class_label, 
                                     filterby = whichparent,
                                     X = subclass_label, 
                                     Y = genes_detected_label, 
                                     level_id = subclass_id,
                                     outname = paste0("genes_detected_",datasetname,"TYPE_"),
                                     nodename = paste0(goodsub[i],"_child_view"),
                                     setnames = whichsubs,
                                     setids = whichids,
                                     setcolors = whichcols)
  }
  
  
  
  if (length(checktypes > 0)) {
    temptypes <- c(as.character(unique(temp$cluster_label)), checktypes)
    whichids <- clupdate$cluster_id[clupdate$cluster_label %in% temptypes]
    whichtypes <- clupdate$cluster_label[clupdate$cluster_label %in% temptypes]
    whichcols <- clupdate$cluster_color[clupdate$cluster_label %in% temptypes]
    ttnames <- whichtypes
    ttids <- clupdate %>% filter(cluster_label %in% temptypes)
    plot_cluster_metrics_nohighlight(data = allmet, 
                                     filter_level = subclass_label, 
                                     filterby = goodsub[i],
                                     X = cluster_label, 
                                     Y = genes_detected_label, 
                                     level_id = cluster_id,
                                     outname = paste0("genes_detected_",datasetname,"_"),
                                     nodename = paste0(goodsub[i],"_child_view"),
                                     setnames = ttnames,
                                     setids = whichids,
                                     setcolors = whichcols)
    
    
  }
  
}

}


make_cluster_metrics(metadata = allmet,
                     nomenclature = clupdate,
                     outdir = outdir,
                     datasetname = "scRNA10xV2")
