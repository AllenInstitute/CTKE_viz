#' This function is a wrapper for the functions in "make_dotplot.R" and allows for the generation of all dot plots for every
#' node in a taxonomy in a single function call.
#' 
#' @param counts - standard count matrix with "sample_name" column in sparse matrix (dgTMatrix) format
#' @param metadata - standard metadata file with "sample_name" column
#' @param taxonomy_file - standard annotated cluster table providing metadata for every node in the taxonomy
#' @param genes - tidy data frame providing genes in desired plotting order. requires the following columns:
#'                 - clusterName - name of each of taxonomy
#'                 - gene - gene symbol
#'                 - NSForestMarker - "y" or "n" indicating whether the gene should be highlighted as a marker on the plot
#' @param outdir - string providing name of directory to save files to
#' @param datasetname - string providing name of dataset. Will be included in output file name



make_all_dots <- function(counts,
                          metadata,
                          taxonomy_file,
                          genes,
                          outdir,
                          datasetname) {
  
   suppressPackageStartupMessages({
     require(plyr)
     require(tidyverse)
     require(scrattch.hicat)
     require(scrattch.vis)
     require(ggfittext)
   })
  print("Making dot plots...depending on the number of clusters in your taxonomy and the size of your data, this may take several minutes.")
  
  is.sequential <- function(x){all(abs(diff(x)) == 1)}
  
  setwd(outdir)
  
  clupdate <- taxonomy_file
  
  
  # Part for extracting relevant things from the cl_update file
  
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
  
  
  
  
  
  
  # Do classes first
  classgenes <- unique(genes$gene[genes$clusterName %in% classnames])
  plotcounts <- cbind(sample_name = colnames(as.matrix(counts)),
                      as.data.frame(t(as.matrix(counts[classgenes,]))))
  
  plotcounts <- plotcounts[plotcounts$sample_name %in% metadata$sample_name,] 
  
  for (i in 1:length(classnames)) {
    print(paste0("Making dot plot for ",classnames[i]))
    geneset <- unique(genes$gene[genes$clusterName %in% classnames[i]])
    tempmark <- genes %>% filter(clusterName==classnames[i]) %>% filter(NSForestMarker=="y")
    tempmark <- tempmark$gene
    
    
    plot_all_dots(data = plotcounts,
                  anno = metadata,
                  genes = geneset,
                  grouping = "class",
                  markers = tempmark,
                  filter_by = classnames,
                  nodename = classnames[i],
                  setids = classids,
                  dataname = datasetname,
                  setnames = classnames)
  }
  
  
  
  
  # Do subclasses next
  
  
  for (i in 1:length(classnames)) {
    tempsubs <- unique(clupdate$subclass_label[clupdate$class_label == classnames[i]])
    tempdat <- unique(clupdate[clupdate$class_label == classnames[i],])
    tempdat <- unique(tempdat %>% select(subclass_id,subclass_label,subclass_color))
    
    subclassgenes <- unique(genes$gene[genes$clusterName %in% tempsubs])
    plotcounts <- cbind(sample_name = colnames(as.matrix(counts)),
                        as.data.frame(t(as.matrix(counts[subclassgenes,]))))
    
    plotcounts <- plotcounts[plotcounts$sample_name %in% metadata$sample_name,] 
    
    for (j in 1:length(tempsubs)){
      print(paste0("Making dot plot for ",tempsubs[j]))
      geneset <- subclassgenes
      tempmark <- genes %>% filter(clusterName==tempsubs[j]) %>% filter(NSForestMarker=="y")
      tempmark <- tempmark$gene
      plot_all_dots(data = plotcounts,
                    anno = metadata,
                    genes = geneset,
                    grouping = "subclass",
                    markers = tempmark,
                    filter_by = subclassnames,
                    nodename = tempsubs[j],
                    setids = tempdat$subclass_id,
                    dataname = datasetname,
                    setnames = tempsubs)
    }
  }
  
  
  # Plot types 
  missingsubs <- setdiff(as.character(subclassnames), unique(as.character(metadata$subclass_label)))
  missingnodes <- setdiff(as.character(typenames), unique(as.character(metadata$cluster_label)))
  
  goodsub <- setdiff(as.character(subclassnames),missingsubs)
  
  for (i in 1:length(goodsub)) {
    temp <- metadata %>% filter(subclass_label %in% goodsub[i])
    temp <- temp[order(temp$cluster_id),]
    checktypes <- setdiff(clupdate$cluster_label[clupdate$subclass_label == goodsub[i]], unique(as.character(temp$cluster_label)))
    subsize <- length(unique(temp$cluster_label[temp$subclass_label == goodsub[i]]))
    
    
    
    
    
    if(length(checktypes) == 0 && subsize > 1) {
      temptypes <- as.character(unique(temp$cluster_label))
      goodsubgenes <- unique(genes$gene[genes$clusterName %in% temptypes])
      plotcounts <- cbind(sample_name = colnames(as.matrix(counts)),
                          as.data.frame(t(as.matrix(counts[goodsubgenes,]))))
      
      plotcounts <- plotcounts[plotcounts$sample_name %in% metadata$sample_name,] 
      
      for (j in 1:length(temptypes)) {
        print(paste0("Making dot plot for ",temptypes[j]))
        ttnames <- temptypes
        ttids <- clupdate %>% filter(cluster_label %in% temptypes)
        geneset <- goodsubgenes
        tempmark <- genes %>% filter(clusterName==temptypes[j]) %>% filter(NSForestMarker=="y")
        tempmark <- tempmark$gene
        plot_all_dots(data = plotcounts,
                      anno = metadata,
                      genes = geneset,
                      grouping = "cluster",
                      markers = tempmark,
                      filter_by = goodsub[i],
                      nodename = temptypes[j],
                      setids = ttids$subclass_id,
                      dataname = datasetname,
                      setnames = ttnames,
                      keep.seq = ifelse(is.sequential(ttids$cluster_id) == FALSE,TRUE,FALSE))
        
      }
    }
    
    if (subsize == 1 & length(checktypes) == 0) {
      temptypes <- as.character(unique(temp$cluster_label))
      print(paste0("Skipping ", temptypes, ", as it is also its own subclass"))
      next
    }
    
    
  }
  
  
  
}
