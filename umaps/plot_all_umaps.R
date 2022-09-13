#' This function is an adaption of Cindy's code for plotting UMAPs for cell type cards. It assumes UMAP coordinates
#' have already been calculated and requires these coordinates and a standard annotation file to work. It creates PNG
#' files of UMAPs for every taxonomy node in an annotation file.
#' 
#' @param umap.coords - dataframe of 2D UMAP coordinates for a taxonomy
#' @param anno - dataframe of annotations in standard format (requires cluster, subclass and class)
#' @param dataset_name - This will be included in the output file name
#' @param outdir - location to save files to


plot_all_umaps <- function(umap.coords,
                           anno,
                           dataset_name,
                           outdir){

  suppressPackageStartupMessages({
    library(scrattch.hicat)
    library(dendextend)
    library(jsonlite)
    library(scrattch.io)
    library(dplyr)
    library(ggplot2)})


  ### Load UMAP coordinates and annotations

  umap.2d <- umap.coords
  anno.df <- anno
  colnames(anno.df)[1] <- "sample_name"
  umap.2d <- umap.2d %>% select(sample_id, UMAP_1, UMAP_2)
  colnames(umap.2d) = c("sample_name","Dim1","Dim2")

  row.names(umap.2d) = umap.2d[[1]]
  umap.df = umap.2d %>% left_join(anno.df)
  umap.2d = umap.2d[-1]
  umap.df$platform_id <- rep(1,nrow(umap.df))
  umap.df$platform_color <- rep("#000001",nrow(umap.df))
  umap.df$platform_label <- rep(dataset_name,nrow(umap.df))
  umap.df$platform_long_label <- rep("xxx")






  ##########################################
  ### UMAP plotting txn levels by platform
  ##########################################
  setwd(outdir)
  class <- unique(umap.df$class_label)
  subclass <- unique(umap.df$subclass_label)

  platform <- unique(umap.df[,c("platform_id","platform_color","platform_label", "platform_long_label")])
  platform <- taRifx::remove.factors(platform)
  platform <- platform[order(platform$platform_id),]
  rownames(platform) <- platform$platform_id

  clust.df <- unique(umap.df[,c("cluster_id", "cluster_label", "cluster_color")])
  clust.df <- clust.df[order(clust.df$cluster_id),]
  clust.df$background <- "grey80"
  rownames(clust.df) <- clust.df$cluster_id

  out.dir <- out.dir

  #highlight cells per class per platform
  cluster_lab <- clust.df$cluster_label
  class <- class
  subclass <- subclass
  
  for(zz  in subclass){
    
    print(zz)
    ### select points to plot in grey
    bg.umap <- umap.df[umap.df$subclass_label != zz,]
  
    bg.umap$color <- "grey80"
    bg.umap$alpha.val <- 0.3
    ###select points to plot in color
    sub.umap <- umap.df[umap.df$subclass_label == zz,]
  
    for(p in 1:nrow(platform)) {
      sel.platform <- platform$platform_label[p]
      print(sel.platform)
      prefix =paste0(zz,"_",sel.platform)
    
      
      bg.umap.2 <- sub.umap[sub.umap$platform_label != sel.platform,]
      if(nrow(bg.umap.2)>0) {
        bg.umap.2$color <- "grey80"
        bg.umap.2$alpha.val <- 0.3
      }
    
      fg.umap <- sub.umap[sub.umap$platform_label == sel.platform,]
    
      if(nrow(fg.umap)>0){
        fg.umap$color <- fg.umap$subclass_color
        fg.umap$alpha.val <- 1
      
        if(nrow(bg.umap.2)>0) {
          plot.umap <- rbind(bg.umap, bg.umap.2, fg.umap)
        }       
        
        else {plot.umap <- rbind(bg.umap, fg.umap)}
      
        g <- ggplot(plot.umap, aes(Dim1, Dim2)) + 
          geom_point(colour=plot.umap$color, size = 0.15, alpha= plot.umap$alpha.val)
        g <- g + theme(panel.background = element_blank(), 
                    axis.line.x = element_line(colour = "black"), 
                    axis.line.y = element_line(colour = "black"))
        g <- g + coord_fixed(ratio=1)
        g <- g + theme_void()
        g <- g + theme(legend.position="none")
        #g = g + ggtitle(prefix)
      
        prefix <- gsub("/","-", prefix)
      
        ggsave(filename = paste0(prefix,"_UMAP2d.png"),
               plot = g,
               width = 4,
               height = 4,
               units = c("in"),
               dpi = 300)
      } 
      
      else { text = paste("\n   No cells")
      g <- ggplot() + 
        annotate("text", x = 4, y = 25, size=8, label = text) + 
        theme_void()
      prefix <- gsub("/","-", prefix)
      ggsave(filename = paste0(prefix,"_UMAP2d.png"),
             plot = g,
             width = 4,
             height = 4,
             units = c("in"),
             dpi = 300)}
    
    
    }
  }



  ######################
  #### CLASS UMAPS #####
  ######################


  for(zz  in class){  
    print(zz)
    ### select points to plot in grey
    bg.umap <- umap.df[umap.df$class_label != zz,]
  
    bg.umap$color <- "grey80"
    bg.umap$alpha.val <- 0.3
    ###select points to plot in color
    sub.umap <- umap.df[umap.df$class_label == zz,]
  
    for(p in 1:nrow(platform)) {
      sel.platform <- platform$platform_label[p]
      print(sel.platform)
      prefix <- paste0(zz,"_",sel.platform)
    
      bg.umap.2 <- sub.umap[sub.umap$platform_label != sel.platform,]
      if(nrow(bg.umap.2)>0) {
        bg.umap.2$color <- "grey80"
        bg.umap.2$alpha.val <- 0.3
      }
    
      fg.umap <- sub.umap[sub.umap$platform_label == sel.platform,]
      
      if(nrow(fg.umap)>0){
        fg.umap$color <- fg.umap$class_color
        fg.umap$alpha.val <- 1
      
        if(nrow(bg.umap.2)>0) {
          plot.umap <- rbind(bg.umap, bg.umap.2, fg.umap)
        }       else {plot.umap <- rbind(bg.umap, fg.umap)}
      
        g <- ggplot(plot.umap, aes(Dim1, Dim2)) + 
          geom_point(colour=plot.umap$color, size = 0.15, alpha= plot.umap$alpha.val)
        g <- g + theme(panel.background = element_blank(), 
                    axis.line.x = element_line(colour = "black"), 
                    axis.line.y = element_line(colour = "black"))
        g <- g + coord_fixed(ratio=1)
        g <- g + theme_void()
        g <- g + theme(legend.position="none")
        #g = g + ggtitle(prefix)
      
        prefix <- gsub("/","-", prefix)
      
        ggsave(filename = paste0(prefix,"_UMAP2d.png"),
               plot = g,
               width = 4,
               height = 4,
               units = c("in"),
               dpi = 300)
      } else { text = paste("\n   No cells")
      g <- ggplot() + 
        annotate("text", x = 4, y = 25, size=8, label = text) + 
        theme_void()
      prefix <- gsub("/","-", prefix)
      ggsave(filename = paste0(prefix,"_UMAP2d.png"),
             plot = g,
             width = 4,
             height = 4,
             units = c("in"),
             dpi = 300)}
    
    
    }
  }


  ######################
  #### CLUSTER UMAPS #####
  ######################


  for(zz  in cluster_lab){ 
    print(zz)
    ### select points to plot in grey
    bg.umap <- umap.df[umap.df$cluster_label != zz,]
  
    bg.umap$color <- "grey80"
    bg.umap$alpha.val <- 0.3
    ###select points to plot in color
    sub.umap <- umap.df[umap.df$cluster_label == zz,]
  
    for(p in 1:nrow(platform)) {
      sel.platform <- platform$platform_label[p]
      print(sel.platform)
      prefix =paste0(zz,"_",sel.platform)
    
      bg.umap.2 <- sub.umap[sub.umap$platform_label != sel.platform,]
        if(nrow(bg.umap.2)>0) {
          bg.umap.2$color <- "grey80"
          bg.umap.2$alpha.val <- 0.3
        }
    
      fg.umap <- sub.umap[sub.umap$platform_label == sel.platform,]
    
      if(nrow(fg.umap)>0){
        fg.umap$color <- fg.umap$cluster_color
        fg.umap$alpha.val <- 1
      
        if(nrow(bg.umap.2)>0) {
          plot.umap <- rbind(bg.umap, bg.umap.2, fg.umap)
        }       else {plot.umap <- rbind(bg.umap, fg.umap)}
      
        g <- ggplot(plot.umap, aes(Dim1, Dim2)) + 
          geom_point(colour=plot.umap$color, size = 0.15, alpha= plot.umap$alpha.val)
        g <- g + theme(panel.background = element_blank(), 
                      axis.line.x = element_line(colour = "black"), 
                      axis.line.y = element_line(colour = "black"))
        g <- g + coord_fixed(ratio=1)
        g <- g + theme_void()
        g <- g + theme(legend.position="none")
        #g = g + ggtitle(prefix)
      
        prefix <- gsub("/","-", prefix)
      
        ggsave(filename = paste0(prefix,"_UMAP2d.png"),
               plot = g,
               width = 4,
               height = 4,
               units = c("in"),
               dpi = 300)
      } else { text = paste("\n   No cells")
      g <- ggplot() + 
        annotate("text", x = 4, y = 25, size=8, label = text) + 
        theme_void()
      prefix <- gsub("/","-", prefix)
      ggsave(filename = paste0(prefix,"_UMAP2d.png"),
             plot = g,
             width = 4,
             height = 4,
             units = c("in"),
             dpi = 300)}
    
    
    }
  }
}
