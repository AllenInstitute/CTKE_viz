#' Wrapper function to generate sub-umaps from larger one
#' 
#' @param select.cl
#' @param prefix
#' @param anno.df
#' @param rd.dat
#' @param rd.cl.means
#' @param sampled.cells
#' @param select.knn.cl.df
#' @param dest.d
#' @param do.init=TRUE
#' 
#' @export






find_umap <- function(select.cl, 
                      prefix, 
                      cl, 
                      anno.df, 
                      rd.dat, 
                      sampled.cells=NULL, 
                      dest.d, 
                      do.init=TRUE,
                      meta.fields=c("cluster"))
{
  
  ##
  if(!is.null(sampled.cells)) {
  tmp.sampled.cells= intersect(sampled.cells, names(cl[cl %in% select.cl]))
  } else {
    tmp.sampled.cells=names(cl[cl %in% select.cl])
  }
  
  tmp.sampled.cells <- intersect(rownames(rd.dat), tmp.sampled.cells)      
  sampled.cl <- cl[names(cl) %in% tmp.sampled.cells]
  tmp.rd.dat <- rd.dat[rownames(rd.dat) %in% tmp.sampled.cells,]
  
  rd.fn= paste0(prefix,".rd.dat.csv")
  data.table::fwrite(as.data.frame(tmp.rd.dat), file=file.path(dest.d,rd.fn), row.names=TRUE)
  umap.fn = paste0(prefix,".umap.2d.sampled.csv")  
  
  if(do.init){
    
    
    rd.cl.means = t(get_cl_means(tmp.rd.dat, sampled.cl))
    
    rd.cl.fn= paste0(prefix,".rd.dat.cl.means.csv")
    write.csv(rd.cl.means[rownames(rd.cl.means) %in% select.cl,], file=rd.cl.fn)
    
    umap.cl.fn = paste0(prefix,".umap.2d.cl.means.csv")  
    #cmd = paste("~/zizhen/bin/run_umap.py -i", rd.cl.fn, "-o", umap.cl.fn, "-n 10")
    cmd = paste("/allen/programs/celltypes/workgroups/mct-t200/Cindy_analysis/Product_release/2021_Cell_type_Cards_Miniatlas/run_umap.py -i", rd.cl.fn, "-o", umap.cl.fn, "-n 10")
    print(cmd)
    system(cmd)
    
    umap.cl.fn = paste0(prefix,".umap.2d.cl.means.csv")  
    cl.center.df <- read.csv(umap.cl.fn)
    colnames(cl.center.df)=c("cl","x","y")
    #cl.center.df$cl = as.character(cl.center.df$cl)
    umap.init = cl.center.df %>% left_join(anno.df) %>% select(sample_name, x, y)
    row.names(umap.init)= umap.init[[1]]
    umap.init = umap.init[,-1]
    umap.init.fn = paste0(prefix,".umap.init.csv")
    data.table::fwrite(umap.init[tmp.sampled.cells,], file= file.path(dest.d,umap.init.fn), row.names=TRUE)
    cmd = paste("/allen/programs/celltypes/workgroups/mct-t200/Cindy_analysis/Product_release/2021_Cell_type_Cards_Miniatlas/run_umap.py -i", rd.fn, "-o", umap.fn, "-t", umap.init.fn)
    }
  else{
    cmd = paste("/allen/programs/celltypes/workgroups/mct-t200/Cindy_analysis/Product_release/2021_Cell_type_Cards_Miniatlas/run_umap.py -i", rd.fn, "-o", umap.fn)
  }
  print(cmd)
  system(cmd)
  umap.fn = paste0(prefix,".umap.2d.sampled.csv")
  alpha=1
  if(length(tmp.sampled.cells) > 50000){
    alpha=0.6
  }
  umap.result = plot_2d_umap_anno(umap.fn, anno.df, dest.d, meta.fields=meta.fields)
  #plot_umap_constellation(umap.result$umap.2d, cl[sampled.cells], cl.df, select.knn.cl.df, dest.d=dest.d)    
  return(umap.result)
}














#' plot_2d_umap_anno
#' 
#' @param umap.fn path to umap coordinates. CSV file containing sample_names and umap x/y coordinates
#' @param anno.df Sample annotations. The first column should be sample_name, and each annotation should have \_id, \_label, and \_color columns. Requires cluster_id which needs to be sequential in order of the dendrogram.
#' @param dest.d path to save plots.
#' @param meta.fields  base name of variables to be represented as bargraphs below dendrogram. Annotation variables need to be represented as \_id, \_label, \_color in anno.
#' @param show.label TRUE or FALSE. To show cluster label on top of plot.
#' @param alpha level of transparency of dots. Value between 0 (fully transparent) and 1 (fully opaque)
#' @param cex size of plotted points. Default = 0.25
#' @param save.format figures can be saved as "png", "pdf" or "both"
#' @param plot.height 
#' @param plot.width
#' @param show.legend TRUE or FALSE. Whether to show legend with plot.  
#'   
#'    
#' @example_data:
#'  
#' load("data/rd_plot_example/example_umap.csv")
#' load("data/rd_plot_example/anno.df.rda")
#' 
#' 
#' @usage plots <- plot_2d_umap_anno(umap.fn="data/rd_plot_example/example_umap.csv",anno.df=anno.df, dest.d="./", meta.fields=c("platform","joint_region"),show.label=FALSE,alpha=0.5, cex=0.15,save.format="both", plot.height=7, plot.width=10, show.legend=TRUE)
#' 
#'  
#'    



plot_2d_umap_anno <- function(umap.fn, 
                              anno.df, 
                              dest.d="./",
                              meta.fields=NULL,
                              show.label=FALSE,
                              alpha=0.65, 
                              cex=0.25,
                              save.format=c("png","pdf","both"),
                              plot.height=7,
                              plot.width=7,
                              show.legend=FALSE)
{
  library(data.table)
  library(dplyr)
  library(ggplot2)
  
  #load umap from csv
  umap.df <- as.data.frame(fread(umap.fn,header=TRUE))
  colnames(umap.df) <- c("sample_name","Dim1","Dim2")
  umap.df <- umap.df[sample(1:nrow(umap.df)),]
  umap.df <- umap.df %>% left_join(anno.df) 
  umap.2d <- umap.df[,c("Dim1","Dim2")]
  row.names(umap.2d)<-umap.df$sample_name
  umap.2d <- umap.2d[sample(1:nrow(umap.2d)),]
  # extract filename for saving
  umap.fn <- basename(umap.fn)
  umap.fn <- gsub(".csv", "",umap.fn)
  
  #setup cluster labels/colors for plotting
  cl <- setNames(umap.df$cl, umap.df$sample_name)
  cl.df <- umap.df %>% select(cluster_id, cluster_label, cluster_color,cl) %>% unique
  cl.color <- setNames(cl.df$cluster_color, cl.df$cl)
  cl.label <- setNames(cl.df$cluster_label, cl.df$cl)
  
  plot.list <- list()
  #plot umap colored by cluster
  if(show.label==TRUE) {
    g <- plot_RD_cl(rd.dat=umap.2d, cl=cl, cl.color = cl.color, cl.label =cl.label,alpha.val=alpha, cex=cex, show.legend = FALSE)
    plot.list$cluster <- g 
  }  
  else {
    if(show.legend==TRUE){
      print("legend")
      g <- plot_RD_cl(rd.dat=umap.2d, cl=cl, cl.color = cl.color, cl.label =cl.label,alpha.val=alpha, cex=cex,label.center=FALSE, show.legend = TRUE)
      g[["labels"]][["colour"]] <- "Cluster"
      legend <- cowplot::get_legend(g)
      
      g <- plot_RD_cl(rd.dat=umap.2d, cl=cl, cl.color = cl.color, cl.label =cl.label,alpha.val=alpha, cex=cex,label.center=FALSE, show.legend = FALSE)
      g <- cowplot::plot_grid(g, legend, ncol=2)
      plot.list$cluster <- g 
    } 
    else{  
      g <- plot_RD_cl(rd.dat=umap.2d, cl=cl, cl.color = cl.color, cl.label =cl.label,alpha.val=alpha, cex=cex,label.center=FALSE, show.legend = FALSE)
      plot.list$cluster <- g 
    }
  }
  
  if(!is.null(meta.fields)) {
    # plot umap colored by other metadata
    for(m in meta.fields){
      tmp.df = umap.df[,paste0(m, c("_id","_label","_color"))] %>% unique
      colnames(tmp.df)=c("id","label","color")
      tmp.df = tmp.df %>% arrange(id)
      tmp.color = setNames(as.character(tmp.df$color), tmp.df$label)
      
      g= plot_RD_meta(rd.dat=umap.2d, 
                      meta=factor(umap.df[,paste0(m, "_label")], 
                             levels=names(tmp.color)),
                      meta.col = tmp.color,
                      alpha=alpha)
      

      if(show.legend==TRUE){
        print("legend")
        g[["labels"]][["colour"]] <- m
        legend <- cowplot::get_legend(g)   
        g = g + theme(axis.title.x=element_blank(), axis.title.y=element_blank())+ 
          theme_void() + 
          theme(legend.position="none") +
          coord_fixed(ratio=1)
        g <- cowplot::plot_grid(g, legend, ncol=2)
        plot.list[[m]] <- g
      }     
      else{
        print("no leg")
        g = g + theme(axis.title.x=element_blank(), axis.title.y=element_blank())+ 
          theme_void() + 
          theme(legend.position="none") +
          coord_fixed(ratio=1)
        plot.list[[m]] <- g
      }
    }
  }
  
  #save list of plots as pdf or png
  if(save.format == "pdf") {
    lapply(names(plot.list), function(nm)
      ggsave(plot=plot.list[[nm]], file=paste0(umap.fn,"_",nm, ".pdf"), useDingbats=FALSE, height=plot.height, width=plot.width   ))
  } else if(save.format == "png"){
    lapply(names(plot.list), function(nm)
      ggsave(plot=plot.list[[nm]], file=paste0(umap.fn,"_",nm, ".png"), height=plot.height, width=plot.width   ))
  } else if(save.format == "both"){
    lapply(names(plot.list), function(nm)
      ggsave(plot=plot.list[[nm]], file=paste0(umap.fn,"_",nm, ".pdf"), useDingbats=FALSE, height=plot.height, width=plot.width ))
    
    lapply(names(plot.list), function(nm)
      ggsave(plot=plot.list[[nm]], file=paste0(umap.fn,"_",nm, ".png"), height=plot.height, width=plot.width  ))
  }
  else{ print("Specify save.format")
  }
  
  return(plot.list)
}



plot_RD_cl <- function(rd.dat, cl, cl.color, cl.label,cex=0.15, fn.size =2, alpha.val=NULL,show.legend=FALSE, legend.size=2, label.center=TRUE, bg="blank",fn.color="black",no.shape=TRUE,ncol=2)
{
  rd.dat=as.data.frame(rd.dat)
  colnames(rd.dat) = paste0("Dim", 1:ncol(rd.dat))
  rd.dat$cl = factor(cl[row.names(rd.dat)])
  if(label.center){
    cl.center = get_RD_cl_center(rd.dat, cl)
  }
  if(!no.shape){
    shape = setNames(1:length(levels(rd.dat$cl)) %% 20 + 1,levels(rd.dat$cl))
    g=ggplot(rd.dat, aes(Dim1, Dim2)) + geom_point(aes(color=cl,shape=cl),size=cex)
    g = g+ scale_shape_manual(values=as.vector(shape[levels(rd.dat$cl)]))      
  }
  else{
    g=ggplot(rd.dat, aes(Dim1, Dim2)) + geom_point(aes(color=cl),size=cex)
  }
  if(!is.null(alpha.val)){
    col = alpha(as.vector(cl.color[levels(rd.dat$cl)]),alpha.val)
  }
  else{
    col = as.vector(cl.color[levels(rd.dat$cl)])
  }
  g = g+ scale_color_manual(values=col,labels=cl.label[levels(rd.dat$cl)])
  if(label.center){
    g = g + geom_point(data=as.data.frame(cl.center), aes(x=x, y=y), size=cex*1.5)
    for(i in 1:nrow(cl.center)){
      g = g +  annotate("text", label=cl.label[row.names(cl.center)[i]], x=cl.center[i,1], y=cl.center[i,2],size=fn.size,color=fn.color)
    }
  }
  if(bg=="blank"){
    g = g + theme_void()
    #g = g + theme(panel.background=element_blank())
    #g = g + theme(axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"))
  }
  else{
    g = g + theme(panel.background= element_rect(fill=bg, color=NA), panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor= element_blank())
  }
  if(show.legend){
    if(no.shape){
      g = g +  guides(colour = guide_legend(override.aes = list(size = legend.size),ncol=ncol))
    }
    else{
      g = g +  guides(colour = guide_legend(override.aes = list(shape = shape[levels(rd.dat$cl)],size = legend.size)),ncol=ncol)
    }
    g = g + theme(legend.position="right")
  }
  else{
    g = g + theme(legend.position="none")
  }
  g = g + coord_fixed(ratio=1)
  return(g)
}






###meta is discretized. 
#' Title
#'
#' @param rd.dat 
#' @param meta 
#' @param meta.col 
#' @param show.legend 
#' @param cex 
#' @param legend.size 
#' @param alpha.val 
#'
#' @return
#' @export
#'
#' @examples
plot_RD_meta <- function(rd.dat, meta, meta.col=NULL,show.legend=TRUE, cex=0.15, legend.size=5,alpha.val=1)
{
  rd.dat = as.data.frame(rd.dat)
  colnames(rd.dat)[1:2] = c("Dim1","Dim2")
  library(ggplot2)
  rd.dat$meta = meta
  p=ggplot(rd.dat, aes(Dim1, Dim2)) + geom_point(aes(color=meta),size=cex)
  if(is.factor(meta)){
    rd.dat = droplevels(rd.dat)
    if(is.null(meta.col)){
      if(length(levels(meta)) > 2){
        meta.col = setNames(jet.colors(length(levels(meta))), levels(meta))
      }
      else{
        meta.col = setNames(c("blue", "orange"), levels(meta))
      }
    }      
    p = p+ scale_color_manual(values=alpha(as.vector(meta.col[levels(rd.dat$meta)]),alpha.val))
    p = p+ theme(panel.background=element_blank(),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"))
  }
  else{
    p = p+ scale_color_gradient(low="blue",high="red")
  }
  if(!show.legend){
    p = p + theme(legend.position="none") 
  }
  else{
    if(is.factor(meta)){
      p = p + guides(colour = guide_legend(override.aes = list(size=legend.size)))
    }
  }
  p = p + coord_fixed(ratio=1)
  return(p)
}






