
library(scrattch.hicat)
library(dendextend)
library(jsonlite)
library(scrattch.io)
library(dplyr)
library(ggplot2)

in.dir <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/Miniatlas/transcriptome_joint"
# 

load("cl.df.rda")
load(file.path(in.dir,"anno.df.rda"))

### add platform to umap

#umap.2d <- read.csv(file.path(in.dir,"umap.2d.md0.3.sampled.csv"), row.names=1)

umap.2d <- read.csv(file.path(in.dir,"umap.2d.all.csv"))
colnames(umap.2d) = c("sample_name","Dim1","Dim2")

#library(stringr)
#umap.2d$platform <- word(umap.2d$sample_name,1,sep = "\\.")
row.names(umap.2d) = umap.2d[[1]]
umap.df = umap.2d %>% left_join(anno.df)
umap.2d = umap.2d[-1]


source('F_find_UMAP.R')



##########################################
### UMAP plotting txn levels by platform
##########################################

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

out.dir <- "umap_by_platform"

#highlight cells per class per platform
cluster_lab <- clust.df$cluster_label#[59:116]
class = class#[2]
subclass = subclass#[4]
#for(zz  in class){
for(zz  in subclass){
#for(zz in cluster_lab[59:116]){    
  print(zz)
  ### select points to plot in grey
  #bg.umap <- umap.df[umap.df$class_label != zz,]
  bg.umap <- umap.df[umap.df$subclass_label != zz,]
  #bg.umap <- umap.df[umap.df$cluster_label != zz,]
  
  bg.umap$color <- "grey80"
  bg.umap$alpha.val <- 0.3
  ###select points to plot in color
  #sub.umap <- umap.df[umap.df$class_label == zz,]
  sub.umap <- umap.df[umap.df$subclass_label == zz,]
  #sub.umap <- umap.df[umap.df$cluster_label == zz,]
  
  for(p in 1:nrow(platform)) {
    sel.platform <- platform$platform_long_label[p]
    print(sel.platform)
    prefix =paste0(zz,"_",sel.platform)
    
    bg.umap.2 <- sub.umap[sub.umap$platform_long_label != sel.platform,]
    if(nrow(bg.umap.2)>0) {
      bg.umap.2$color <- "grey80"
      bg.umap.2$alpha.val <- 0.3
    }
    
    fg.umap <- sub.umap[sub.umap$platform_long_label == sel.platform,]
    
    if(nrow(fg.umap)>0){
    #fg.umap$color <- fg.umap$class_color
    fg.umap$color <- fg.umap$subclass_color
    #fg.umap$color <- fg.umap$cluster_color
    fg.umap$alpha.val <- 1
    
      if(nrow(bg.umap.2)>0) {
    plot.umap <- rbind(bg.umap, bg.umap.2, fg.umap)
    }       else {plot.umap <- rbind(bg.umap, fg.umap)}
    
    g = ggplot(plot.umap, aes(Dim1, Dim2)) + 
      geom_point(colour=plot.umap$color, size = 0.15, alpha= plot.umap$alpha.val)
    g = g + theme(panel.background = element_blank(), 
                  axis.line.x = element_line(colour = "black"), 
                  axis.line.y = element_line(colour = "black"))
    g = g + coord_fixed(ratio=1)
    g = g + theme_void()
    g = g + theme(legend.position="none")
    #g = g + ggtitle(prefix)
    
    prefix <- gsub("/","-", prefix)
    
    ggsave(filename = file.path(out.dir,paste0(prefix,"_UMAP2d.png")),
           plot = g,
           width = 4,
           height = 4,
           units = c("in"),
           dpi = 300)
    
    ggsave(filename = file.path(out.dir,paste0(prefix,"_UMAP2d.svg")),
           plot = g,
           width = 4,
           height = 4,
           units = c("in"),
           dpi = 300)
    } else { text = paste("\n   No cells")
        g=ggplot() + 
            annotate("text", x = 4, y = 25, size=8, label = text) + 
            theme_void()
        prefix <- gsub("/","-", prefix)
        ggsave(filename = file.path(out.dir,paste0(prefix,"_UMAP2d.png")),
               plot = g,
               width = 4,
               height = 4,
               units = c("in"),
               dpi = 300)}
    
    
  }
}
    

##########################################
### UMAP plotting gene expression by platform
##########################################

library(ggnewscale)

dat.dir <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT8.0_Miniatlas"


platform <- unique(umap.df[,c("platform_id","platform_color","platform_label", "platform_long_label")])
platform <- taRifx::remove.factors(platform)
platform <- platform[order(platform$platform_id),]
rownames(platform) <- platform$platform_id


markers <- read.csv("yao_markers.csv")

genes <- markers %>% tidyr::pivot_longer(cols=c(GABAergic, Glutamatergic, Non.Neuronal.Non.Neural),names_to = "class", values_to = "gene")
genes <- genes[genes$gene != "",]
genes <- genes$gene
genes <- unique(genes)
out.dir <- "umap_Yao.markers_by_platform"



markers <- read.csv("NSForestMarkers_genesymbols.csv")
genes <- markers %>% tidyr::pivot_longer(cols=c(M1_gene, M2_gene, M3_gene, M4_gene, M5_gene),names_to = "class", values_to = "gene", values_drop_na = TRUE)
genes <- genes %>% select(gene, class)
#genes <- genes[genes$gene != "",]
#genes <- genes %>% filter_all(any_vars(!is.na(.)))
genes <- genes$gene
genes <- unique(genes)
out.dir <- "umap_NSF.markers_by_platform"



library(foreach)  # for parallelization
library(doMC)
registerDoMC(25)



for(p in 2:nrow(platform)) {  
  sel.platform <- platform$platform_long_label[p]
  print(sel.platform)
  
  bg.umap <- umap.df[umap.df$platform_long_label != sel.platform,]
  bg.umap$color <- "grey90"
  bg.umap$alpha.val <- 0.3
  
  fg.umap <- umap.df[umap.df$platform_long_label == sel.platform,]
  
  print("loading data")
  load(file.path(dat.dir, paste0(sel.platform,"/norm.dat.rda")))
  
  genes <- intersect(genes, rownames(norm.dat))
  
  colnames(norm.dat) <- paste0(sel.platform,".", colnames(norm.dat))
  norm.dat <- norm.dat[genes,colnames(norm.dat) %in% fg.umap$sample_name]
  
  norm.dat <- as.matrix(norm.dat)
  #norm.dat <- as.data.frame(t(2^norm.dat-1))
  norm.dat <- as.data.frame(t(norm.dat))
  norm.dat$sample_name <- rownames(norm.dat)
  
  fg.umap <- left_join(fg.umap, norm.dat)
  
  print("start plotting genes")
  foreach(i=1:length(genes)) %dopar% {
    #for(g  in genes){
    library(ggplot2)
    g=genes[i]
    print(g)
    
    prefix =paste0(g,"_",sel.platform)
    
    
    pic = ggplot() + 
      geom_point(data=bg.umap, 
                 aes(x=Dim1, y=Dim2,colour="grey85"),
                 size = 0.15,
                 alpha=0.3) +
      scale_color_identity() +
      new_scale_color() +
      geom_point(data=fg.umap, 
                 aes(x=Dim1, y=Dim2,colour=.data[[g]]), 
                 size = 0.15,
                 alpha=1) +
      scale_color_gradient(low = "gray70", high = "red")
    
    pic = pic + theme(panel.background = element_blank(), 
                      axis.line.x = element_line(colour = "black"), 
                      axis.line.y = element_line(colour = "black"))
    pic = pic + coord_fixed(ratio=1)
    pic = pic + theme_void()
    pic = pic + theme(legend.position="none")
    
    
    ggsave(filename = file.path(out.dir,paste0(prefix,"_UMAP2d.png")),
           plot = pic,
           width = 4,
           height = 4,
           units = c("in"),
           dpi = 300)
    
    ggsave(filename = file.path(out.dir,paste0(prefix,"_UMAP2d.svg")),
           plot = pic,
           width = 4,
           height = 4,
           units = c("in"),
           dpi = 300)
  } 
  
  
}


parallel::stopCluster(cl = cl)

  





