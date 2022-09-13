library(scrattch.hicat)
library(dendextend)
library(jsonlite)
library(scrattch.io)
library(dplyr)
library(ggplot2)

in.dir <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/Miniatlas/transcriptome_joint"
# 
 load(file.path(in.dir, "cl.final.rda"))
# 
# View(cl.df)
# #write.csv(cl.df, file="cl.df.csv")
# 
# #add neighborhoods
# 
# cl.df <- read.csv("cl.df.csv")
# cl.df <- cl.df[,-1]
# 
# save(cl, cl.df, file="cl.df.rda")

cl.df <- read.csv("cell_type_to_NH.csv", row.names=1)


#########################
### Dendrogram
load(file.path(in.dir,"anno.df.rda"))
anno <- anno.df[anno.df$class_label != "Low Quality",]


load(file.path(in.dir,"dend.labeled.rda"))

ggsave("Plots/dend.nr.pdf",plot_dend(dend.labeled), height=5,width=10)
ggsave("Plots/dend.pdf",plot(dend.labeled), height=5,width=10)
ggsave("Plots/Global_dend.png",plot(dend.labeled), height=5,width=10)

clades <- unique(cl.df$neighborhood_label)
for(i in clades) {
  
  rm.labels <- cl.df$cluster_label[cl.df$neighborhood_label != i ]
  n.dend <- prune_dend(dend = dend.labeled,
             rm.labels = rm.labels,
             top.level = TRUE)  
  ggsave(paste0("Plots/",i,"_dend.pdf"),plot(n.dend), height=5,width=10)
  ggsave(paste0("Plots/",i,"_dend.png"),plot(n.dend), height=5,width=10)
  
}


## adding new node labels  
nomenclature_information <- build_nomenclature_table( dend.labeled)			


allMatrix <- read.csv("node.info.csv")
dend=dend.labeled

# Convert "label" to "original_label"

dendro_data = as.ggdend(dend)
dendro_data$nodes$label = get_nodes_attr(dend, "label")
dendro_data$nodes = dendro_data$nodes[is.na(dendro_data$nodes$leaf),]

node_data = dendro_data$nodes
segment_data <- dendro_data$segments
label_data <- dendro_data$labels

node_data$new.lab <- allMatrix$label[match(node_data$label, allMatrix$original_label)]
node_data$new.lab[node_data$new.lab == node_data$label] <- ""

node_size = 2; r = c(-0.1, 1)
node_data$node_color = "black"
p=  ggplot() + 
    geom_text(data = node_data, aes(x = x, y = y,label = new.lab, color = node_color), size = node_size, vjust = 1) + 
    geom_segment(data = segment_data, aes(x = x,xend = xend, y = y, yend = yend), color = "gray50") + 
    geom_text(data = label_data, aes(x = x, y = -0.01, label = label, color = col), size = node_size, angle = 90, hjust = 1) + 
  scale_color_identity() + theme_dendro() + scale_y_continuous(limits = r)

ggsave("Plots/dend.labs.pdf",p, height=5,width=10)
ggsave("Plots/dend.labs.png",p, height=5,width=10)


#########################
### Sunburst
source("//allen/programs/celltypes/workgroups/mct-t200/Cindy_analysis/my_Tools/R_utils/Plotting/Plot_sunburst.R")


library(dplyr)
library(plotly)

sunburstDF <- as.sunburstDF(cl.df,
                            levels=c("class","neighborhood","subclass","cluster"),
                            valueCol = "cluster_size",
                            rootname="MOp")


p <- plot_ly() %>%
      add_trace(ids = sunburstDF$ids,
                labels = sunburstDF$labels,
                parents =sunburstDF$parent,
                values = sunburstDF$values,
                type = 'sunburst',
                marker = list(colors = sunburstDF$color),
                domain = list(column = 1),
                branchvalues = 'total'
                )%>%
      layout(grid = list(columns =1, rows = 1),
              margin = list(l = 0, r = 0, b = 0, t = 0)
      )

p


library(htmlwidgets)
saveWidget(p, "sunburst_cluster.html", selfcontained = F, libdir = "lib_burst")
saveWidget(p, "sunburst_cell.html", selfcontained = F, libdir = "lib_burst")





#########################
### UMAP

### Global UMAP

## files used in manuscript:
#umap.2d.md0.3.sampled.pdf
#rd.dat.sampled.csv
in.dir <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/Miniatlas/transcriptome_joint"

umap.2d <- read.csv(file.path(in.dir,"umap.2d.md0.3.sampled.csv"), row.names=1)




### Neighborhood UMAPs

source('F_find_UMAP.R')

load("cl.df.rda")

load(file.path(in.dir,"anno.df.rda"))
load(file.path(in.dir,"sampled.cells.rda"))
rd.dat <- read.csv(file.path(in.dir,"rd.dat.csv"), row.names=1)

rd.dat <- as.matrix(rd.dat)

clades <- unique(cl.df$neighborhood_label)

#n="IT-PT"
for(n in clades){
  print(n)
  
  select.cl = cl.df$cl[cl.df$neighborhood_label==n]
  
  umap.2d <- find_umap(select.cl=select.cl, 
                       prefix=n, 
                       cl=cl, 
                       anno.df=anno.df, 
                       rd.dat=rd.dat, 
                       sampled.cells= NULL,
                       dest.d="./",
                       do.init=TRUE,
                       meta.fields=NULL)
  }


### Neighborhood UMAPs starting from global umap centroids

cl.center.umap.2d <- as.data.frame(get_RD_cl_center(umap.2d, cl))

n=clades[3]

for(n in clades){
  
  prefix=n
  print(prefix)
  
  select.cl = cl.df$cl[cl.df$neighborhood_label==n]
  select.cl.center <- cl.center.umap.2d[rownames(cl.center.umap.2d) %in% select.cl,]
  select.cl.center$cl <- rownames(select.cl.center) 
  
  tmp.sampled.cells <- names(cl[cl %in% select.cl])
  
  select.cl.center$cl = as.numeric(select.cl.center$cl)
  umap.init = select.cl.center %>% left_join(anno.df) %>% select(sample_name, x, y)
  
  row.names(umap.init)= umap.init[[1]]
  umap.init = umap.init[,-1]
  umap.init.fn = paste0(prefix,".umap.init.csv")
  data.table::fwrite(umap.init, file= umap.init.fn, row.names=TRUE)
  
  tmp.rd.dat <- rd.dat[rownames(umap.init),]
  
  rd.fn= paste0(prefix,".rd.dat.csv")
  data.table::fwrite(as.data.frame(tmp.rd.dat), file=rd.fn, row.names=TRUE)
  umap.fn = paste0(prefix,".umap.2d.sampled.csv")  
  
  
  cmd = paste("run_umap.py -i", rd.fn, "-o", umap.fn, "-t", umap.init.fn)

  print(cmd)
  system(cmd)
  
  plot_2d_umap_anno(umap.fn, 
                    anno.df, 
                    dest.d="./",
                    meta.fields=NULL,
                    show.label=FALSE,
                    alpha=0.5, 
                    cex=0.25,
                    save.format="both",
                    plot.height=7,
                    plot.width=7,
                    show.legend=TRUE)
}

#########################
### UMAP plotting std

neighborhoods <- unique(cl.df$neighborhood_label)
in.dir= paste0(getwd(),"/Plots")
out.dir= paste0(getwd(),"/Neighborhood_plots")

for(hood in neighborhoods) {
  print(hood)
  umap.2d <- read.csv(file.path(in.dir,paste0(hood,".umap.2d.sampled.csv")))
  colnames(umap.2d) = c("sample_name","Dim1","Dim2")
  row.names(umap.2d) = umap.2d[[1]]
  umap.df = umap.2d %>% left_join(anno.df)
  umap.2d = umap.2d[-1]
  # umap with all cluster colors
  clus.cl = with(umap.df, setNames(cluster_label,sample_name))
  clus.col = setNames(cl.df$cluster_color, cl.df$cluster_label)
  g= plot_RD_meta(umap.2d, meta = factor(clus.cl), meta.col = clus.col,alpha=0.5)
  g = g + coord_fixed(ratio=1)
  g = g + theme_void()
  g = g + theme(legend.position="none")
  ggsave(filename = file.path(out.dir,paste0(hood,"_UMAP2d_cluster.png")),
         plot = g,
         width = 4,
         height = 4,
         units = c("in"),
         dpi = 300)
  ##### plot umap with 1 cluster highlighted
  clust.df <- unique(umap.df[,c("cluster_id", "cluster_label", "cluster_color")])
  rownames(clust.df) <- 1:nrow(clust.df)
  for(i in 1:nrow(clust.df)) {
    prefix =paste(clust.df$cluster_id[i], clust.df$cluster_label[i])
    print(prefix)
    prefix <- sub(" ", "_", prefix)
    prefix <- sub("/","-", prefix)
    clust.df$cl.col <- ifelse(clust.df$cluster_color == clust.df$cluster_color[i], as.character(clust.df$cluster_color), "grey80")
    clus.col = setNames(clust.df$cl.col, clust.df$cluster_label)
    g= plot_RD_meta(umap.2d, meta = factor(clus.cl), meta.col = clus.col,alpha=0.5)
    g = g + coord_fixed(ratio=1)
    g = g + theme_void()
    g = g + theme(legend.position="none")
    ggsave(filename = file.path(out.dir,paste0(prefix,"_",hood,"_UMAP2d_cluster.png")),
           plot = g,
           width = 4,
           height = 4,
           units = c("in"),
           dpi = 300)
  }
}

#########################
### Constellation


load("cl.df.rda")
# in.dir <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/Miniatlas/transcriptome_joint"
# rd.dat <- read.csv(file.path(in.dir,"rd.dat.csv"), row.names=1)
# 
# rd.dat <- as.matrix(rd.dat)
# 
# knn.result <- get_knn_graph(rd.dat, cl=cl, cl.df=cl.df, k=15, knn.outlier.th=2, outlier.frac.th=0.5)
# save(knn.result, file="knn.result.rda")

load("knn.result.rda")

### Global 

knn.cl.df <- knn.result[["knn.cl.df"]]

umap.2d <- read.csv(file.path(in.dir,"umap.2d.md0.3.sampled.csv"), row.names=1)

colnames(umap.2d) = c("Lim1","Lim2")

select.cl <- cl[names(cl) %in% rownames(umap.2d)]
cl.center.umap.2d <- get_RD_cl_center(umap.2d, select.cl)
cl.center.umap.2d <- as.data.frame(cl.center.umap.2d)
cl.center.umap.2d <- tibble::rownames_to_column(cl.center.umap.2d, "cl")
cl.center.umap.2d$cl <- as.numeric(as.character(cl.center.umap.2d$cl))
cl.center.df <- left_join(cl.center.umap.2d, select(cl.df, cl, cluster_id, cluster_label, cluster_color, size), by="cl") 
names(cl.center.df)[names(cl.center.df) == 'size'] <- 'cluster_size'

knn.cl.df <- knn.cl.df[knn.cl.df$Freq >=3,]

plot_constellation(knn.cl.df = knn.cl.df, cl.center.df = cl.center.df, out.dir = "Plots", plot.parts = FALSE, node.dodge=TRUE, exxageration=2 )


### Neighborhood constellation

clades <- unique(cl.df$neighborhood_label)

#additional 
cl.df <- cl.df[,-14]

for(n in clades){
  
    print(n)

    umap.fn = paste0(n,".umap.2d.sampled.csv")
    umap.2d <- read.csv(file.path("Plots",umap.fn),row.names=1)
    colnames(umap.2d) = c("Lim1","Lim2")
    
    print("UMAP centroids")
    cl.center.umap.2d <- get_RD_cl_center(umap.2d, cl)
    
    cl.center.umap.2d <- as.data.frame(cl.center.umap.2d)
    cl.center.umap.2d <- tibble::rownames_to_column(cl.center.umap.2d, "cl")
    cl.center.umap.2d$cl <- as.integer(as.character(cl.center.umap.2d$cl))
    cl.center.df <- left_join(cl.center.umap.2d, cl.df, by="cl") 
    #cl.center.df <- left_join(cl.center.umap.2d, select(cl.df, cl, cluster_id, cluster_label, cluster_color, cluster_size, clade_id, clade_label, clade_color), by="cl") 
    save(cl.center.df, file=file.path("Plots", paste0(n,".cl.center.df.rda")))
    #write.csv(cl.center.umap.2d, file=paste0(prefix, ".cl.center.umap.2d.csv"))
    print("plotting constellation")
    plot_constellation(knn.cl.df = knn.cl.df, 
                       cl.center.df = cl.center.df, 
                       out.dir = "Plots", 
                       node.label="cluster_id", 
                       exxageration=4, 
                       curved = TRUE, 
                       plot.parts=FALSE, 
                       plot.height=25, 
                       plot.width=25, 
                       node.dodge=TRUE, 
                       label.size=2, 
                       max_size=10 )
    }
  




#########################
### Gene expression dotplot
library(tidyverse)
library(scrattch.vis)

source("molgen-shiny/distillery/annocomp_functions.R")

mouse_meta <- as.data.frame(feather::read_feather( "Mouse_M1_10xV3_Metadata.feather"))
mouse_data <- readRDS("Mouse_M1_10xV3_Matrix.RDS")



CGE.genes <- c("Gad1, Lamp5, Prox1, Ndnf, Sncg, Vip, Lhx6, Pax6, Chat, Egln3, Serpinf1")
tmp <-stringi::stri_extract_all_words(CGE.genes)
CGE.genes <- unlist(tmp)
#CGE.genes[6] <- "Nkx2-1"


MGE.genes <- c("Gad1","Lhx6, Sst, Chodl, Nos1, Myh8, Crh, C1ql3, Reln, Th, Calb1 Pvalb, Vipr2")
tmp <-stringi::stri_extract_all_words(MGE.genes)
MGE.genes <- unlist(tmp)

IT_ET.genes <- c("Slc17a7, Slc30a3, Stard8, Otof, Rorb, Rspo1, Ddit4l, S100b, Sulf1, Rxfp1, Fam84b,Bcl6")
tmp <-stringi::stri_extract_all_words(IT_ET.genes)
IT_ET.genes <- unlist(tmp)

NP_CT_L6b.genes <- c("Slc17a7, Fezf2, Sla2, Grik1, Nxph1, Sulf1, Syt6, Foxp2, Nxph4, Cplx3")
tmp <-stringi::stri_extract_all_words(NP_CT_L6b.genes)
NP_CT_L6b.genes <- unlist(tmp)

Other.genes <- c("Meis2, Mbp, Mog, Pdgfra, Gfap, Csf1r, Slc6a13, Vwf")
tmp <-stringi::stri_extract_all_words(Other.genes)
Other.genes <- unlist(tmp)


marker.genes <- list(CGE=CGE.genes,
                     MGE=MGE.genes,
                     IT_ET=IT_ET.genes,
                     NP_CT_L6b=NP_CT_L6b.genes, 
                     Other=Other.genes)

clades <- unique(cl.df$neighborhood_label)

anno <- anno.df %>% filter(platform_label == "snRNA 10X v3 B")
anno.cl <- cl[names(cl) %in% anno$sample_name] 
anno <- as.data.frame(anno.cl)
colnames(anno) <- "cl"
anno$sample_name <- sub("10X_nuclei_v3_Broad.", "", rownames(anno))	
anno <- anno[,c(2,1)]                      
anno <- anno %>% left_join(cl.df)
#mouse_data <- cpm(mouse_data)

for(n in clades){
  
  print(n)
  
  tmp.anno = anno[anno$neighborhood_label==n,]
  
  prefix <- sub("-", "_", n)
  prefix <- sub("-", "_", prefix)
  genes_use = marker.genes[[prefix]]
  
  tmp.dat <- mouse_data[rownames(mouse_data) %in% genes_use,colnames(mouse_data) %in% tmp.anno$sample_name ]
  tmp.dat <- as.data.frame(t(as.matrix(tmp.dat)))
  tmp.dat$sample_name <- rownames(tmp.dat)
  
  print("plotting")
  g=group_dot_plot(tmp.dat, 
                   tmp.anno, 
                   genes=genes_use, 
                   grouping="cluster", 
                   log_scale=TRUE, 
                   normalize_rows=TRUE, 
                   fill_stat="mean",
                   label_height=25, 
                   show_counts=FALSE, 
                   colorset=c("dodgerblue","gray80","orange","orangered"))

  ggsave(file=paste0(prefix,"_dotplot_geneexpression.png"), 
         plot = g,
         width = 4,
         height = 6,
         units = c("in"),
         dpi = 300)
  
}




























