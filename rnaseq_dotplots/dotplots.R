# Script to make dotplots for all taxons. This script plots the Broad data, but can modified by
# simply loading in another dataset (done in other scripts not included here but in my //allen directory). The basic steps are:
#
# - 1. Load in annotations and expression matrix
# - 2. Loop through all taxons at each level of tree (class, subclass and cluster). Again, maybe not the most elegant 
#      solution here, but it works for now.
# - 3. Plot and write to svg/csv

library(plyr)
library(tidyverse)
library(Matrix)
setwd("//allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr/")
source("load_metadata.R")
setwd("//allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr/snRNA 10xV3 Broad")
broadcounts <- readRDS("Mouse_M1_10xV3_Matrix.RDS")

classgenes <- genelist$GABAergic


plotcounts <- cbind(sample_name = colnames(as.matrix(broadcounts)),
                    as.data.frame(t(as.matrix(broadcounts[classgenes,]))))
plotcounts <- plotcounts[plotcounts$sample_name %in% snRNA_10xv3_B_met$sample_name,]


# Colorset for dot plot
colorset <- c("#fee8c8","#fdbb84","#e34a33","orangered")

setwd("//allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr/snRNA 10xV3 Broad/dot plots")


# Function using modified_dot.R, does some filtering of the data and then writes both SVG of the plot
# and CSV containing the data and annotations necessary to reproduce it
source("modified_dot.R")
plot_all_dots <- function(data,
                          metadata,
                          genes,
                          grouping,
                          X,
                          colorset,
                          filter_level,
                          filter_by,
                          nodename,
                          setids,
                          ids,
                          dataname,
                          setnames) {
  
  X <- enquo(X)
  
  filter_level <- enquo(filter_level)
  
  anno <- metadata %>%
    filter(!!filter_level %in% filter_by)
  
  anno <- anno %>%
    mutate(X = as_factor(!!X)) 
  
  anno <- anno %>%
    arrange(!!X)
  
  anno$X <- as.factor(anno$X)
  levels(anno$X) <- setnames
  catorder <- setnames
  anno$X <- factor(anno$X, catorder)
  
  
  
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
                       font_size = 10,
                       max_width = 17,
                       rotate_counts = TRUE)
  print(plt)
  
  
  # Write SVG and CSV files to local directory
  nodename <- gsub("/","-",nodename)
  ggsave(filename = paste0("gene expression_",dataname,nodename,".svg"))
  alldat <- cbind.data.frame(anno,data[data$sample_name %in% anno$sample_name,])
  alldat <- alldat[!duplicated(as.list(alldat))]
  write.table(alldat, paste0("gene expression_",dataname, nodename, "_data.csv"), sep=",", 
              row.names = F)
  
}




# Class plots
for (i in 1:length(classnames)) {
  
  plot_all_dots(data = plotcounts,
                metadata = snRNA_10xv3_B_met,
                genes = genelist[[i]],
                grouping = "class",
                X = class_label,
                colorset = colorset,
                filter_level = class_label,
                filter_by = classnames,
                nodename = classnames[i],
                setids = classids,
                ids = classids$class_id,
                dataname = "snRNA10Xv3B_",
                setnames = classnames)
  
  
  
}

# Glut subclass plots

glutsubgenes <- genelist[,c(8:15)]
glutsubgenes$L6ITCar3 <- genelist$`L6 IT Car3`

glutsubgenes <- glutsubgenes[,c(1,2,3,7,9,3,6,8,5)]
glutsubgenesvec <- unique(as.vector(t(glutsubgenes)))

plotcounts <- cbind(sample_name = colnames(as.matrix(broadcounts)),
                    as.data.frame(t(as.matrix(broadcounts[glutsubgenesvec,]))))
plotcounts <- plotcounts[plotcounts$sample_name %in% snRNA_10xv3_B_met$sample_name,]



for (i in 1:length(glutsubnames)) {
  
  plot_all_dots(data = plotcounts,
                metadata = snRNA_10xv3_B_met,
                genes = glutsubgenes[[i]],
                grouping = "subclass",
                X = subclass_label,
                colorset = colorset,
                filter_level = class_label,
                filter_by = "Glutamatergic",
                nodename = glutsubnames[i],
                ids = glutsubids$subclass_id,
                dataname = "snRNA10Xv3B_",
                setnames = glutsubnames)
  
  
  
}


#############

# GABA subclass plots

gabasubgenes <- genelist[,gabasubnames]
gabasubgenesvec <- unique(as.vector(t(gabasubgenes)))

plotcounts <- cbind(sample_name = colnames(as.matrix(broadcounts)),
                    as.data.frame(t(as.matrix(broadcounts[gabasubgenesvec,]))))
plotcounts <- plotcounts[plotcounts$sample_name %in% snRNA_10xv3_B_met$sample_name,]


for (i in 1:length(gabasubnames)) {
  
  plot_all_dots(data = plotcounts,
                metadata = snRNA_10xv3_B_met,
                genes = gabasubgenes[[i]],
                grouping = "subclass",
                X = subclass_label,
                colorset = colorset,
                filter_level = class_label,
                filter_by = "GABAergic",
                nodename = gabasubnames[i],
                setids = gabasubids,
                ids = gabasubids$subclass_id,
                dataname = "snRNA10Xv3B_",
                setnames = gabasubnames)
  
  
  
}


# Non-neuronal subclass plots

colnames(genelist)[which(names(genelist) == "Astro")] <- "Astrocyte"
colnames(genelist)[which(names(genelist) == "Oligo")] <- "Oligodendrocyte"

nonneuronalgenes <- genelist[,nonneuronalsubnames]

nonneuronalgenesvec <- unique(as.vector(t(nonneuronalgenes)))

plotcounts <- cbind(sample_name = colnames(as.matrix(broadcounts)),
                    as.data.frame(t(as.matrix(broadcounts[nonneuronalgenesvec,]))))
plotcounts <- plotcounts[plotcounts$sample_name %in% snRNA_10xv3_B_met$sample_name,]




for (i in 1:length(nonneuronalsubnames)) {
  
  plot_all_dots(data = plotcounts,
                metadata = snRNA_10xv3_B_met,
                genes = nonneuronalgenes[[i]],
                grouping = "subclass",
                X = subclass_label,
                colorset = colorset,
                filter_level = class_label,
                filter_by = "Non-Neuronal",
                nodename = nonneuronalsubnames[i],
                setids = nonneuronalsubids,
                ids = nonneuronalsubids$subclass_id,
                dataname = "snRNA10Xv3B_",
                setnames = nonneuronalsubnames)
  
  
  
}
#######

# Non-neural subclass plots

colnames(genelist)[which(names(genelist) == "Endo")] <- "Endothelial"
colnames(genelist)[which(names(genelist) == "Micro-PVM")] <- "Microglia-PVM"
colnames(genelist)[which(names(genelist) == "Peri")] <- "Pericyte"

nonneuralgenes <- genelist[,nonneuralsubnames]

nonneuralgenesvec <- unique(as.vector(t(nonneuralgenes)))

plotcounts <- cbind(sample_name = colnames(as.matrix(broadcounts)),
                    as.data.frame(t(as.matrix(broadcounts[nonneuralgenesvec,]))))
plotcounts <- plotcounts[plotcounts$sample_name %in% snRNA_10xv3_B_met$sample_name,]




for (i in 1:length(nonneuralsubnames)) {
  
  plot_all_dots(data = plotcounts,
                metadata = snRNA_10xv3_B_met,
                genes = nonneuralgenes[[i]],
                grouping = "subclass",
                X = subclass_label,
                colorset = colorset,
                filter_level = class_label,
                filter_by = "Non-Neural",
                nodename = nonneuralsubnames[i],
                setids = nonneuralsubids,
                ids = nonneuralsubids$subclass_id,
                dataname = "snRNA10Xv3B_",
                setnames = nonneuralsubnames)
  
  
  
}


setwd("//allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr/scRNA 10xV2 AIBS/dot plots/test")

# Figuring out how to plot types

typegenes <- genelist
gabaterms <- c("Sncg","Vip","Pvalb","Meis2","Pvalb","Lamp5")


gabatypegenes <- genelist %>% 
  select_at(vars(matches(paste(gabaterms, collapse="|")))) 
gabatypegenes <- gabatypegenes[,7:ncol(gabatypegenes)]



# Doing this in not the best way
# Lamp5
#####
lamp5genes <- gabatypegenes %>% 
  select_at(vars(matches(paste("Lamp5", collapse="|")))) 

lamp5genes <- cbind(lamp5genes, genelist$`Lamp5 Egln3_1`)


lamp5genesvec <- unique(as.vector(t(lamp5genes)))

plotcounts <- cbind(sample_name = colnames(as.matrix(broadcounts)),
                    as.data.frame(t(as.matrix(broadcounts[lamp5genesvec,]))))
plotcounts <- plotcounts[plotcounts$sample_name %in% snRNA_10xv3_B_met$sample_name,]

lamp5ids <- typeids %>% 
  filter(str_detect(cluster_label, "Lamp5"))
lamp5names <- lamp5ids$cluster_label
colnames(lamp5genes)[ncol(lamp5genes)] <- "Lamp5 Egln3_1"
lamp5genes <- lamp5genes[lamp5names]


for (i in 1:length(lamp5names)) {
  
  plot_all_dots(data = plotcounts,
                metadata = snRNA_10xv3_B_met,
                genes = lamp5genes[[i]],
                grouping = "cluster",
                X = cluster_label,
                colorset = colorset,
                filter_level = subclass_label,
                filter_by = "Lamp5",
                nodename = lamp5names[i],
                setids = lamp5ids,
                ids = lamp5ids$cluster_id,
                dataname = "snRNA10Xv3B_",
                setnames = lamp5names)
  
  
  
}

# Sncg
#######
# Doing this in not the best way
sncggenes <- gabatypegenes %>% 
  select_at(vars(matches(paste("Sncg", collapse="|")))) 




sncggenesvec <- unique(as.vector(t(sncggenes)))

plotcounts <- cbind(sample_name = colnames(as.matrix(broadcounts)),
                    as.data.frame(t(as.matrix(broadcounts[sncggenesvec,]))))
plotcounts <- plotcounts[plotcounts$sample_name %in% snRNA_10xv3_B_met$sample_name,]

sncgids <- typeids %>% 
  filter(str_detect(cluster_label, "Sncg"))
sncgnames <- sncgids$cluster_label
sncggenes <- sncggenes[sncgnames]


for (i in 1:length(sncgnames)) {
  
  plot_all_dots(data = plotcounts,
                metadata = snRNA_10xv3_B_met,
                genes = sncggenes[[i]],
                grouping = "cluster",
                X = cluster_label,
                colorset = colorset,
                filter_level = subclass_label,
                filter_by = "Sncg",
                nodename = sncgnames[i],
                setids = sncgids,
                ids = sncgids$cluster_id,
                dataname = "snRNA10Xv3B_",
                setnames = sncgnames)
  
  
  
}
# Vip
########
# Vip

vipgenes <- gabatypegenes %>% 
  select_at(vars(matches(paste("Vip", collapse="|")))) 

vipids <- clupdate %>% 
  filter(subclass_label == "Vip") %>%
  select(cluster_id, cluster_label, cluster_color)

vipnames <- vipids$cluster_label
colnames(vipgenes)[colnames(vipgenes) == "Vip Htr1f_1"] <- "Vip Htr1f"
vipgenes <- vipgenes[vipnames]


vipgenesvec <- unique(as.vector(t(vipgenes)))

plotcounts <- cbind(sample_name = colnames(as.matrix(broadcounts)),
                    as.data.frame(t(as.matrix(broadcounts[vipgenesvec,]))))
plotcounts <- plotcounts[plotcounts$sample_name %in% snRNA_10xv3_B_met$sample_name,]




for (i in 1:length(vipnames)) {
  
  plot_all_dots(data = plotcounts,
                metadata = snRNA_10xv3_B_met,
                genes = vipgenes[[i]],
                grouping = "cluster",
                X = cluster_label,
                colorset = colorset,
                filter_level = subclass_label,
                filter_by = "Vip",
                nodename = vipnames[i],
                setids = vipids,
                ids = vipids$cluster_id,
                dataname = "snRNA10Xv3B_",
                setnames = vipnames)
  
  
  
}


######
# sst

sstgenes <- genelist %>% 
  select_at(vars(matches(paste("Sst", collapse="|")))) 

sstids <- clupdate %>% 
  filter(subclass_label == "Sst") %>%
  select(cluster_label,cluster_id, cluster_color)

sstnames <- sstids$cluster_label
colnames(sstgenes)[colnames(sstgenes) == "Sst_Pappa"] <- "Sst Pappa"
colnames(sstgenes)[colnames(sstgenes) == "Sst Etv1_1"] <- "Sst Etv1"
sstgenes <- sstgenes[sstnames]


sstgenesvec <- unique(as.vector(t(sstgenes)))

plotcounts <- cbind(sample_name = colnames(as.matrix(broadcounts)),
                    as.data.frame(t(as.matrix(broadcounts[sstgenesvec,]))))
plotcounts <- plotcounts[plotcounts$sample_name %in% snRNA_10xv3_B_met$sample_name,]

#######


for (i in 1:length(sstnames)) {
  
  plot_all_dots(data = plotcounts,
                metadata = snRNA_10xv3_B_met,
                genes = sstgenes[[i]],
                grouping = "cluster",
                X = cluster_label,
                colorset = colorset,
                filter_level = subclass_label,
                filter_by = "Sst",
                nodename = sstnames[i],
                setids = sstids,
                ids = sstids$cluster_id,
                dataname = "snRNA10Xv3B_",
                setnames = sstnames)
  
  
  
}


# Meis2

meis2genes <- gabatypegenes %>% 
  select_at(vars(matches(paste("Meis2", collapse="|")))) 

meis2ids <- clupdate %>% 
  filter(subclass_label == "Meis2") %>%
  select(cluster_label,cluster_id, cluster_color)

meis2names <- meis2ids$cluster_label
colnames(meis2genes)[1] <- "Meis2"


meis2genesvec <- unique(as.vector(t(meis2genes)))

plotcounts <- cbind(sample_name = colnames(as.matrix(broadcounts)),
                    as.data.frame(t(as.matrix(broadcounts[meis2genesvec,]))))
plotcounts <- plotcounts[plotcounts$sample_name %in% snRNA_10xv3_B_met$sample_name,]


for (i in 1:length(meis2names)) {
  
  plot_all_dots(data = plotcounts,
                metadata = snRNA_10xv3_B_met,
                genes = meis2genes[[i]],
                grouping = "cluster",
                X = cluster_label,
                colorset = colorset,
                filter_level = subclass_label,
                filter_by = "Meis2",
                nodename = meis2names[i],
                setids = meis2ids,
                ids = meis2ids$cluster_id,
                dataname = "snRNA10Xv3B_",
                setnames = meis2names)
  
  
  
}
#####


# Pvalb

pvalbgenes <- gabatypegenes %>% 
  select_at(vars(matches(paste("Pvalb", collapse="|")))) 

pvalbids <- clupdate %>% 
  filter(subclass_label == "Pvalb") %>%
  select(cluster_label,cluster_id, cluster_color)

pvalbnames <- pvalbids$cluster_label
pvalbgenes <- pvalbgenes[pvalbnames]
pvalbgenesvec <- unique(as.vector(t(pvalbgenes)))

plotcounts <- cbind(sample_name = colnames(as.matrix(broadcounts)),
                    as.data.frame(t(as.matrix(broadcounts[pvalbgenesvec,]))))
plotcounts <- plotcounts[plotcounts$sample_name %in% snRNA_10xv3_B_met$sample_name,]

pvalbgenes2 <- pvalbgenes
pvalbgenes2[pvalbgenes2 == "Nkx2-1"] <- "Nkx21"
colnames(plotcounts)[colnames(plotcounts) == "Nkx2-1"] <- "Nkx21"


for (i in 1:length(pvalbnames)) {
  
  plot_all_dots(data = plotcounts,
                metadata = snRNA_10xv3_B_met,
                genes = pvalbgenes2[[i]],
                grouping = "cluster",
                X = cluster_label,
                colorset = colorset,
                filter_level = subclass_label,
                filter_by = "Pvalb",
                nodename = pvalbnames[i],
                setids = pvalbids,
                ids = pvalbids$cluster_id,
                dataname = "snRNA10Xv3B_",
                setnames = pvalbnames)
  
  
  
}



######

# L2/3 IT

L23ITgenes <- genelist %>% 
  select_at(vars(matches(paste("L2-3 IT", collapse="|")))) 
L23ITgenes <- L23ITgenes[2:4]

L23ITids <- clupdate %>% 
  filter(subclass_label == "L2/3 IT") %>%
  select(cluster_label,cluster_id, cluster_color)

L23ITnames <- L23ITids$cluster_label



L23ITgenesvec <- unique(as.vector(t(L23ITgenes)))

plotcounts <- cbind(sample_name = colnames(as.matrix(broadcounts)),
                    as.data.frame(t(as.matrix(broadcounts[L23ITgenesvec,]))))
plotcounts <- plotcounts[plotcounts$sample_name %in% snRNA_10xv3_B_met$sample_name,]



for (i in 1:length(L23ITnames)) {
  
  plot_all_dots(data = plotcounts,
                metadata = snRNA_10xv3_B_met,
                genes = L23ITgenes[[i]],
                grouping = "cluster",
                X = cluster_label,
                colorset = colorset,
                filter_level = subclass_label,
                filter_by = "L2/3 IT",
                nodename = L23ITnames[i],
                setids = L23ITids,
                ids = L23ITids$cluster_id,
                dataname = "snRNA10Xv3B_",
                setnames = L23ITnames)
  
  
  
}


#########

# L4/5 IT

L45ITgenes <- genelist %>% 
  select_at(vars(matches(paste("L4-5 IT", collapse="|")))) 
L45ITgenes <- L45ITgenes[2:3]

L45ITids <- clupdate %>% 
  filter(subclass_label == "L4/5 IT") %>%
  select(cluster_label,cluster_id, cluster_color)

L45ITnames <- L45ITids$cluster_label
colnames(L45ITgenes)[colnames(L45ITgenes) == "L4-5 IT_1"] <- "L4/5 IT_1"
colnames(L45ITgenes)[colnames(L45ITgenes) == "L4-5 IT_2"] <- "L4/5 IT_2"
L45ITgenes <- L45ITgenes[L45ITnames]


L45ITgenesvec <- unique(as.vector(t(L45ITgenes)))

plotcounts <- cbind(sample_name = colnames(as.matrix(broadcounts)),
                    as.data.frame(t(as.matrix(broadcounts[L45ITgenesvec,]))))
plotcounts <- plotcounts[plotcounts$sample_name %in% snRNA_10xv3_B_met$sample_name,]


for (i in 1:length(L45ITnames)) {
  
  plot_all_dots(data = plotcounts,
                metadata = snRNA_10xv3_B_met,
                genes = L45ITgenes[[i]],
                grouping = "cluster",
                X = cluster_label,
                colorset = colorset,
                filter_level = subclass_label,
                filter_by = "L4/5 IT",
                nodename = L45ITnames[i],
                setids = L45ITids,
                ids = L45ITids$cluster_id,
                dataname = "snRNA10Xv3B_",
                setnames = L45ITnames)
  
  
  
}
######
# L5 IT

L5ITgenes <- genelist %>% 
  select_at(vars(matches(paste("L5 IT", collapse="|")))) 
L5ITgenes <- L5ITgenes[2:5]

L5ITids <- clupdate %>% 
  filter(subclass_label == "L5 IT") %>%
  select(cluster_label,cluster_id, cluster_color)

L5ITnames <- L5ITids$cluster_label
L5ITgenes <- L5ITgenes[L5ITnames]


L5ITgenesvec <- unique(as.vector(t(L5ITgenes)))

plotcounts <- cbind(sample_name = colnames(as.matrix(broadcounts)),
                    as.data.frame(t(as.matrix(broadcounts[L5ITgenesvec,]))))
plotcounts <- plotcounts[plotcounts$sample_name %in% snRNA_10xv3_B_met$sample_name,]


for (i in 1:length(L5ITnames)) {
  
  plot_all_dots(data = plotcounts,
                metadata = snRNA_10xv3_B_met,
                genes = L5ITgenes[[i]],
                grouping = "cluster",
                X = cluster_label,
                colorset = colorset,
                filter_level = subclass_label,
                filter_by = "L5 IT",
                nodename = L5ITnames[i],
                setids = L5ITids,
                ids = L5ITids$cluster_id,
                dataname = "snRNA10Xv3B_",
                setnames = L5ITnames)
  
  
  
}

# L6 IT
############
# L6 IT

L6ITgenes <- genelist %>% 
  select_at(vars(matches(paste("L6 IT", collapse="|")))) 
L6ITgenes <- L6ITgenes[3:4]

L6ITids <- clupdate %>% 
  filter(subclass_label == "L6 IT") %>%
  select(cluster_label,cluster_id, cluster_color)

L6ITnames <- L6ITids$cluster_label
L6ITgenes <- L6ITgenes[L6ITnames]


L6ITgenesvec <- unique(as.vector(t(L6ITgenes)))

plotcounts <- cbind(sample_name = colnames(as.matrix(broadcounts)),
                    as.data.frame(t(as.matrix(broadcounts[L6ITgenesvec,]))))
plotcounts <- plotcounts[plotcounts$sample_name %in% snRNA_10xv3_B_met$sample_name,]


for (i in 1:length(L6ITnames)) {
  
  plot_all_dots(data = plotcounts,
                metadata = snRNA_10xv3_B_met,
                genes = L6ITgenes[[i]],
                grouping = "cluster",
                X = cluster_label,
                colorset = colorset,
                filter_level = subclass_label,
                filter_by = "L6 IT",
                nodename = L6ITnames[i],
                setids = L6ITids,
                ids = L6ITids$cluster_id,
                dataname = "snRNA10Xv3B_",
                setnames = L6ITnames)
  
  
  
}

# L5 ET
#######
# L5 ET

L5ETgenes <- genelist %>% 
  select_at(vars(matches(paste("L5 ET", collapse="|")))) 
L5ETgenes <- L5ETgenes[2:5]

L5ETids <- clupdate %>% 
  filter(subclass_label == "L5 ET") %>%
  select(cluster_label,cluster_id, cluster_color)

L5ETnames <- L5ETids$cluster_label
L5ETgenes <- L5ETgenes[L5ETnames]


L5ETgenesvec <- unique(as.vector(t(L5ETgenes)))

plotcounts <- cbind(sample_name = colnames(as.matrix(broadcounts)),
                    as.data.frame(t(as.matrix(broadcounts[L5ETgenesvec,]))))
plotcounts <- plotcounts[plotcounts$sample_name %in% snRNA_10xv3_B_met$sample_name,]


for (i in 1:length(L5ETnames)) {
  
  plot_all_dots(data = plotcounts,
                metadata = snRNA_10xv3_B_met,
                genes = L5ETgenes[[i]],
                grouping = "cluster",
                X = cluster_label,
                colorset = colorset,
                filter_level = subclass_label,
                filter_by = "L5 ET",
                nodename = L5ETnames[i],
                setids = L5ETids,
                ids = L5ETids$cluster_id,
                dataname = "snRNA10Xv3B_",
                setnames = L5ETnames)
  
  
  
}

# L5 ET

# L5/6 NP
#######
# L5/6 NP

L56NPgenes <- genelist %>% 
  select_at(vars(matches(paste("L5-6 NP", collapse="|")))) 
L56NPgenes <- L56NPgenes[2:5]

colnames(L56NPgenes) <- gsub("-","/",colnames(L56NPgenes))


L56NPids <- clupdate %>% 
  filter(subclass_label == "L5/6 NP") %>%
  select(cluster_label,cluster_id, cluster_color)

L56NPnames <- L56NPids$cluster_label
L56NPgenes <- L56NPgenes[L56NPnames]


L56NPgenesvec <- unique(as.vector(t(L56NPgenes)))

plotcounts <- cbind(sample_name = colnames(as.matrix(broadcounts)),
                    as.data.frame(t(as.matrix(broadcounts[L56NPgenesvec,]))))
plotcounts <- plotcounts[plotcounts$sample_name %in% snRNA_10xv3_B_met$sample_name,]


for (i in 1:length(L56NPnames)) {
  
  plot_all_dots(data = plotcounts,
                metadata = snRNA_10xv3_B_met,
                genes = L56NPgenes[[i]],
                grouping = "cluster",
                X = cluster_label,
                colorset = colorset,
                filter_level = subclass_label,
                filter_by = "L5/6 NP",
                nodename = L56NPnames[i],
                setids = L56NPids,
                ids = L56NPids$cluster_id,
                dataname = "snRNA10Xv3B_",
                setnames = L56NPnames)
  
  
  
}

# L6 CT
#######
# L6 CT

L6CTgenes <- genelist %>% 
  select_at(vars(matches(paste("CT", collapse="|")))) 
L6CTgenes <- L6CTgenes[2:8]

colnames(L6CTgenes) <- gsub("-","/",colnames(L6CTgenes))


L6CTids <- clupdate %>% 
  filter(subclass_label == "L6 CT") %>%
  select(cluster_label,cluster_id, cluster_color)

L6CTnames <- L6CTids$cluster_label
L6CTgenes <- L6CTgenes[L6CTnames]


L6CTgenesvec <- unique(as.vector(t(L6CTgenes)))

plotcounts <- cbind(sample_name = colnames(as.matrix(broadcounts)),
                    as.data.frame(t(as.matrix(broadcounts[L6CTgenesvec,]))))
plotcounts <- plotcounts[plotcounts$sample_name %in% snRNA_10xv3_B_met$sample_name,]


for (i in 1:length(L6CTnames)) {
  
  plot_all_dots(data = plotcounts,
                metadata = snRNA_10xv3_B_met,
                genes = L6CTgenes[[i]],
                grouping = "cluster",
                X = cluster_label,
                colorset = colorset,
                filter_level = subclass_label,
                filter_by = "L6 CT",
                nodename = L6CTnames[i],
                setids = L6CTids,
                ids = L6CTids$cluster_id,
                dataname = "snRNA10Xv3B_",
                setnames = L6CTnames)
  
  
  
}

# L6 CT

# L6b
#######
# L6b

L6bgenes <- genelist %>% 
  select_at(vars(matches(paste("L6b", collapse="|")))) 
L6bgenes <- L6bgenes[2:6]

colnames(L6bgenes) <- gsub("-","/",colnames(L6bgenes))


L6bids <- clupdate %>% 
  filter(subclass_label == "L6b") %>%
  select(cluster_label,cluster_id, cluster_color)

L6bnames <- L6bids$cluster_label
L6bgenes <- L6bgenes[L6bnames]


L6bgenesvec <- unique(as.vector(t(L6bgenes)))

plotcounts <- cbind(sample_name = colnames(as.matrix(broadcounts)),
                    as.data.frame(t(as.matrix(broadcounts[L6bgenesvec,]))))
plotcounts <- plotcounts[plotcounts$sample_name %in% snRNA_10xv3_B_met$sample_name,]


for (i in 1:length(L6bnames)) {
  
  plot_all_dots(data = plotcounts,
                metadata = snRNA_10xv3_B_met,
                genes = L6bgenes[[i]],
                grouping = "cluster",
                X = cluster_label,
                colorset = colorset,
                filter_level = subclass_label,
                filter_by = "L6b",
                nodename = L6bnames[i],
                setids = L6bids,
                ids = L6bids$cluster_id,
                dataname = "snRNA10Xv3B_",
                setnames = L6bnames)
  
  
  
}


# Astro
#######

Astrogenes <- genelist %>% 
  select_at(vars(matches(paste("Astro", collapse="|")))) 
Astrogenes <- Astrogenes[2:4]

Astroids <- clupdate %>% 
  filter(subclass_label == "Astrocyte") %>%
  select(cluster_label,cluster_id, cluster_color)

Astronames <- Astroids$cluster_label
Astrogenes <- Astrogenes[Astronames]


Astrogenesvec <- unique(as.vector(t(Astrogenes)))

plotcounts <- cbind(sample_name = colnames(as.matrix(broadcounts)),
                    as.data.frame(t(as.matrix(broadcounts[Astrogenesvec,]))))
plotcounts <- plotcounts[plotcounts$sample_name %in% snRNA_10xv3_B_met$sample_name,]


for (i in 1:length(Astronames)) {
  
  plot_all_dots(data = plotcounts,
                metadata = snRNA_10xv3_B_met,
                genes = Astrogenes[[i]],
                grouping = "cluster",
                X = cluster_label,
                colorset = colorset,
                filter_level = subclass_label,
                filter_by = "Astrocyte",
                nodename = Astronames[i],
                setids = Astroids,
                ids = Astroids$cluster_id,
                dataname = "snRNA10Xv3B_",
                setnames = Astronames)
  
  
  
}


# Oligo
############
# Oligo

oligogenes <- genelist %>% 
  select_at(vars(matches(paste("Oligo", collapse="|")))) 
oligogenes <- oligogenes[2:9]




oligoids <- clupdate %>% 
  filter(subclass_label == "Oligodendrocyte") %>%
  select(cluster_label,cluster_id, cluster_color)

oligonames <- oligoids$cluster_label
oligogenes <- oligogenes[oligonames]


oligogenesvec <- unique(as.vector(t(oligogenes)))

plotcounts <- cbind(sample_name = colnames(as.matrix(broadcounts)),
                    as.data.frame(t(as.matrix(broadcounts[oligogenesvec,]))))
plotcounts <- plotcounts[plotcounts$sample_name %in% snRNA_10xv3_B_met$sample_name,]


for (i in 1:length(oligonames)) {
  
  plot_all_dots(data = plotcounts,
                metadata = snRNA_10xv3_B_met,
                genes = oligogenes[[i]],
                grouping = "cluster",
                X = cluster_label,
                colorset = colorset,
                filter_level = subclass_label,
                filter_by = "Oligodendrocyte",
                nodename = oligonames[i],
                setids = oligoids,
                ids = oligoids$cluster_id,
                dataname = "snRNA10Xv3B_",
                setnames = oligonames)
  
  
  
}


# VLMC
############

vlmcgenes <- genelist %>% 
  select_at(vars(matches(paste("VLMC", collapse="|")))) 
vlmcgenes <- vlmcgenes[2:8]




vlmcids <- clupdate %>% 
  filter(subclass_label == "VLMC") %>%
  select(cluster_label,cluster_id, cluster_color)

vlmcnames <- vlmcids$cluster_label
vlmcgenes <- vlmcgenes[vlmcnames]


vlmcgenesvec <- unique(as.vector(t(vlmcgenes)))

plotcounts <- cbind(sample_name = colnames(as.matrix(broadcounts)),
                    as.data.frame(t(as.matrix(broadcounts[vlmcgenesvec,]))))
plotcounts <- plotcounts[plotcounts$sample_name %in% snRNA_10xv3_B_met$sample_name,]


for (i in 1:length(vlmcnames)) {
  
  plot_all_dots(data = plotcounts,
                metadata = snRNA_10xv3_B_met,
                genes = vlmcgenes[[i]],
                grouping = "cluster",
                X = cluster_label,
                colorset = colorset,
                filter_level = subclass_label,
                filter_by = "VLMC",
                nodename = vlmcnames[i],
                setids = vlmcids,
                ids = vlmcids$cluster_id,
                dataname = "snRNA10Xv3B_",
                setnames = vlmcnames)
  
  
  
}


# Micro-PVM
#######

pvmgenes <- genelist %>% 
  select_at(vars(matches(paste("PVM", collapse="|")))) 

microgenes <- genelist %>%
  select_at(vars(matches(paste("Micro", collapse = "|"))))

micropvmgenes <- cbind(microgenes, pvmgenes)
micropvmgenes <- micropvmgenes[, c(2,4,5,6)]

micropvmids <- clupdate %>% 
  filter(subclass_label == "Microglia-PVM") %>%
  select(cluster_label,cluster_id, cluster_color)

micropvmnames <- micropvmids$cluster_label
micropvmgenes <- micropvmgenes[micropvmnames]
micropvmgenesvec <- unique(as.vector(t(micropvmgenes)))

plotcounts <- cbind(sample_name = colnames(as.matrix(broadcounts)),
                    as.data.frame(t(as.matrix(broadcounts[micropvmgenesvec,]))))
plotcounts <- plotcounts[plotcounts$sample_name %in% snRNA_10xv3_B_met$sample_name,]


for (i in 1:length(micropvmnames)) {
  
  plot_all_dots(data = plotcounts,
                metadata = snRNA_10xv3_B_met,
                genes = micropvmgenes[[i]],
                grouping = "cluster",
                X = cluster_label,
                colorset = colorset,
                filter_level = subclass_label,
                filter_by = "Microglia-PVM",
                nodename = micropvmnames[i],
                setids = micropvmids,
                ids = micropvmids$cluster_id,
                dataname = "snRNA10Xv3B_",
                setnames = micropvmnames)
  
  
  
}