# This script takes the gene-UMAP files created by Cindy and creates the "Panel_Image" CSV for the product team,
# essentially outlining on which cell type cards each of the UMAP files should be included. 

setwd("/allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr")
MOpgenes <- read_csv("MOP_genelist.csv")
MOpgenes <- MOpgenes[1:10,]


setwd("/allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr/umap_Yao.markers_by_platform")
indir <- "/allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr/umap_Yao.markers_by_platform"

filenames <- list.files(indir, pattern = ".png")
genes <- unique(sub("_.*","",filenames))

colnames(MOpgenes) <- gsub("-","/",colnames(MOpgenes))
colnames(MOpgenes) <- gsub("Meis2...17","Meis2",colnames(MOpgenes))
colnames(MOpgenes) <- gsub("Meis2...67","Meis2_type",colnames(MOpgenes))
colnames(MOpgenes) <- gsub("Non/Neuronal","Non-Neuronal",colnames(MOpgenes))
colnames(MOpgenes) <- gsub("Non/Neural","Non-Neural",colnames(MOpgenes))

taxons <- list()
temp <- c()
for (i in 1:length(genes)) {
  for (j in 1:ncol(MOpgenes)) {
    if (any(genes[i] %in% MOpgenes[[j]]) == TRUE) {

      temp <- c(temp,colnames(MOpgenes)[j])
      
    }
    taxons[[i]] <- temp
  }
  temp <- c()
}

taxons
test <- taxons
names(test) <- genes
test[sapply(test, is.null)] <- NULL

gene_list <- rep(names(test), times = lapply(test,length))
gene_list <- rep(gene_list,each = 7)
taxon_list <- unlist(unname(test))
taxon_list <- rep(taxon_list,each = 7)


tempdf <- cbind.data.frame(Taxon = taxon_list, Gene_Name = gene_list)
filereps <- as.numeric(lapply(unname(test),length))
files2rep <- split(filenames,ceiling(seq_along(filenames)/7))
files2rep <- files2rep[c(1:4,6:13,15:18,20:26,28,29,31:length(files2rep))]

empty <- c()

for (i in 1:length(filereps)) {
  temp <- rep(files2rep[i],filereps[i])
  empty <- c(empty,temp)
}

completefilenames <- unname(unlist(empty))

tempdf$UMAP_Image_Filename <- completefilenames
tempdf$Taxon[tempdf$Taxon == "Micro/PVM"] <- "Microglia-PVM"
tempdf$Taxon[tempdf$Taxon == "Peri"] <- "Pericyte"
tempdf$Taxon[tempdf$Taxon == "Astro"] <- "Astrocyte"
tempdf$Taxon[tempdf$Taxon == "Oligo"] <- "Oligodendrocyte"
tempdf$Taxon[tempdf$Taxon == "Endo"] <- "Endothelial"
tempdf$Taxon[tempdf$Taxon == "OPC Pdgfra"] <- "OPC"

cellsetnom2 <- cellsetnom[!is.na(cellsetnom$`cell set preferred alias`),]

emptyvec <- c()
for (i in 1:nrow(tempdf)) {
  if (tempdf$Taxon[i] %in% cellsetnom$`cell set preferred alias` == TRUE){
    temp <- cellsetnom2$`cell set accession`[tempdf$Taxon[i] == cellsetnom2$`cell set preferred alias`]
    emptyvec <- c(emptyvec,temp)
}
}

tempdf$Accession <- emptyvec

allds <- c("scRNA 10xV2 AIBS","scRNA 10xV3 AIBS", "snRNA 10xV2 AIBS","snRNA 10xV3 AIBS",
           "snRNA 10xV3 Broad","scRNA SMART","snRNA SMART")
allds <- rep(allds,nrow(tempdf)/7)
tempdf$DatasetID <- allds
write.table(tempdf,"YaoMarkersUMAP.csv",sep=",",row.names=F)



setwd("/allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr/UMAPGenesCSV")
dataset <- ldply(list.files(), read.csv, header=TRUE)
dataset <- unique(dataset)
write.table(dataset, "UMAP-Gene Mapping.csv",sep=",",row.names = F)

###########

# NSF Markers

# Starting to try to get the UMAP-gene mapping together

setwd("/allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr")
MOpgenes <- read_csv("MOP_genelist.csv")
MOpgenes <- MOpgenes[1:10,]


setwd("/allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr/umap_NSF.markers_by_platform")
indir <- "/allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr/umap_NSF.markers_by_platform"

filenames <- list.files(indir, pattern = ".png")
genes <- unique(sub("_.*","",filenames))

colnames(MOpgenes) <- gsub("-","/",colnames(MOpgenes))
colnames(MOpgenes) <- gsub("Meis2...17","Meis2",colnames(MOpgenes))
colnames(MOpgenes) <- gsub("Meis2...67","Meis2_type",colnames(MOpgenes))
colnames(MOpgenes) <- gsub("Non/Neuronal","Non-Neuronal",colnames(MOpgenes))
colnames(MOpgenes) <- gsub("Non/Neural","Non-Neural",colnames(MOpgenes))

taxons <- list()
temp <- c()
for (i in 1:length(genes)) {
  for (j in 1:ncol(MOpgenes)) {
    if (any(genes[i] %in% MOpgenes[[j]]) == TRUE) {
      
      temp <- c(temp,colnames(MOpgenes)[j])
      
    }
    taxons[[i]] <- temp
  }
  temp <- c()
}

taxons
test <- taxons
names(test) <- genes
test[sapply(test, is.null)] <- NULL

gene_list <- rep(names(test), times = lapply(test,length))
gene_list <- rep(gene_list,each = 7)
taxon_list <- unlist(unname(test))
taxon_list <- rep(taxon_list,each = 7)


tempdf <- cbind.data.frame(Taxon = taxon_list, Gene_Name = gene_list)
filereps <- as.numeric(lapply(unname(test),length))
files2rep <- split(filenames,ceiling(seq_along(filenames)/7))
files2rep <- files2rep[c(1,4:226,228:length(files2rep))]

empty <- c()

for (i in 1:length(filereps)) {
  temp <- rep(files2rep[i],filereps[i])
  empty <- c(empty,temp)
}

completefilenames <- unname(unlist(empty))

tempdf$UMAP_Image_Filename <- completefilenames
tempdf$Taxon[tempdf$Taxon == "Micro/PVM"] <- "Microglia-PVM"
tempdf$Taxon[tempdf$Taxon == "Peri"] <- "Pericyte"
tempdf$Taxon[tempdf$Taxon == "Astro"] <- "Astrocyte"
tempdf$Taxon[tempdf$Taxon == "Oligo"] <- "Oligodendrocyte"
tempdf$Taxon[tempdf$Taxon == "Endo"] <- "Endothelial"
tempdf$Taxon[tempdf$Taxon == "OPC Pdgfra"] <- "OPC"

cellsetnom2 <- cellsetnom[!is.na(cellsetnom$`cell set preferred alias`),]

emptyvec <- c()
for (i in 1:nrow(tempdf)) {
  if (tempdf$Taxon[i] %in% cellsetnom$`cell set preferred alias` == TRUE){
    temp <- cellsetnom2$`cell set accession`[tempdf$Taxon[i] == cellsetnom2$`cell set preferred alias`]
    emptyvec <- c(emptyvec,temp)
  }
}

tempdf$Accession <- emptyvec

allds <- c("scRNA 10xV2 AIBS","scRNA 10xV3 AIBS", "snRNA 10xV2 AIBS","snRNA 10xV3 AIBS",
           "snRNA 10xV3 Broad","scRNA SMART","snRNA SMART")
allds <- rep(allds,nrow(tempdf)/7)
tempdf$DatasetID <- allds
write.table(tempdf,"NSFMarkersUMAP.csv",sep=",",row.names=F)
