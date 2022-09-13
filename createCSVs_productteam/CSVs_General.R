# Auto-generate CSVs with critical metadata needed for the product team to populate cell type cards

# This script defines the steps taken to create the CSVs requested by the product team for each panel, here using 
# the dotplots I generated for the Broad dataset as an example. The columns of the final output CSV include: Panel_Image 
# (a text description of the image), Panel_Title (in this case, "Gene Expression"), Taxon_Accession_ID, Dataset_ID (in this case
# the Broad dataset), Image_File_Name, and Order (the order each image should appear, from left to right or top to bottom if it's
# for one of the dropdown menus).
#
# This approach was used for all datasets for all panels:
# - 1. Load Accession IDs and metadata
# - 2. Clean names as necessary to match files
# - 3. Read in files and create vectors of the easy stuff (Dataset_ID for example)
# - 4. Match Accession IDs and cbind everything to a dataframe
# - 5. Write dataframe to CSV
# - 6. Rinse & repat

# CSVs can then be combined using a simple script (load and cbind or use something like vroom() )

setwd("/allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr/")
#clupdate <- read_csv("cl_updated.df.csv")
nomenclature <- read_csv("Accesion_IDs_MOp.csv")
nomenclature$`cell set preferred alias`[119] <- "Non-Neural"
nomenclature$`cell set preferred alias`[136] <- "OPC"
nomenclature$`cell set preferred alias`[130] <- "Astrocyte"
nomenclature$`cell set preferred alias`[131] <- "Oligodendrocyte"
nomenclature$`cell set preferred alias`[133] <- "Microglia-PVM"
nomenclature$`cell set preferred alias`[103] <- "Endothelial"
nomenclature$`cell set preferred alias`[112] <- "Pericyte"
mousenomtypes <- nomenclature[1:116,]
mousenomother <- nomenclature[117:140,]

dir <- "/allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr/snRNA 10xV3 Broad/dot plots/test"
#setwd(dir)
filenames <- list.files(path = dir, pattern = ".svg")
paneltitles <- rep("Gene Expression", length(filenames))
dataset <- rep("snRNA 10xV3 Broad", length(filenames))
order <- rep(1,length(filenames))


taxons <- sub(".*v3B_","",filenames)
taxons <- gsub(".svg","",taxons)
taxons <- gsub("-","/",taxons)


overlapnom <- mousenomtypes[mousenomtypes$`cell set preferred alias` %in% taxons,]
overlapnom <- overlapnom[order(overlapnom$`cell set preferred alias`),]


descriptions <- paste0("Dot plot of ", taxons, " marker gene expression, snRNA 10xV3 Broad")

BroadCSV <- cbind.data.frame(descriptions, paneltitles, overlapnom$`cell set accession`, dataset, filenames, order)
colnames(BroadCSV) <- c("Panel_Image","Panel_Title","Taxon_Accession_ID","Dataset_ID","Image_File_Name","Order")


## classes and subclasses

dir2 <- "/allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr/snRNA 10xV3 Broad/dot plots/"
#setwd(dir2)
filenames <- list.files(path = dir2,pattern = ".svg")
paneltitles <- rep("Gene Expression", length(filenames))
dataset <- rep("snRNA 10xV3 Broad", length(filenames))
order <- rep(1,length(filenames))

taxons <- sub(".*v3B_","",filenames)
taxons <- gsub(".svg","",taxons)
taxons[c(5,6,9)] <- gsub("-","/",taxons[c(5,6,9)])

overlapnom <- nomenclature[nomenclature$`cell set preferred alias` %in% taxons,]
overlapnom <- overlapnom[order(overlapnom$`cell set preferred alias`),]
overlapnom <- overlapnom[c(1:14,16:29),]
overlapnom <- overlapnom[order(taxons),]


descriptions <- paste0("Dot plot of ", taxons, " marker gene expression, snRNA 10xV3 Broad ")

BroadCSV2 <- cbind.data.frame(descriptions, paneltitles, overlapnom$`cell set accession`, dataset, filenames, order)
colnames(BroadCSV2) <- c("Panel_Image","Panel_Title","Taxon_Accession_ID","Dataset_ID","Image_File_Name","Order")
alldotBroad <- rbind.data.frame(BroadCSV,BroadCSV2)
write.table(alldotBroad, "GeneExpression_snRNA10xV3B_ImageInfo.csv",sep=",",row.names = F)
