##### Semi-automated short text descriptions for cell type cards. The basic steps are as follows:
#
# - 1. Read in and clean up data and metadata from the RNA-seq, MERFISH, and Patch-seq datasets
# - 2. Loop through all subclasses and for each subclass, calculate some basic metrics if they exist
# - 3. Paste metrics with predefined text based on Jeremy's template
# - 4. Write description to data frame
#
# The final output includes columns for taxon preferred alias, accession ID, and description. For the actual descriptions
# going into cell type cards, I did some manual editing of the descriptions output by this script, but this was minimal and 
# took about 5-10 minutes. 
#
# This same procedure can be repeated for classes and types with some minimal modification of the script. 

# Loading and cleaning data (not the most elegant but it works)
#######
# Setting up script to do automatic short text descriptions
setwd("/allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/merfish/")
merfishmeta <- read_csv("cell_metadata_final.csv")
merfishmeta$subclass[merfishmeta$subclass == "L5_PT"] <- "L5 ET"
merfishmeta$class_label[merfishmeta$subclass == "Endothelial" | 
                          merfishmeta$subclass == "Microglia-PVM"| merfishmeta$subclass == "SMC"| 
                          merfishmeta$subclass == "VLMC" | merfishmeta$subclass == "Pericyte"] <- "Non-Neural"
merfishmeta$class_label[merfishmeta$subclass == "Astrocyte" | 
                          merfishmeta$subclass == "Oligodendrocyte"| merfishmeta$subclass == "OPC"] <- "Non-Neuronal"

# Electrophysiology
setwd("/allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/ephys/data")
ephysfeats <- read_csv("Tolias_m1_patchseq_ephys_features.csv")
colnames(ephysfeats)[1] <- "Cell"

patchmeta <- read_csv("Tolias_m1_patchseq_meta_data.csv")


glut_col <- clupdate %>%
  filter(class_label == "Glutamatergic")
glut_col <- data.frame(glut_col,stringsAsFactors = FALSE)
glut_col$subclass_label <- as.character(glut_col$subclass_label)
glut_col$subclass_color <- as.character(glut_col$subclass_color)


toMatch <- c("ET","IT","CT","NP")

sub_metadata <- patchmeta %>%
  filter(grepl(paste(toMatch,collapse="|"), patchmeta$`RNA family`)==TRUE) 

glutephys <- merge(sub_metadata,ephysfeats, by = "Cell")
glutephys$`RNA family` <- gsub("ET","L5 ET", glutephys$`RNA family`)
glutephys$`RNA family`[glutephys$`RNA type` == "L5/6 NP CT"] <- "L6 CT"
glutephys$`RNA family` <- gsub("NP","L5/6 NP", glutephys$`RNA family`)
glutephys$`RNA family`[grepl("L5 IT",glutephys$`RNA type`)] <- "L5 IT"
glutephys$`RNA family`[grepl("L4/5 IT",glutephys$`RNA type`)] <- "L4/5 IT"
glutephys$`RNA family`[grepl("L2/3 IT",glutephys$`RNA type`)] <- "L2/3 IT"
glutephys$`RNA family`[grepl("L6b",glutephys$`RNA type`)] <- "L6b"
glutephys$`RNA family`[grepl("L6 CT",glutephys$`RNA type`)] <- "L6 CT"
glutephys$`RNA family`[grepl("L6 IT",glutephys$`RNA type`)] <- "L6 IT"


gaba_col <- clupdate %>%
  filter(class_label == "GABAergic")
gaba_col <- data.frame(gaba_col,stringsAsFactors = FALSE)
gaba_col$subclass_label <- as.character(gaba_col$subclass_label)
gaba_col$subclass_color <- as.character(gaba_col$subclass_color)


toMatch <- c("Pvalb","Vip","Sncg","Sst","Sst Chodl","Meis2","Lamp5")

gaba_metadata <- patchmeta %>%
  filter(grepl(paste(toMatch,collapse="|"), patchmeta$`RNA family`)==TRUE) 

gabaephys <- merge(gaba_metadata,ephysfeats, by = "Cell")
gabaephys$`RNA family`[gabaephys$`RNA type` == "Sst Chodl"] <- "Sst Chodl"


morphmet <- read_csv("Tolias_m1_patchseq_morph_features.csv")
colnames(morphmet)[1] <- "Cell"
sub_metadata <- rbind(sub_metadata,gaba_metadata)
allmorph <- merge(sub_metadata,morphmet,by="Cell")

allmorph$`RNA family` <- gsub("ET","L5 ET", allmorph$`RNA family`)
allmorph$`RNA family`[allmorph$`RNA type` == "L5/6 NP CT"] <- "L6 CT"
allmorph$`RNA family` <- gsub("NP","L5/6 NP", allmorph$`RNA family`)
allmorph$`RNA family`[grepl("L5 IT",allmorph$`RNA type`)] <- "L5 IT"
allmorph$`RNA family`[grepl("L4/5 IT",allmorph$`RNA type`)] <- "L4/5 IT"
allmorph$`RNA family`[grepl("L2/3 IT",allmorph$`RNA type`)] <- "L2/3 IT"
allmorph$`RNA family`[grepl("L6b",allmorph$`RNA type`)] <- "L6b"
allmorph$`RNA family`[grepl("L6 CT",allmorph$`RNA type`)] <- "L6 CT"
allmorph$`RNA family`[grepl("L6 IT",allmorph$`RNA type`)] <- "L6 IT"


# Getting everything set up
setwd("/allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr/")
nomenclature <- read_csv("Accesion_IDs_MOp.csv")

source("load_metadata.R")
cellsetnom <- read_csv("cellsetnomenclature.csv")
cellsetnom$`cell set preferred alias`[157] <- "Microglia-PVM"
cellsetnom$`cell set preferred alias`[140] <- "L4/5 IT"
cellsetnom$`cell set preferred alias`[152] <- "Oligodendrocyte"
cellsetnom$`cell set preferred alias`[150] <- "Astrocyte"
cellsetnom$`cell set preferred alias`[103] <- "Endothelial"
cellsetnom$`cell set preferred alias`[112] <- "Pericyte"
cellsetnom$`cell set preferred alias`[91] <- "OPC"

merfishmeta$subclass <- gsub("cytes","cyte",merfishmeta$subclass)
merfishmeta$subclass <- gsub("_"," ",merfishmeta$subclass)
merfishmeta$subclass <- gsub("23","2/3",merfishmeta$subclass)
merfishmeta$subclass <- gsub("45","4/5",merfishmeta$subclass)
merfishmeta$subclass <- gsub("56","5/6",merfishmeta$subclass)
merfishmeta$subclass <- gsub("Microglia|PVM","Microglia-PVM",merfishmeta$subclass)

alignedalias <- read_csv("alignedaliases.csv")
clupdate <- read_csv("cl_updated.df.csv")
nsfmarkers <- read_csv("NSForestMarkers_genesymbols.csv")
nsfmarkers <- nsfmarkers %>%
  select(-contains("ens"))

subnom <- cellsetnom[cellsetnom$`cell set preferred alias` %in% subnames,]
subnom <- subnom[c(1,2,4:25),]

glutsubnom <- subnom[subnom$`cell set preferred alias` %in% glutsubnames,]
glutsubalign <- alignedalias[alignedalias$`aligned aliases` %in% glutsubnames,]
glutsubnom <- glutsubnom[match(glutsubnames, glutsubnom$`cell set preferred alias`),]
glutsubalign <- glutsubalign[match(glutsubnames, glutsubalign$`aligned aliases`),]

metarep <- paste0(cellsetnom$species[1]," ",cellsetnom$`cell set structure`[1], " (", cellsetnom$`taxonomy id`[1], "), ")


gabasubnom <- subnom[subnom$`cell set preferred alias` %in% gabasubnames,]
gabasubalign <- alignedalias[alignedalias$`aligned aliases` %in% gabasubnames,]
gabasubnom <- gabasubnom[match(gabasubnames, gabasubnom$`cell set preferred alias`),]
gabasubalign <- gabasubalign[match(gabasubnames, gabasubalign$`aligned aliases`),]

nonneuralsubnom <- subnom[subnom$`cell set preferred alias` %in% nonneuralsubnames,]
nonneuralsubalign <- alignedalias[alignedalias$`aligned aliases` %in% nonneuralsubnames,]
nonneuralsubnom <- nonneuralsubnom[match(nonneuralsubnames, nonneuralsubnom$`cell set preferred alias`),]
nonneuralsubalign <- nonneuralsubalign[match(nonneuralsubnames, nonneuralsubalign$`aligned aliases`),]

nonneuronalsubnom <- subnom[subnom$`cell set preferred alias` %in% nonneuronalsubnames,]
nonneuronalsubalign <- alignedalias[alignedalias$`aligned aliases` %in% nonneuronalsubnames,]
nonneuronalsubnom <- nonneuronalsubnom[match(nonneuronalsubnames, nonneuronalsubnom$`cell set preferred alias`),]
nonneuronalsubalign <- nonneuronalsubalign[match(nonneuronalsubnames, nonneuronalsubalign$`aligned aliases`),]



######
# Loop through all taxons grouped by siblings. This is not the most elegant way to do this but made it easier to do direct
# comparisons with siblings. The file "cl.updated.df.csv" could probably be leveraged to get the sibling relationships and avoid 
# this method, reducing the amount of code necessary, but I was in a hurry. 

########
# Test loop
temp <- c()
for (i in 1:nrow(glutsubnom)) {
  
  # Pull taxon name and determine which data types are present
  taxon <- glutsubnom$`cell set preferred alias`[i]
  hasephys <- glutsubnom$hasEphys[glutsubnom$`cell set preferred alias` == taxon]
  hasmorph <- glutsubnom$hasMorpho[glutsubnom$`cell set preferred alias`== taxon]
  hasspatial <- glutsubnom$hasSpatial[glutsubnom$`cell set preferred alias`== taxon]
  
  # General info, genetic markers, and cluster metrics
  children <- unique(clupdate$cluster_label[clupdate$subclass_label == taxon])
  
  # Pull additional aliases for descriptions
  patchtypes <- NA
  if (!is.na(glutsubnom$`cell set  additional alias`[i]) == TRUE & glutsubnom$`cell set  additional alias`[i] != taxon) {
    patchtypes <- glutsubnom$`cell set  additional alias`[i]
  }
  
  # Pull NS-Forest markers, average genes detected, and proportion of cells in dataset
  markers <- as.character(nsfmarkers[nsfmarkers$Accession == glutsubnom$`cell set accession`[i],3:7])
  markers <- paste(markers[!is.na(markers) == TRUE], collapse = " and ")
  avggenes <- mean(snRNA_10xv3_B_met$genes_detected_label[snRNA_10xv3_B_met$subclass_label == taxon],na.rm=TRUE)
  avggenes <- round(avggenes)
  checkgenenums <- avggenes - avgall
  moreorless <- ifelse(checkgenenums < 0, " less than ", " more than ")
  avgall <- mean(snRNA_10xv3_B_met$genes_detected_label[snRNA_10xv3_B_met$class_label == "Glutamatergic"],na.rm=TRUE)
  avgall <- round(avgall)
  
  propcells <- sum(clupdate$cluster_size[clupdate$subclass_label == taxon])/sum(clupdate$cluster_size[clupdate$class_label == "Glutamatergic"])
  propcells <- round(propcells*100, 2)
  
  
  # Piece together text based on what data exists for each taxon
  if (is.na(patchtypes) == FALSE) {
  
  geninfo <- paste0(taxon, " is: ", glutsubalign$`aligned alias description`[i],
                        ". In ", metarep, taxon, " has ", 
                        length(children)," child types: ", paste(children[1:length(children)-1], collapse=", "),
                        " and ", children[length(children)], ". ", taxon,
                        " includes the Patch-seq and connectivity types ", paste(patchtypes, collapse = " and "),
                        ". ",
                        "The minimal set of markers required to distinguish this cell type from other cell types in the primary motor cortex is ",
                        markers, ", based on the NS Forest method of marker gene selection. ",
                        taxon, " cells represent ", propcells,
                        "% of all ", unique(clupdate$class_label[clupdate$subclass_label == taxon]),
                        " cells collected in the snRNA-seq 10X v3 B dataset. ", taxon, " cells have an average of ",
                        avggenes, " genes detected, which is ", abs(avggenes - avgall), moreorless, "the average for all Glutamatergic cells (",
                        avgall, "). ")
  }
  
  else {
    geninfo <- paste0(taxon, " is: ", glutsubalign$`aligned alias description`[i],
                      ". In ", metarep, taxon, " has ", 
                      length(children)," child types: ", paste(children[1:length(children)-1], collapse=", "),
                      " and ", children[length(children)], ". ",
                      "The minimal set of markers required to distinguish this cell type from other cell types in the primary motor cortex is ",
                      markers, ", based on the NS Forest method of marker gene selection. ",
                      taxon, " cells represent ", propcells,
                      "% of all ", unique(clupdate$class_label[clupdate$subclass_label == taxon]),
                      " cells collected in the snRNA-seq 10X v3 B dataset. ", taxon, " cells have an average of ",
                      avggenes, " genes detected, which is ", abs(avggenes - avgall), moreorless, "the average for all Glutamatergic cells (",
                      avgall, "). ")
  }
  
  tempspat <- ""
  
  if (hasspatial == "y") {
    avgdepth <- round(mean(merfishmeta$normalized_depth[merfishmeta$subclass == taxon], na.rm = TRUE),3)
    alldepth <- round(mean(merfishmeta$normalized_depth[merfishmeta$class_label == "Glutamatergic"], na.rm = TRUE), 3)
    
    spatlocinfo <- paste0("The average normalized cortical depth of ", taxon, " cells is ",
                          avgdepth, ", compared to an average of ", alldepth, " for all Glutamatergic cells. ")
    geninfo <- paste0(geninfo, spatlocinfo)
  }
  
  
  
  
  
  if (hasephys == "y") {
    
    firing <- round(mean(glutephys$`Max number of APs`[glutephys$`RNA family` == taxon],na.rm = TRUE)/0.6,2)
    inputres <- round(mean(glutephys$`Input resistance (MOhm)`[glutephys$`RNA family` == taxon],na.rm = TRUE),2)
    updown <- round(mean(glutephys$`Upstroke-to-downstroke ratio`[glutephys$`RNA family` == taxon],na.rm = TRUE),2)
    
    allfiring <- round(mean(glutephys$`Max number of APs`,na.rm = TRUE)/0.6,2)
    allinputres <- round(mean(glutephys$`Input resistance (MOhm)`,na.rm = TRUE),2)
    allupdown <- round(mean(glutephys$`Upstroke-to-downstroke ratio`,na.rm = TRUE),2)
    
    
    ephysinfo <- paste0(taxon, " cells have an average firing rate of ", firing,
                        " Hz, compared to ", allfiring, 
                        " Hz for all Glutamatergic cells recorded. Other distinguishing features include an input resistance of ",
                        inputres, " MOhms (", allinputres, " MOhms in all Glutamatergic cells), and upstroke:downstroke ratio of ",
                        updown, " (", allupdown, " for all Glutamatergic cells). ")
    geninfo <- paste0(geninfo,ephysinfo)
  }
  
  
  
  
  
  if (hasmorph == "y") {
    morphdepth <- round(mean(allmorph$`normalized depth`[allmorph$`RNA family` == taxon],na.rm = TRUE),2)
    axonlength <- round(mean(allmorph$`axon total length`[allmorph$`RNA family` == taxon],na.rm = TRUE),2)
    dendritelength <- round(mean(allmorph$`dendrite total length`[allmorph$`RNA family` == taxon],na.rm = TRUE),2)
    stems <- round(mean(allmorph$stems[allmorph$`RNA family` == taxon],na.rm = TRUE),2)
    
    
    allmorphdepth <- round(mean(allmorph$`normalized depth`[allmorph$`cell class`=="exc"],na.rm = TRUE),2)
    allaxonlength <- round(mean(allmorph$`axon total length`[allmorph$`cell class`=="exc"],na.rm = TRUE),2)
    alldendritelength <- round(mean(allmorph$`dendrite total length`[allmorph$`cell class`=="exc"],na.rm = TRUE),2)
    allstems <- round(mean(allmorph$stems[allmorph$`cell class`=="exc"],na.rm = TRUE),2)
    
    
    morphinfo <- paste0("The most distinguishing morpholgical features of ", taxon,
                        ", as determined by reconstructions from Patch-seq data, are normalized cortical depth (", morphdepth, " vs. ", allmorphdepth,
                        " average for all Glutamatergic cells), axon total length (", axonlength, " vs. ", allaxonlength, 
                        " microns for all Glutamatergic cells), dendrite total length (", dendritelength, " vs. ", alldendritelength,
                        " microns for all Glutamatergic cells), and stems (", stems, " vs. ", allstems, 
                        " for all Glutamatergic cells). ")
    
    geninfo <- paste0(geninfo, morphinfo)
  }
  
  
  
  print(geninfo)
  temp <- c(temp,geninfo)
  

  
  
  
}

glutdesc <- cbind.data.frame(glutsubnom$`cell set preferred alias`,glutsubnom$`cell set accession`,temp)
colnames(glutdesc) <- c("Alias", "Accession", "Summary")



temp <- c()

for (i in 1:nrow(gabasubnom)) {
  
  taxon <- gabasubnom$`cell set preferred alias`[i]
  hasephys <- gabasubnom$hasEphys[gabasubnom$`cell set preferred alias` == taxon]
  hasmorph <- gabasubnom$hasMorpho[gabasubnom$`cell set preferred alias`== taxon]
  hasspatial <- gabasubnom$hasSpatial[gabasubnom$`cell set preferred alias`== taxon]
  
  # General info, genetic markers, and cluster metrics
  children <- unique(clupdate$cluster_label[clupdate$subclass_label == taxon])
  
  patchtypes <- NA
  
  if (!is.na(gabasubnom$`cell set  additional alias`[i]) == TRUE & gabasubnom$`cell set  additional alias`[i] != taxon) {
    patchtypes <- gabasubnom$`cell set  additional alias`[i]
  }
  
  
  markers <- as.character(nsfmarkers[nsfmarkers$Accession == gabasubnom$`cell set accession`[i],3:7])
  markers <- paste(markers[!is.na(markers) == TRUE], collapse = " and ")
  avggenes <- mean(snRNA_10xv3_B_met$genes_detected_label[snRNA_10xv3_B_met$subclass_label == taxon],na.rm=TRUE)
  avggenes <- round(avggenes)
  checkgenenums <- avggenes - avgall
  moreorless <- ifelse(checkgenenums < 0, " less than ", " more than ")
  avgall <- mean(snRNA_10xv3_B_met$genes_detected_label[snRNA_10xv3_B_met$class_label == "GABAergic"],na.rm=TRUE)
  avgall <- round(avgall)
  
  propcells <- sum(clupdate$cluster_size[clupdate$subclass_label == taxon])/sum(clupdate$cluster_size[clupdate$class_label == "GABAergic"])
  propcells <- round(propcells*100, 2)
  
  
  if (is.na(patchtypes) == FALSE) {
    
    geninfo <- paste0(taxon, " is: ", gabasubalign$`aligned alias description`[i],
                      ". In ", metarep, taxon, " has ", 
                      length(children)," child types: ", paste(children[1:length(children)-1], collapse=", "),
                      " and ", children[length(children)], ". ", taxon,
                      " includes the Patch-seq and connectivity types ", paste(patchtypes, collapse = " and "),
                      ". ",
                      "The minimal set of markers required to distinguish this cell type from other cell types in the primary motor cortex is ",
                      markers, ", based on the NS Forest method of marker gene selection. ",
                      taxon, " cells represent ", propcells,
                      "% of all ", unique(clupdate$class_label[clupdate$subclass_label == taxon]),
                      " cells collected in the snRNA-seq 10X v3 B dataset. ", taxon, " cells have an average of ",
                      avggenes, " genes detected, which is ", abs(avggenes - avgall), moreorless, "the average for all GABAergic cells (",
                      avgall, "). ")
  }
  
  else {
    geninfo <- paste0(taxon, " is: ", gabasubalign$`aligned alias description`[i],
                      ". In ", metarep, taxon, " has ", 
                      length(children)," child types: ", paste(children[1:length(children)-1], collapse=", "),
                      " and ", children[length(children)], ". ",
                      "The minimal set of markers required to distinguish this cell type from other cell types in the primary motor cortex is ",
                      markers, ", based on the NS Forest method of marker gene selection. ",
                      taxon, " cells represent ", propcells,
                      "% of all ", unique(clupdate$class_label[clupdate$subclass_label == taxon]),
                      " cells collected in the snRNA-seq 10X v3 B dataset. ", taxon, " cells have an average of ",
                      avggenes, " genes detected, which is ", abs(avggenes - avgall), moreorless, "the average for all GABAergic cells (",
                      avgall, "). ")
  }
  
  tempspat <- ""
  
  if (hasspatial == "y") {
    avgdepth <- round(mean(merfishmeta$normalized_depth[merfishmeta$subclass == taxon], na.rm = TRUE),3)
    alldepth <- round(mean(merfishmeta$normalized_depth[merfishmeta$class_label == "GABAergic"], na.rm = TRUE), 3)
    
    spatlocinfo <- paste0("The average normalized cortical depth of ", taxon, " cells is ",
                          avgdepth, ", compared to an average of ", alldepth, " for all GABAergic cells. ")
    geninfo <- paste0(geninfo, spatlocinfo)
  }
  
  
  
  
  
  if (hasephys == "y") {
    
    firing <- round(mean(gabaephys$`Max number of APs`[gabaephys$`RNA family` == taxon],na.rm = TRUE)/0.6,2)
    inputres <- round(mean(gabaephys$`Input resistance (MOhm)`[gabaephys$`RNA family` == taxon],na.rm = TRUE),2)
    updown <- round(mean(gabaephys$`Upstroke-to-downstroke ratio`[gabaephys$`RNA family` == taxon],na.rm = TRUE),2)
    
    allfiring <- round(mean(gabaephys$`Max number of APs`,na.rm = TRUE)/0.6,2)
    allinputres <- round(mean(gabaephys$`Input resistance (MOhm)`,na.rm = TRUE),2)
    allupdown <- round(mean(gabaephys$`Upstroke-to-downstroke ratio`,na.rm = TRUE),2)
    
    
    ephysinfo <- paste0(taxon, " cells have an average firing rate of ", firing,
                        " Hz, compared to ", allfiring, 
                        " Hz for all GABAergic cells recorded. Other distinguishing features include an input resistance of ",
                        inputres, " MOhms (", allinputres, " MOhms in all GABAergic cells), and upstroke:downstroke ratio of ",
                        updown, " (", allupdown, " for all GABAergic cells). ")
    geninfo <- paste0(geninfo,ephysinfo)
  }
  
  
  
  
  
  if (hasmorph == "y") {
    morphdepth <- round(mean(allmorph$`normalized depth`[allmorph$`RNA family` == taxon],na.rm = TRUE),2)
    axonlength <- round(mean(allmorph$`axon total length`[allmorph$`RNA family` == taxon],na.rm = TRUE),2)
    dendritelength <- round(mean(allmorph$`dendrite total length`[allmorph$`RNA family` == taxon],na.rm = TRUE),2)
    stems <- round(mean(allmorph$stems[allmorph$`RNA family` == taxon],na.rm = TRUE),2)
    
    
    allmorphdepth <- round(mean(allmorph$`normalized depth`[allmorph$`cell class`=="inh"],na.rm = TRUE),2)
    allaxonlength <- round(mean(allmorph$`axon total length`[allmorph$`cell class`=="inh"],na.rm = TRUE),2)
    alldendritelength <- round(mean(allmorph$`dendrite total length`[allmorph$`cell class`=="inh"],na.rm = TRUE),2)
    allstems <- round(mean(allmorph$stems[allmorph$`cell class`=="exc"],na.rm = TRUE),2)
    
    
    morphinfo <- paste0("The most distinguishing morpholgical features of ", taxon,
                        ", as determined by reconstructions from Patch-seq data, are normalized cortical depth (", morphdepth, " vs. ", allmorphdepth,
                        " average for all GABAergic cells), axon total length (", axonlength, " vs. ", allaxonlength, 
                        " microns for all GABAergic cells), dendrite total length (", dendritelength, " vs. ", alldendritelength,
                        " microns for all GABAergic cells), and stems (", stems, " vs. ", allstems, 
                        " for all GABAergic cells). ")
    
    geninfo <- paste0(geninfo, morphinfo)
  }
  
  
  
  print(geninfo)
  temp <- c(temp,geninfo)
  
  
  
  
  
}

gabadesc <- cbind.data.frame(gabasubnom$`cell set preferred alias`,gabasubnom$`cell set accession`,temp)
colnames(gabadesc) <- c("Alias", "Accession", "Summary")





temp <- c()
for (i in 1:nrow(nonneuralsubnom)) {
  
  taxon <- nonneuralsubnom$`cell set preferred alias`[i]
  hasspatial <- nonneuralsubnom$hasSpatial[nonneuralsubnom$`cell set preferred alias`== taxon]
  
  # General info, genetic markers, and cluster metrics
  children <- unique(clupdate$cluster_label[clupdate$subclass_label == taxon])
  
  patchtypes <- NA
  
  if (!is.na(nonneuralsubnom$`cell set  additional alias`[i]) == TRUE & nonneuralsubnom$`cell set  additional alias`[i] != taxon) {
    patchtypes <- nonneuralsubnom$`cell set  additional alias`[i]
  }
  
  
  markers <- as.character(nsfmarkers[nsfmarkers$Accession == nonneuralsubnom$`cell set accession`[i],3:7])
  markers <- paste(markers[!is.na(markers) == TRUE], collapse = " and ")
  avggenes <- mean(snRNA_10xv3_B_met$genes_detected_label[snRNA_10xv3_B_met$subclass_label == taxon],na.rm=TRUE)
  avggenes <- round(avggenes)
  
  avgall <- mean(snRNA_10xv3_B_met$genes_detected_label[snRNA_10xv3_B_met$class_label == "Non-Neural"],na.rm=TRUE)
  avgall <- round(avgall)
  checkgenenums <- avggenes - avgall
  moreorless <- ifelse(checkgenenums < 0, " less than ", " more than ")
  propcells <- sum(clupdate$cluster_size[clupdate$subclass_label == taxon])/sum(clupdate$cluster_size[clupdate$class_label == "Non-Neural"])
  propcells <- round(propcells*100, 2)
  
  
  if (is.na(patchtypes) == FALSE) {
    
    geninfo <- paste0(taxon, " is: ", nonneuralsubalign$`aligned alias description`[i],
                      ". In ", metarep, taxon, " has ", 
                      length(children)," child types: ", paste(children[1:length(children)-1], collapse=", "),
                      " and ", children[length(children)], ". ", taxon,
                      " includes the Patch-seq and connectivity types ", paste(patchtypes, collapse = " and "),
                      ". ",
                      "The minimal set of markers required to distinguish this cell type from other cell types in the primary motor cortex is ",
                      markers, ", based on the NS Forest method of marker gene selection. ",
                      taxon, " cells represent ", propcells,
                      "% of all ", unique(clupdate$class_label[clupdate$subclass_label == taxon]),
                      " cells collected in the snRNA-seq 10X v3 B dataset. ", taxon, " cells have an average of ",
                      avggenes, " genes detected, which is ", abs(avggenes - avgall), " more than the average for all Non-Neural cells (",
                      avgall, "). ")
  }
  
  else {
    geninfo <- paste0(taxon, " is: ", nonneuralsubalign$`aligned alias description`[i],
                      ". In ", metarep, taxon, " has ", 
                      length(children)," child types: ", paste(children[1:length(children)-1], collapse=", "),
                      " and ", children[length(children)], ". ",
                      "The minimal set of markers required to distinguish this cell type from other cell types in the primary motor cortex is ",
                      markers, ", based on the NS Forest method of marker gene selection. ",
                      taxon, " cells represent ", propcells,
                      "% of all ", unique(clupdate$class_label[clupdate$subclass_label == taxon]),
                      " cells collected in the snRNA-seq 10X v3 B dataset. ", taxon, " cells have an average of ",
                      avggenes, " genes detected, which is ", abs(avggenes - avgall), moreorless, "the average for all Non-Neural cells (",
                      avgall, "). ")
  }
  
  tempspat <- ""
  
  if (hasspatial == "y") {
    avgdepth <- round(mean(merfishmeta$normalized_depth[merfishmeta$subclass == taxon], na.rm = TRUE),3)
    alldepth <- round(mean(merfishmeta$normalized_depth[merfishmeta$class_label == "Non-Neural"], na.rm = TRUE), 3)
    
    spatlocinfo <- paste0("The average normalized cortical depth of ", taxon, " cells is ",
                          avgdepth, ", compared to an average of ", alldepth, " for all Non-Neural cells. ")
    geninfo <- paste0(geninfo, spatlocinfo)
  }
  
  
  

  
  
  
  print(geninfo)
  temp <- c(temp,geninfo)
  
  
  
  
}

nonneuraldesc <- cbind.data.frame(nonneuralsubnom$`cell set preferred alias`,nonneuralsubnom$`cell set accession`,temp)
colnames(nonneuraldesc) <- c("Alias", "Accession", "Summary")








temp <- c()
for (i in 1:nrow(nonneuronalsubnom)) {
  
  taxon <- nonneuronalsubnom$`cell set preferred alias`[i]
  hasspatial <- nonneuronalsubnom$hasSpatial[nonneuronalsubnom$`cell set preferred alias`== taxon]
  
  # General info, genetic markers, and cluster metrics
  children <- unique(clupdate$cluster_label[clupdate$subclass_label == taxon])
  
  
  patchtypes <- NA
  if (!is.na(nonneuronalsubnom$`cell set  additional alias`[i]) == TRUE & nonneuronalsubnom$`cell set  additional alias`[i] != taxon) {
    patchtypes <- nonneuronalsubnom$`cell set  additional alias`[i]
  }
  
  
  markers <- as.character(nsfmarkers[nsfmarkers$Accession == nonneuronalsubnom$`cell set accession`[i],3:7])
  markers <- paste(markers[!is.na(markers) == TRUE], collapse = " and ")
  avggenes <- mean(snRNA_10xv3_B_met$genes_detected_label[snRNA_10xv3_B_met$subclass_label == taxon],na.rm=TRUE)
  avggenes <- round(avggenes)
  
  avgall <- mean(snRNA_10xv3_B_met$genes_detected_label[snRNA_10xv3_B_met$class_label == "Non-Neuronal"],na.rm=TRUE)
  avgall <- round(avgall)
  checkgenenums <- avggenes - avgall
  moreorless <- ifelse(checkgenenums < 0, " less than ", " more than ")
  
  propcells <- sum(clupdate$cluster_size[clupdate$subclass_label == taxon])/sum(clupdate$cluster_size[clupdate$class_label == "Non-Neuronal"])
  propcells <- round(propcells*100, 2)
  
  
  if (is.na(patchtypes) == FALSE) {
    
    geninfo <- paste0(taxon, " is: ", nonneuronalsubalign$`aligned alias description`[i],
                      ". In ", metarep, taxon, " has ", 
                      length(children)," child types: ", paste(children[1:length(children)-1], collapse=", "),
                      " and ", children[length(children)], ". ", taxon,
                      " includes the Patch-seq and connectivity types ", paste(patchtypes, collapse = " and "),
                      ". ",
                      "The minimal set of markers required to distinguish this cell type from other cell types in the primary motor cortex is ",
                      markers, ", based on the NS Forest method of marker gene selection. ",
                      taxon, " cells represent ", propcells,
                      "% of all ", unique(clupdate$class_label[clupdate$subclass_label == taxon]),
                      " cells collected in the snRNA-seq 10X v3 B dataset. ", taxon, " cells have an average of ",
                      avggenes, " genes detected, which is ", abs(avggenes - avgall), moreorless, "the average for all Non-Neuronal cells (",
                      avgall, "). ")
  }
  
  else {
    geninfo <- paste0(taxon, " is: ", nonneuronalsubalign$`aligned alias description`[i],
                      ". In ", metarep, taxon, " has ", 
                      length(children)," child types: ", paste(children[1:length(children)-1], collapse=", "),
                      " and ", children[length(children)], ". ",
                      "The minimal set of markers required to distinguish this cell type from other cell types in the primary motor cortex is ",
                      markers, ", based on the NS Forest method of marker gene selection. ",
                      taxon, " cells represent ", propcells,
                      "% of all ", unique(clupdate$class_label[clupdate$subclass_label == taxon]),
                      " cells collected in the snRNA-seq 10X v3 B dataset. ", taxon, " cells have an average of ",
                      avggenes, " genes detected, which is ", abs(avggenes - avgall), moreorless, "the average for all Non-Neuronal cells (",
                      avgall, "). ")
  }
  
  tempspat <- ""
  
  if (hasspatial == "y") {
    avgdepth <- round(mean(merfishmeta$normalized_depth[merfishmeta$subclass == taxon], na.rm = TRUE),3)
    alldepth <- round(mean(merfishmeta$normalized_depth[merfishmeta$class_label == "Non-Neuronal"], na.rm = TRUE), 3)
    
    spatlocinfo <- paste0("The average normalized cortical depth of ", taxon, " cells is ",
                          avgdepth, ", compared to an average of ", alldepth, " for all Non-Neuronal cells. ")
    geninfo <- paste0(geninfo, spatlocinfo)
  }
  
  
  
  
  
  
  
  print(geninfo)
  temp <- c(temp,geninfo)
  
  
 
  
  
}

nonneuronaldesc <- cbind.data.frame(nonneuronalsubnom$`cell set preferred alias`,nonneuronalsubnom$`cell set accession`,temp)
colnames(nonneuronaldesc) <- c("Alias", "Accession", "Summary")


subdesc <- rbind(glutdesc, gabadesc, nonneuraldesc, nonneuronaldesc)

######
# After completing all descriptions, you can write to CSV or convert to JSON (preferred format of product team) as is done below
######
# convert to JSON

setwd("/allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr")

finaldesc <- read_csv("All Descriptions.csv")
finaldesc <- finaldesc[,2:5]


test <- finaldesc[,2:4]
test <- split(test, seq(nrow(finaldesc)))
names(test) <- finaldesc$Accession
descjson <- toJSON(x = test, dataframe = 'columns', pretty = T)

write(descjson,"All Descriptions.json")
