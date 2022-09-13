# Doing types now


typenom <- cellsetnom[cellsetnom$`cell set preferred alias` %in% typenames,]
typenom <- typenom[c(1:27,29:68,70:108,110:113),]

typealign <- alignedalias[alignedalias$`aligned aliases` %in% typenames,]



typealign <- typealign[match(typenames, typenom$`cell set preferred alias`),]

combineephys <- rbind(combineephys,gabaephys)

typenom <- typenom[57:61,]
temp <- c()
for (i in 1:nrow(typenom)) {
  
  taxon <- typenom$`cell set preferred alias`[i]
  hasephys <- typenom$hasEphys[typenom$`cell set preferred alias` == taxon]
  hasmorph <- typenom$hasMorpho[typenom$`cell set preferred alias`== taxon]
  hasspatial <- typenom$hasSpatial[typenom$`cell set preferred alias`== taxon]
  alignedaliasdesc <- typealign$`aligned alias description`[i]
  
  # General info, genetic markers, and cluster metrics
  parent <- typenom$parent[i]
  
  patchtypes <- NA
  
  if (!is.na(typenom$`cell set  additional alias`[i]) == TRUE & typenom$`cell set  additional alias`[i] != taxon) {
    patchtypes <- typenom$`cell set  additional alias`[i]
  }
  
  
  markers <- as.character(nsfmarkers[nsfmarkers$Accession == typenom$`cell set accession`[i],3:7])
  markers <- paste(markers[!is.na(markers) == TRUE], collapse = " and ")
  avggenes <- mean(snRNA_10xv3_B_met$genes_detected_label[snRNA_10xv3_B_met$cluster_label == taxon],na.rm=TRUE)
  avggenes <- round(avggenes)
  checkgenenums <- avggenes - avgall
  moreorless <- ifelse(checkgenenums < 0, " less than ", " more than ")
  avgall <- mean(snRNA_10xv3_B_met$genes_detected_label[snRNA_10xv3_B_met$subclass_label == parent],na.rm=TRUE)
  avgall <- round(avgall)
  
  propcells <- sum(clupdate$cluster_size[clupdate$cluster_label == taxon])/sum(clupdate$cluster_size[clupdate$subclass_label == parent])
  propcells <- round(propcells*100, 2)
  
  
  if (is.na(patchtypes) == FALSE) {
    
    geninfo <- paste0("In ", metarep, taxon, " is a member of the ", 
                      parent," subclass. ", taxon,
                      " includes the Patch-seq and connectivity types ", paste(patchtypes, collapse = " and "),
                      ". ",
                      "The minimal set of markers required to distinguish this cell type from other cell types in the primary motor cortex is ",
                      markers, ", based on the NS Forest method of marker gene selection. ",
                      taxon, " cells represent ", propcells,
                      "% of all ", unique(clupdate$subclass_label[clupdate$cluster_label == taxon]),
                      " cells collected in the snRNA-seq 10X v3 B dataset. ", taxon, " cells have an average of ",
                      avggenes, " genes detected, which is ", abs(avggenes - avgall), moreorless, "the average for all ",parent ," cells (",
                      avgall, "). ")
  }
  
  else {
    geninfo <- paste0("In ", metarep, taxon, " is a member of the ", 
                      parent," subclass. ",
                      "The minimal set of markers required to distinguish this cell type from other cell types in the primary motor cortex is ",
                      markers, ", based on the NS Forest method of marker gene selection. ",
                      taxon, " cells represent ", propcells,
                      "% of all ", unique(clupdate$subclass_label[clupdate$cluster_label == taxon]),
                      " cells collected in the snRNA-seq 10X v3 B dataset. ", taxon, " cells have an average of ",
                      avggenes, " genes detected, which is ", abs(avggenes - avgall), moreorless, "the average for all ",parent ," cells (",
                      avgall, "). ")
  }
  

  
  if (hasspatial == "y") {
    avgdepth <- round(mean(merfishmeta$normalized_depth[merfishmeta$label == taxon], na.rm = TRUE),3)
    alldepth <- round(mean(merfishmeta$normalized_depth[merfishmeta$subclass == parent], na.rm = TRUE), 3)
    
    spatlocinfo <- paste0("The average normalized cortical depth of ", taxon, " cells is ",
                          avgdepth, ", compared to an average of ", alldepth, " for all ", parent, " cells. ")
    geninfo <- paste0(geninfo, spatlocinfo)
  }
  
  
  
  
  
  if (hasephys == "y") {
    
    firing <- round(mean(combineephys$`Max number of APs`[combineephys$`RNA type` == taxon],na.rm = TRUE)/0.6,2)
    inputres <- round(mean(combineephys$`Input resistance (MOhm)`[combineephys$`RNA type` == taxon],na.rm = TRUE),2)
    updown <- round(mean(combineephys$`Upstroke-to-downstroke ratio`[combineephys$`RNA type` == taxon],na.rm = TRUE),2)
    
    allfiring <- round(mean(combineephys$`Max number of APs`[combineephys$`RNA family` == parent],na.rm = TRUE)/0.6,2)
    allinputres <- round(mean(combineephys$`Input resistance (MOhm)`[combineephys$`RNA family` == parent],na.rm = TRUE),2)
    allupdown <- round(mean(combineephys$`Upstroke-to-downstroke ratio`[combineephys$`RNA family` == parent],na.rm = TRUE),2)
    
    
    ephysinfo <- paste0(taxon, " cells have an average firing rate of ", firing,
                        " Hz, compared to ", allfiring, 
                        " Hz for all ", parent, " cells recorded. Other distinguishing features include an input resistance of ",
                        inputres, " MOhms (", allinputres, " MOhms in all ", parent, " cells), and upstroke:downstroke ratio of ",
                        updown, " (", allupdown, " for all ", parent, " cells). ")
    geninfo <- paste0(geninfo,ephysinfo)
  }
  
  
  
  
  
  if (hasmorph == "y") {
    morphdepth <- round(mean(allmorph$`normalized depth`[allmorph$`RNA type` == taxon],na.rm = TRUE),2)
    axonlength <- round(mean(allmorph$`axon total length`[allmorph$`RNA type` == taxon],na.rm = TRUE),2)
    dendritelength <- round(mean(allmorph$`dendrite total length`[allmorph$`RNA type` == taxon],na.rm = TRUE),2)
    stems <- round(mean(allmorph$stems[allmorph$`RNA type` == taxon],na.rm = TRUE),2)
    
    
    allmorphdepth <- round(mean(allmorph$`normalized depth`[allmorph$`RNA family` == parent],na.rm = TRUE),2)
    allaxonlength <- round(mean(allmorph$`axon total length`[allmorph$`RNA family` == parent],na.rm = TRUE),2)
    alldendritelength <- round(mean(allmorph$`dendrite total length`[allmorph$`RNA family` == parent],na.rm = TRUE),2)
    allstems <- round(mean(allmorph$stems[allmorph$`RNA family` == parent],na.rm = TRUE),2)
    
    
    morphinfo <- paste0("The most distinguishing morpholgical features of ", taxon,
                        ", as determined by reconstructions from Patch-seq data, are normalized cortical depth (", morphdepth, " vs. ", allmorphdepth,
                        " average for all ", parent, " cells), axon total length (", axonlength, " vs. ", allaxonlength, 
                        " microns for all ", parent, " cells), dendrite total length (", dendritelength, " vs. ", alldendritelength,
                        " microns for all ", parent, " cells), and stems (", stems, " vs. ", allstems, 
                        " for all ", parent, " cells). ")
    
    geninfo <- paste0(geninfo, morphinfo)
  }
  
  
  
  print(geninfo)
  temp <- c(temp,geninfo)
  
  
  
  
  
}

typedesc <- cbind.data.frame(typenom$`cell set preferred alias`,typenom$`cell set accession`,temp)
colnames(typedesc) <- c("Alias", "Accession", "Summary")
