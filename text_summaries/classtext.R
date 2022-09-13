# Doing classs now

cellsetnom$`cell set preferred alias`[153] <- "Non-Neural"
classnom <- cellsetnom[cellsetnom$`cell set preferred alias` %in% classnames,]

alignedalias$`aligned aliases`[30:33] <- c("Non-Neuronal", "GABAergic", "Glutamatergic", "Non-Neural")

classalign <- alignedalias[alignedalias$`aligned aliases` %in% classnames,]



classalign <- classalign[match(classnames, classnom$`cell set preferred alias`),]
classalign <- classalign[c(2,4,3,1),]



temp <- c()
for (i in 1:nrow(classnom)) {
  
  taxon <- classnom$`cell set preferred alias`[i]
  hasephys <- classnom$hasEphys[classnom$`cell set preferred alias` == taxon]
  hasmorph <- classnom$hasMorpho[classnom$`cell set preferred alias`== taxon]
  hasspatial <- classnom$hasSpatial[classnom$`cell set preferred alias`== taxon]
  alignedaliasdesc <- classalign$`aligned alias description`[i]
  
  # General info, genetic markers, and cluster metrics

  
  markers <- as.character(nsfmarkers[nsfmarkers$Accession == classnom$`cell set accession`[i],3:7])
  markers <- paste(markers[!is.na(markers) == TRUE], collapse = " and ")
  avggenes <- mean(snRNA_10xv3_B_met$genes_detected_label[snRNA_10xv3_B_met$class_label == taxon],na.rm=TRUE)
  avggenes <- round(avggenes)
  checkgenenums <- avggenes - avgall
  moreorless <- ifelse(checkgenenums < 0, " less than ", " more than ")
  avgall <- mean(snRNA_10xv3_B_met$genes_detected_label,na.rm=TRUE)
  avgall <- round(avgall)
  
  propcells <- sum(clupdate$cluster_size[clupdate$class_label == taxon])/sum(clupdate$cluster_size)
  propcells <- round(propcells*100, 2)
  
  
  geninfo <- paste0(taxon, " is: ", classalign$`aligned alias description`[i], ". ",
                      "The minimal set of markers required to distinguish this cell class from other cell classs in the primary motor cortex is ",
                      markers, ", based on the NS Forest method of marker gene selection. ",
                      taxon, " cells represent ", propcells,
                      "% of all cells collected in the snRNA-seq 10X v3 B dataset. ", taxon, " cells have an average of ",
                      avggenes, " genes detected, which is ", abs(avggenes - avgall), moreorless, "the average for allcells (",
                      avgall, "). ")
  
  
  
  
  if (hasspatial == "y") {
    avgdepth <- round(mean(merfishmeta$normalized_depth[merfishmeta$class_label == taxon], na.rm = TRUE),3)
    alldepth <- round(mean(merfishmeta$normalized_depth, na.rm = TRUE), 3)
    
    spatlocinfo <- paste0("The average normalized cortical depth of ", taxon, " cells is ",
                          avgdepth, ", compared to an average of ", alldepth, " for all cells. ")
    geninfo <- paste0(geninfo, spatlocinfo)
  }
  
  
  
  
  
  if (hasephys == "y") {
    
    firing <- round(mean(combineephys$`Max number of APs`[combineephys$`RNA family` %in% unique(clupdate$subclass_label[clupdate$class_label == taxon])],
                         na.rm = TRUE)/0.6,2)
    inputres <- round(mean(combineephys$`Input resistance (MOhm)`[combineephys$`RNA family` %in% unique(clupdate$subclass_label[clupdate$class_label == taxon])],na.rm = TRUE),2)
    updown <- round(mean(combineephys$`Upstroke-to-downstroke ratio`[combineephys$`RNA family` %in% unique(clupdate$subclass_label[clupdate$class_label == taxon])],na.rm = TRUE),2)
    
    allfiring <- round(mean(combineephys$`Max number of APs`,na.rm = TRUE)/0.6,2)
    allinputres <- round(mean(combineephys$`Input resistance (MOhm)`,na.rm = TRUE),2)
    allupdown <- round(mean(combineephys$`Upstroke-to-downstroke ratio`,na.rm = TRUE),2)
    
    
    ephysinfo <- paste0(taxon, " cells have an average firing rate of ", firing,
                        " Hz, compared to ", allfiring, 
                        " Hz for all neurons recorded. Other distinguishing features include an input resistance of ",
                        inputres, " MOhms (", allinputres, " MOhms in all neurons), and upstroke:downstroke ratio of ",
                        updown, " (", allupdown, " for all neurons). ")
    geninfo <- paste0(geninfo,ephysinfo)
  }
  
  
  
  
  
  if (hasmorph == "y") {
    morphdepth <- round(mean(allmorph$`normalized depth`[allmorph$`RNA family` %in% unique(clupdate$subclass_label[clupdate$class_label == taxon])],na.rm = TRUE),2)
    axonlength <- round(mean(allmorph$`axon total length`[allmorph$`RNA family` %in% unique(clupdate$subclass_label[clupdate$class_label == taxon])],na.rm = TRUE),2)
    dendritelength <- round(mean(allmorph$`dendrite total length`[allmorph$`RNA family` %in% unique(clupdate$subclass_label[clupdate$class_label == taxon])],na.rm = TRUE),2)
    stems <- round(mean(allmorph$stems[allmorph$`RNA family` %in% unique(clupdate$subclass_label[clupdate$class_label == taxon])],na.rm = TRUE),2)
    
    
    allmorphdepth <- round(mean(allmorph$`normalized depth`,na.rm = TRUE),2)
    allaxonlength <- round(mean(allmorph$`axon total length`,na.rm = TRUE),2)
    alldendritelength <- round(mean(allmorph$`dendrite total length`,na.rm = TRUE),2)
    allstems <- round(mean(allmorph$stems,na.rm = TRUE),2)
    
    
    morphinfo <- paste0("The most distinguishing morpholgical features of ", taxon,
                        ", as determined by reconstructions from Patch-seq data, are normalized cortical depth (", morphdepth, " vs. ", allmorphdepth,
                        " average for all neurons), axon total length (", axonlength, " vs. ", allaxonlength, 
                        " microns for all neurons), dendrite total length (", dendritelength, " vs. ", alldendritelength,
                        " microns for all neurons), and stems (", stems, " vs. ", allstems, 
                        " for all neurons). ")
    
    geninfo <- paste0(geninfo, morphinfo)
  }
  
  
  
  print(geninfo)
  temp <- c(temp,geninfo)
  
  
  
  
  
}

classdesc <- cbind.data.frame(classnom$`cell set preferred alias`,classnom$`cell set accession`,temp)
colnames(classdesc) <- c("Alias", "Accession", "Summary")
