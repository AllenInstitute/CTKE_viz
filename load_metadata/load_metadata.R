setwd("//allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr")
dir <- "\\\\allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr"
allmet <- read_csv("final_met_all.csv")

clupdate <- read_csv("cl_updated.df.csv")
genelist <- read_csv("MOP_genelist.csv")
genelist$SMC[9] <- "Kcnj8"
genelist <- genelist[1:10,]


scRNA_10xv2_A_met <- allmet %>% filter(platform_label %in% "scRNA 10X v2 A")
scRNA_10xv2_A_met$sample_name <- gsub("10X_cells_v2_AIBS.","",scRNA_10xv2_A_met$sample_name)
scRNA_10xv2_A_met$cluster_id[scRNA_10xv2_A_met$cluster_label == "Sst Pvalb Calb2"] <- 41
scRNA_10xv2_A_met$cluster_id[scRNA_10xv2_A_met$cluster_label == "Astro Aqp4_Gfap"] <- 93
scRNA_10xv2_A_met$cluster_id[scRNA_10xv2_A_met$cluster_label == "Astro_Top2a"] <- 92
scRNA_10xv2_A_met$cluster_id[scRNA_10xv2_A_met$cluster_label == "Astro Aqp4_Slc7a10"] <- 94

scRNA_10xv3_A_met <- allmet %>% filter(platform_label %in% "scRNA 10X v3 A")
scRNA_10xv3_A_met$sample_name <- gsub("10X_cells_v3_AIBS.","",scRNA_10xv3_A_met$sample_name)
scRNA_10xv3_A_met$cluster_id[scRNA_10xv3_A_met$cluster_label == "Sst Pvalb Calb2"] <- 41
scRNA_10xv3_A_met$cluster_id[scRNA_10xv3_A_met$cluster_label == "Sst Calb2"] <- 39
scRNA_10xv3_A_met$cluster_id[scRNA_10xv3_A_met$cluster_label == "L5/6 NP CT"] <- 76
scRNA_10xv3_A_met$cluster_id[scRNA_10xv3_A_met$cluster_label == "L6 CT Gpr139"] <- 77
scRNA_10xv3_A_met$cluster_id[scRNA_10xv3_A_met$cluster_label == "L6 CT Cpa6"] <- 78
scRNA_10xv3_A_met$cluster_id[scRNA_10xv3_A_met$cluster_label == "L6 CT Grp"] <- 79
scRNA_10xv3_A_met$cluster_id[scRNA_10xv3_A_met$cluster_label == "L6 CT Pou3f2"] <- 80
scRNA_10xv3_A_met$cluster_id[scRNA_10xv3_A_met$cluster_label == "L6 CT Kit_1"] <- 81
scRNA_10xv3_A_met$cluster_id[scRNA_10xv3_A_met$cluster_label == "L6 CT Kit_2"] <- 82
scRNA_10xv3_A_met$cluster_id[scRNA_10xv3_A_met$cluster_label == "Astro_Top2a"] <- 92
scRNA_10xv3_A_met$cluster_id[scRNA_10xv3_A_met$cluster_label == "Astro Aqp4_Gfap"] <- 93
scRNA_10xv3_A_met$cluster_id[scRNA_10xv3_A_met$cluster_label == "Astro Aqp4_Slc7a10"] <- 94

snRNA_10xv2_A_met <- allmet %>% filter(platform_label %in% "snRNA 10X v2 A")
snRNA_10xv2_A_met$sample_name <- gsub("10X_nuclei_v2_AIBS.","",snRNA_10xv2_A_met$sample_name)
snRNA_10xv2_A_met$cluster_id[snRNA_10xv2_A_met$cluster_label == "Sst Pvalb Calb2"] <- 41
snRNA_10xv2_A_met$cluster_id[snRNA_10xv2_A_met$cluster_label == "Sst Calb2"] <- 39
snRNA_10xv2_A_met$cluster_id[snRNA_10xv2_A_met$cluster_label == "Astro_Top2a"] <- 92
snRNA_10xv2_A_met$cluster_id[snRNA_10xv2_A_met$cluster_label == "Astro Aqp4_Gfap"] <- 93
snRNA_10xv2_A_met$cluster_id[snRNA_10xv2_A_met$cluster_label == "Astro Aqp4_Slc7a10"] <- 94

snRNA_10xv3_A_met <- allmet %>% filter(platform_label %in% "snRNA 10X v3 A")
snRNA_10xv3_A_met$sample_name <- gsub("10X_nuclei_v3_AIBS.","",snRNA_10xv3_A_met$sample_name)
snRNA_10xv3_A_met$cluster_id[snRNA_10xv3_A_met$cluster_label == "Sst Pvalb Calb2"] <- 41
snRNA_10xv3_A_met$cluster_id[snRNA_10xv3_A_met$cluster_label == "Sst Calb2"] <- 39
snRNA_10xv3_A_met$cluster_id[snRNA_10xv3_A_met$cluster_label == "Meis2"] <- 58
snRNA_10xv3_A_met$cluster_id[snRNA_10xv3_A_met$cluster_label == "L5/6 NP CT"] <- 76
snRNA_10xv3_A_met$cluster_id[snRNA_10xv3_A_met$cluster_label == "L6 CT Gpr139"] <- 77
snRNA_10xv3_A_met$cluster_id[snRNA_10xv3_A_met$cluster_label == "L6 CT Cpa6"] <- 78
snRNA_10xv3_A_met$cluster_id[snRNA_10xv3_A_met$cluster_label == "L6 CT Grp"] <- 79
snRNA_10xv3_A_met$cluster_id[snRNA_10xv3_A_met$cluster_label == "L6 CT Pou3f2"] <- 80
snRNA_10xv3_A_met$cluster_id[snRNA_10xv3_A_met$cluster_label == "L6 CT Kit_1"] <- 81
snRNA_10xv3_A_met$cluster_id[snRNA_10xv3_A_met$cluster_label == "L6 CT Kit_2"] <- 82
snRNA_10xv3_A_met$cluster_id[snRNA_10xv3_A_met$cluster_label == "Astro_Top2a"] <- 92
snRNA_10xv3_A_met$cluster_id[snRNA_10xv3_A_met$cluster_label == "Astro Aqp4_Gfap"] <- 93
snRNA_10xv3_A_met$cluster_id[snRNA_10xv3_A_met$cluster_label == "Astro Aqp4_Slc7a10"] <- 94

snRNA_10xv3_B_met <- allmet %>% filter(platform_label %in% "snRNA 10X v3 B")
snRNA_10xv3_B_met$sample_name <- gsub("10X_nuclei_v3_Broad.","",snRNA_10xv3_B_met$sample_name)


smartcellmet <- allmet %>% filter(platform_label %in% "scRNA SMART")
smartcellmet$sample_name <- gsub("SmartSeq_cells_AIBS.","",smartcellmet$sample_name)
smartcellmet$sample_name <- gsub("-",".",smartcellmet$sample_name)
smartcellmet$cluster_id[smartcellmet$cluster_label == "Sst Pvalb Calb2"] <- 41
smartcellmet$cluster_id[smartcellmet$cluster_label == "L2/3 IT_1"] <- 60
smartcellmet$cluster_id[smartcellmet$cluster_label == "Sst Calb2"] <- 39
smartcellmet$cluster_id[smartcellmet$cluster_label == "L2/3 IT_2"] <- 61
smartcellmet$cluster_id[smartcellmet$cluster_label == "L2/3 IT_3"] <- 62
smartcellmet$cluster_id[smartcellmet$cluster_label == "L5 ET_1"] <- 72
smartcellmet$cluster_id[smartcellmet$cluster_label == "L5 ET_2"] <- 73
smartcellmet$cluster_id[smartcellmet$cluster_label == "L5 ET_3"] <- 74
smartcellmet$cluster_id[smartcellmet$cluster_label == "L5 ET_4"] <- 75
smartcellmet$cluster_id[smartcellmet$cluster_label == "L5/6 NP CT"] <- 76
smartcellmet$cluster_id[smartcellmet$cluster_label == "L6 CT Gpr139"] <- 77
smartcellmet$cluster_id[smartcellmet$cluster_label == "L6 CT Cpa6"] <- 78
smartcellmet$cluster_id[smartcellmet$cluster_label == "L6 CT Grp"] <- 79
smartcellmet$cluster_id[smartcellmet$cluster_label == "L6 CT Pou3f2"] <- 80
smartcellmet$cluster_id[smartcellmet$cluster_label == "L6 CT Kit_1"] <- 81
smartcellmet$cluster_id[smartcellmet$cluster_label == "L6 CT Kit_2"] <- 82
smartcellmet$cluster_id[smartcellmet$cluster_label == "Astro_Top2a"] <- 92
smartcellmet$cluster_id[smartcellmet$cluster_label == "Astro Aqp4_Gfap"] <- 93
smartcellmet$cluster_id[smartcellmet$cluster_label == "Astro Aqp4_Slc7a10"] <- 94

smartnucleimet <- allmet %>% filter(platform_label %in% "snRNA SMART")
smartnucleimet$sample_name <- gsub("SmartSeq_nuclei_AIBS.","",smartnucleimet$sample_name)
smartnucleimet$sample_name <- gsub("-",".",smartnucleimet$sample_name)
smartnucleimet$cluster_id[smartnucleimet$cluster_label == "Sst Pvalb Calb2"] <- 41
smartnucleimet$cluster_id[smartnucleimet$cluster_label == "L2/3 IT_1"] <- 60
smartnucleimet$cluster_id[smartnucleimet$cluster_label == "Sst Calb2"] <- 39
smartnucleimet$cluster_id[smartnucleimet$cluster_label == "L2/3 IT_2"] <- 61
smartnucleimet$cluster_id[smartnucleimet$cluster_label == "L2/3 IT_3"] <- 62
smartnucleimet$cluster_id[smartnucleimet$cluster_label == "L5 ET_1"] <- 72
smartnucleimet$cluster_id[smartnucleimet$cluster_label == "L5 ET_2"] <- 73
smartnucleimet$cluster_id[smartnucleimet$cluster_label == "L5 ET_3"] <- 74
smartnucleimet$cluster_id[smartnucleimet$cluster_label == "L5 ET_4"] <- 75
smartnucleimet$cluster_id[smartnucleimet$cluster_label == "L5/6 NP CT"] <- 76
smartnucleimet$cluster_id[smartnucleimet$cluster_label == "L6 CT Gpr139"] <- 77
smartnucleimet$cluster_id[smartnucleimet$cluster_label == "L6 CT Cpa6"] <- 78
smartnucleimet$cluster_id[smartnucleimet$cluster_label == "L6 CT Grp"] <- 79
smartnucleimet$cluster_id[smartnucleimet$cluster_label == "L6 CT Pou3f2"] <- 80
smartnucleimet$cluster_id[smartnucleimet$cluster_label == "L6 CT Kit_1"] <- 81
smartnucleimet$cluster_id[smartnucleimet$cluster_label == "L6 CT Kit_2"] <- 82
smartnucleimet$cluster_id[smartnucleimet$cluster_label == "Astro_Top2a"] <- 92
smartnucleimet$cluster_id[smartnucleimet$cluster_label == "Astro Aqp4_Gfap"] <- 93
smartnucleimet$cluster_id[smartnucleimet$cluster_label == "Astro Aqp4_Slc7a10"] <- 94


typenames <- unique(allmet$cluster_label)
typeids <- cbind.data.frame(typenames, unique(allmet$cluster_id), unique(allmet$cluster_color))
colnames(typeids) <- c("cluster_label","cluster_id","cluster_color")
typeids <- typeids[order(typeids$cluster_id),]
typenames <- typeids$cluster_label

subclassnames <- unique(allmet$subclass_label)
subids <- cbind.data.frame(subclassnames, unique(allmet$subclass_id),unique(allmet$subclass_color))
colnames(subids) <- c("subclass_label","subclass_id","subclass_color")
subids <- subids[order(subids$subclass_id),]
subnames <-subids$subclass_label

glutsub <- allmet %>% filter(class_label %in% "Glutamatergic")
glutsubnames <- unique(glutsub$subclass_label)
glutsubids <- cbind.data.frame(glutsubnames, unique(glutsub$subclass_id),unique(glutsub$subclass_color))
colnames(glutsubids) <- c("subclass_label","subclass_id","subclass_color")
glutsubids <- glutsubids[order(glutsubids$subclass_id),]
glutsubnames <- glutsubids$subclass_label

gabasub <- allmet %>% filter(class_label %in% "GABAergic")
gabasubnames <- unique(gabasub$subclass_label)
gabasubids <- cbind.data.frame(gabasubnames, unique(gabasub$subclass_id),unique(gabasub$subclass_color))
colnames(gabasubids) <- c("subclass_label","subclass_id","subclass_color")
gabasubids <- gabasubids[order(gabasubids$subclass_id),]
gabasubnames <- gabasubids$subclass_label

nonneuralsub <- allmet %>% filter(class_label %in% "Non-Neural")
nonneuralsubnames <- unique(nonneuralsub$subclass_label)
nonneuralsubids <- cbind.data.frame(nonneuralsubnames, unique(nonneuralsub$subclass_id),unique(nonneuralsub$subclass_color))
colnames(nonneuralsubids) <- c("subclass_label","subclass_id","subclass_color")
nonneuralsubids <- nonneuralsubids[order(nonneuralsubids$subclass_id),]
nonneuralsubnames <- nonneuralsubids$subclass_label

nonneuronalsub <- allmet %>% filter(class_label %in% "Non-Neuronal")
nonneuronalsubnames <- unique(nonneuronalsub$subclass_label)
nonneuronalsubids <- cbind.data.frame(nonneuronalsubnames, unique(nonneuronalsub$subclass_id),unique(nonneuronalsub$subclass_color))
colnames(nonneuronalsubids) <- c("subclass_label","subclass_id","subclass_color")
nonneuronalsubids <- nonneuronalsubids[order(nonneuronalsubids$subclass_id),]
nonneuronalsubnames <- nonneuronalsubids$subclass_label


classnames <- unique(allmet$class_label)
classids <- cbind.data.frame(classnames, unique(allmet$class_id),unique(allmet$class_color))
colnames(classids) <- c("class_label","class_id","class_color")
