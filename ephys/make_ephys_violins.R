#' Plotting some electrophysiology features from the Tolias lab Patch-seq data. Steps as follows
#'
#'  - 1. Load taxonomy metadata and patch-seq data
#'  - 2. Reconcile metadata from each dataset
#'  - 3. Write function to plot grid of violin plots comparing different ephys features
#'  - 4. Plot glutamatergic and GABAergic stuff separately
#'  - 5. Save plots as SVGs and data subset as CSVs
#'
#' @param taxonomy_file - standard annotated cluster table as referenced in many other scripts
#' @param patch_metadata - metadata file of all cells recorded from as part of Tolias lab Patch-seq experiments, needs cell ID
#' @param ephys_feats - data frame of ephys summary features obtained from Tolias lab Patch-seq experiments 
#' 

setwd("//allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/genexpr/taxonomy_and_genefiles")
clupdate <- read_csv("cl_updated.df.csv")

dat_dir <- "//allen/programs/celltypes/workgroups/humancelltypes/RayS/CTCard_Content/ephys/data"
setwd(dat_dir)
patch_metadata <- read.csv("Tolias_m1_patchseq_meta_data.csv",sep=",")
ephys_feats <- read.csv("Tolias_m1_patchseq_ephys_features.csv",sep=",")

make_ephys_plots <-  function(taxonomy_file,
                              patch_metadata,
                              ephys_feats){
  
  clupdate < taxonomy_file
  neuronmet <- clupdate %>% filter(class_label %in% "Glutamatergic" | class_label %in% "GABAergic")
  
  
  # Subset desired cell set and features
  glut_col <- clupdate %>%
    filter(class_label == "Glutamatergic")
  glut_col <- data.frame(glut_col,stringsAsFactors = FALSE)
  glut_col$subclass_label <- as.character(glut_col$subclass_label)
  glut_col$subclass_color <- as.character(glut_col$subclass_color)
  
  
  toMatch <- c("ET","IT","CT","NP")
  
  sub_metadata <- patch_metadata %>%
    filter(grepl(paste(toMatch,collapse="|"), patch_metadata$RNA.family)==TRUE) %>%
    rename(cell.id = Cell)
  
  
  
  merged_dat <- merge(sub_metadata,ephys_feats,by.y = "cell.id")
  
  merged_dat <- merged_dat %>%
    select(cell.id, RNA.family, RNA.type, AP.amplitude..mV., AP.width..ms.,
           Input.resistance..MOhm., Max.number.of.APs,
           Membrane.time.constant..ms., Upstroke.to.downstroke.ratio)
  
  # Match to mouse taxonomy labels/colors
  merged_dat$cluster_color <-  glut_col$cluster_color[match(merged_dat$RNA.type, glut_col$cluster_label)]
  merged_dat$cluster_label <-  glut_col$cluster_label[match(merged_dat$RNA.type, glut_col$cluster_label)]
  merged_dat$cluster_id <-  glut_col$cluster_id[match(merged_dat$RNA.type, glut_col$cluster_label)]
  merged_dat$subclass_label <-  glut_col$subclass_label[match(merged_dat$RNA.type, glut_col$cluster_label)]
  merged_dat$subclass_color <-  glut_col$subclass_color[match(merged_dat$RNA.type, glut_col$cluster_label)]
  merged_dat$subclass_id <-  glut_col$subclass_id[match(merged_dat$RNA.type, glut_col$cluster_label)]
  merged_dat$class_label <- rep("Glutamatergic",length(merged_dat$cell.id))
  merged_dat <- data.frame(merged_dat, stringsAsFactors = T)
  
  # Cleaned up labels with desired features - manual (although one could manually
  # change the variable names in the dataset)
  colnames(merged_dat)[4:9] <- c("AP Amplitude (mv)", 
                                 "AP Width (ms)", 
                                 "Input Resistance (MOhms)",
                                 "Max # of APs", 
                                 "Membrane Time Constant (ms)",
                                 "Upstroke:Downstroke ratio")
  
  # Function to plot all features
  plot_ephys_feats <- function(data,
                               X_plot,
                               level_id,
                               filter_level,
                               filter_by,
                               outname,
                               nodename,
                               setnames,
                               setids,
                               setcols) {
    
    X_plot <- enquo(X_plot)
    filter_level <- enquo(filter_level)
    level_id <- enquo(level_id)
    
    data <- data %>%
      filter(!!filter_level %in% filter_by)
    
    data <- data %>%
      mutate(X_plot = as_factor(!!X_plot)) 
    
    data <- data %>%
      arrange(!!level_id)
    
    data <- data %>%
      mutate(level_id = as.numeric(!!level_id)) %>%
      mutate(level_id = cut(!!level_id, breaks = length(setnames)))
    
    levels(data$level_id) <- setids
    catorder <- setnames
    data$X_plot <- factor(data$X_plot, catorder)
    
    data <- melt(data)
    
    data <- data[data$variable == "AP Amplitude (mv)" | data$variable == "AP Width (ms)" | data$variable == "Input Resistance (MOhms)" | 
                   data$variable =="Max # of APs" | data$variable == "Membrane Time Constant (ms)" | data$variable == "Upstroke:Downstroke ratio",]
    data <- data %>%
      filter(!!filter_level %in% filter_by)
    
    data <- data %>%
      mutate(X_plot = as_factor(!!X_plot)) 
    
    data <- data %>%
      arrange(level_id)
    
    
    levels(data$level_id) <- setids
    
    data$X_plot <- factor(data$X_plot, catorder)
    lineloc <- which(levels(as.factor(data$X_plot)) == nodename)
    
    p <- ggplot(data, aes(X_plot, value, fill = X_plot)) + 
      geom_rect(data = data,
                xmin = lineloc-0.5,
                xmax = lineloc+0.5,
                ymin = 0,
                ymax = Inf,
                fill = "grey91",
                alpha=0.11
      )+
      geom_violin(alpha=0.9,scale = "width", width = 1) + 
      scale_fill_manual(values = setcols, drop = FALSE) + 
      scale_x_discrete(drop = FALSE) +
      scale_color_manual(values = setcols) + 
      theme_light() +
      geom_jitter(fill = "black", width = 0.2, size=0.5) +
      theme(text = element_text(size=15),
            aspect.ratio = 0.29, legend.position = "none",
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.title.y.left = element_blank(),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, vjust = 0.85, hjust=1),
            axis.title.x = element_blank(),
            strip.text.y = element_text(
              size = 9, color = "black", face = "bold"
            ),
            panel.spacing = unit(0.85, "lines"),
            strip.background = element_rect(fill = "white", colour = "white", size = 0.8),
            strip.text = element_text(colour = "grey30",face = "bold")) 
    
    
    
    p <- p + facet_wrap(~variable, ncol=2,scales="free",shrink = TRUE,
                        labeller = labeller(value = c("AP Amplitude (mV)", "AP width (ms)",
                                                      "Input Resistance (MOhms)","Max # of APs",
                                                      "Membrane Time Constant (ms)", "Upstroke:downstroke")))
    p
    
    nodename <- gsub("/","-",nodename)
    svg(paste0(outname,nodename,"_characteristics.svg"),width=9.3, height=9.8)
    print(p)
    dev.off()
  }
  
  
  # Glutamatergic
  ######

  
  glutsubnames <- as.character(unique(glut_col$subclass_label))
  glutsubids <- unique(glut_col$subclass_id)
  glutsubcols <- unique(glut_col$subclass_color)
  
  for (i in 1:length(glutsubnames)) {
    
    plot_ephys_feats(
      data = merged_dat, 
      X_plot = subclass_label, 
      level_id = subclass_id,
      filter_level = class_label, 
      filter_by = "Glutamatergic", 
      outname = "electrophysiology_",
      nodename = glutsubnames[i],
      setnames = glutsubnames,
      setids = glutsubids,
      setcols = glutsubcols 
    )
    
  }
  
  gluttypenames <- as.character(unique(glut_col$cluster_label))
  gluttypeids <- unique(glut_col$cluster_id)
  gluttypecols <- unique(glut_col$cluster_color)
  
  missingsubs <- "L6 IT Car3"
  missingnodes <- setdiff(as.character(gluttypenames), unique(as.character(merged_dat$cluster_label)))
  
  patchglutsub <- glutsubnames[c(1:4,6:9)]
  
  
  
  for (i in 1:length(patchglutsub)) {
    temp <- merged_dat %>% filter(subclass_label %in% patchglutsub[i])
    temp <- temp[order(temp$cluster_id),]
    checktypes <- setdiff(clupdate$cluster_label[clupdate$subclass_label == patchglutsub[i]], unique(as.character(temp$cluster_label)))
    subsize <- length(unique(temp$cluster_label[temp$subclass_label == patchglutsub[i]]))
    
    if(length(checktypes) == 0 && subsize > 1) {
      temptypes <- as.character(unique(temp$cluster_label))
      for (j in 1:length(temptypes)) {
        ttnames <- temptypes
        ttids <- clupdate %>% filter(cluster_label %in% temptypes)
        plot_ephys_feats(data = merged_dat, 
                         filter_level = subclass_label, 
                         filter_by = patchglutsub[i],
                         X_plot = cluster_label,
                         level_id = cluster_id,
                         outname = "electrophysiology_",
                         nodename = ttnames[j],
                         setnames = ttnames,
                         setids = ttids$cluster_id,
                         setcols = ttids$cluster_color)
        
      }
    }
    
    if (subsize == 1 & length(checktypes) == 0) {
      temptypes <- as.character(unique(temp$cluster_label))
      whichparent <- unique(clupdate$class_label[clupdate$subclass_label == patchglutsub[i]])
      whichsubs <- unique(clupdate$subclass_label[clupdate$class_label == whichparent])
      whichids <- unique(clupdate$subclass_id[clupdate$class_label == whichparent])
      whichcols <- unique(clupdate$subclass_color[clupdate$class_label == whichparent])
      plot_ephys_feats(data = merged_dat, 
                       filter_level = class_label, 
                       filter_by = whichparent,
                       X_plot = subclass_label, 
                       level_id = subclass_id,
                       outname = "electrophysiology_TYPE_",
                       nodename = patchglutsub[i],
                       setnames = whichsubs,
                       setids = whichids,
                       setcols = whichcols)
    }
    
    
    
    if (length(checktypes > 0)) {
      temptypes <- c(as.character(unique(temp$cluster_label)), checktypes)
      whichids <- clupdate$cluster_id[clupdate$cluster_label %in% temptypes]
      whichtypes <- clupdate$cluster_label[clupdate$cluster_label %in% temptypes]
      whichcols <- clupdate$cluster_color[clupdate$cluster_label %in% temptypes]
      
      for (j in 1:length(temptypes)) {
        ttnames <- whichtypes
        ttids <- clupdate %>% filter(cluster_label %in% temptypes)
        plot_ephys_feats(data = merged_dat, 
                         filter_level = subclass_label, 
                         filter_by = patchglutsub[i],
                         X_plot = cluster_label, 
                         level_id = cluster_id,
                         outname = "electrophysiology_",
                         nodename = ttnames[j],
                         setnames = ttnames,
                         setids = whichids,
                         setcols = whichcols)
        
      }
    }
    
    
  }
  
  
  
  missingglut <- c(missingsubs,missingnodes)
  
  

  
  
  #######
  gaba_col <- clupdate %>%
    filter(class_label == "GABAergic")
  gaba_col <- data.frame(gaba_col,stringsAsFactors = FALSE)
  gaba_col$subclass_label <- as.character(gaba_col$subclass_label)
  gaba_col$subclass_color <- as.character(gaba_col$subclass_color)
  
  toMatch <- c("Lamp5","Pvalb","Sncg","Sst","Vip")
  
  sub_metadata <- patch_metadata %>%
    filter(grepl(paste(toMatch,collapse="|"), patch_metadata$RNA.family)==TRUE) %>%
    rename(cell.id = Cell)
  
  merged_dat <- merge(sub_metadata,ephys_feats,by.y = "cell.id")
  
  merged_dat <- merged_dat %>%
    select(cell.id, RNA.family, RNA.type, AP.amplitude..mV., AP.width..ms.,
           Input.resistance..MOhm., Max.number.of.APs,
           Membrane.time.constant..ms., Upstroke.to.downstroke.ratio)
  
  # Match to mouse taxonomy labels/colors
  merged_dat$cluster_color <-  gaba_col$cluster_color[match(merged_dat$RNA.type, gaba_col$cluster_label)]
  merged_dat$cluster_label <-  gaba_col$cluster_label[match(merged_dat$RNA.type, gaba_col$cluster_label)]
  merged_dat$cluster_id <-  gaba_col$cluster_id[match(merged_dat$RNA.type, gaba_col$cluster_label)]
  merged_dat$subclass_label <-  gaba_col$subclass_label[match(merged_dat$RNA.type, gaba_col$cluster_label)]
  merged_dat$subclass_color <-  gaba_col$subclass_color[match(merged_dat$RNA.type, gaba_col$cluster_label)]
  merged_dat$subclass_id <-  gaba_col$subclass_id[match(merged_dat$RNA.type, gaba_col$cluster_label)]
  merged_dat$class_label <- rep("GABAergic",length(merged_dat$cell.id))
  merged_dat <- data.frame(merged_dat, stringsAsFactors = T)
  
  
  colnames(merged_dat)[4:9] <- c("AP Amplitude (mv)", 
                                 "AP Width (ms)", 
                                 "Input Resistance (MOhms)",
                                 "Max # of APs", 
                                 "Membrane Time Constant (ms)",
                                 "Upstroke:Downstroke ratio")
  
  # GABAergic
  ######

  gabasubnames <- as.character(unique(gaba_col$subclass_label))
  gabasubids <- unique(gaba_col$subclass_id)
  gabasubcols <- unique(gaba_col$subclass_color)
  
  for (i in 1:length(gabasubnames)) {
    
    plot_ephys_feats(
      data = merged_dat, 
      X_plot = subclass_label, 
      level_id = subclass_id,
      filter_level = class_label, 
      filter_by = "GABAergic", 
      outname = "electrophysiology_",
      nodename = gabasubnames[i],
      setnames = gabasubnames,
      setids = gabasubids,
      setcols = gabasubcols 
    )
    
  }
  
  gabatypenames <- as.character(unique(gaba_col$cluster_label))
  gabatypeids <- unique(gaba_col$cluster_id)
  gabatypecols <- unique(gaba_col$cluster_color)
  
  missingsubs <- "Meis2"
  missingnodes <- setdiff(as.character(gabatypenames), unique(as.character(merged_dat$cluster_label)))
  
  patchgabasub <- gabasubnames
  
  
  for (i in 1:length(patchgabasub)) {
    temp <- merged_dat %>% filter(subclass_label %in% patchgabasub[i])
    temp <- temp[order(temp$cluster_id),]
    checktypes <- setdiff(clupdate$cluster_label[clupdate$subclass_label == patchgabasub[i]], unique(as.character(temp$cluster_label)))
    subsize <- length(unique(temp$cluster_label[temp$subclass_label == patchgabasub[i]]))
    
    if(length(checktypes) == 0 && subsize > 1) {
      temptypes <- as.character(unique(temp$cluster_label))
      for (j in 1:length(temptypes)) {
        ttnames <- temptypes
        ttids <- clupdate %>% filter(cluster_label %in% temptypes)
        plot_ephys_feats(data = merged_dat, 
                         filter_level = subclass_label, 
                         filter_by = patchgabasub[i],
                         X_plot = cluster_label,
                         level_id = cluster_id,
                         outname = "electrophysiology_",
                         nodename = ttnames[j],
                         setnames = ttnames,
                         setids = ttids$cluster_id,
                         setcols = ttids$cluster_color)
        
      }
    }
    
    if (subsize == 1 & length(checktypes) == 0) {
      temptypes <- as.character(unique(temp$cluster_label))
      whichparent <- unique(clupdate$class_label[clupdate$subclass_label == patchgabasub[i]])
      whichsubs <- unique(clupdate$subclass_label[clupdate$class_label == whichparent])
      whichids <- unique(clupdate$subclass_id[clupdate$class_label == whichparent])
      whichcols <- unique(clupdate$subclass_color[clupdate$class_label == whichparent])
      plot_ephys_feats(data = merged_dat, 
                       filter_level = class_label, 
                       filter_by = whichparent,
                       X_plot = subclass_label, 
                       level_id = subclass_id,
                       outname = "electrophysiology_TYPE_",
                       nodename = patchgabasub[i],
                       setnames = whichsubs,
                       setids = whichids,
                       setcols = whichcols)
    }
    
    
    
    if (length(checktypes > 0)) {
      temptypes <- c(as.character(unique(temp$cluster_label)), checktypes)
      whichids <- clupdate$cluster_id[clupdate$cluster_label %in% temptypes]
      whichtypes <- clupdate$cluster_label[clupdate$cluster_label %in% temptypes]
      whichcols <- clupdate$cluster_color[clupdate$cluster_label %in% temptypes]
      
      for (j in 1:length(temptypes)) {
        ttnames <- whichtypes
        ttids <- clupdate %>% filter(cluster_label %in% temptypes)
        plot_ephys_feats(data = merged_dat, 
                         filter_level = subclass_label, 
                         filter_by = patchgabasub[i],
                         X_plot = cluster_label, 
                         level_id = cluster_id,
                         outname = "electrophysiology_",
                         nodename = ttnames[j],
                         setnames = ttnames,
                         setids = whichids,
                         setcols = whichcols)
        
      }
    }
    
    
  }
  
  
}

make_ephys_plots(taxonomy_file = clupdate,
                 patch_metadata = patch_metadata, 
                 ephys_feats = ephys_feats)
