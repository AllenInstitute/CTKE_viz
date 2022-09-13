library(tidyverse)
library(Seurat)

outs_dir <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/Nik/Projects_for_other_people/Ray_Sanchez/M1_distance_mats/"
dir.create(outs_dir)


seurat_obj_exc <- readRDS("//allen/programs/celltypes/workgroups/rnaseqanalysis/Nik/Analyses_for_m1_paper/New_Human_mouse_marmoset_integration/Integration_12_09_19/Step2/Exc/sample.combined.cells")
seurat_obj_exc <- readRDS("//allen/programs/celltypes/workgroups/rnaseqanalysis/Nik/Analyses_for_m1_paper/New_Human_mouse_marmoset_integration/Integration_12_09_19/Step2/Inh/sample.combined.cells")
seurat_obj_exc <- readRDS("//allen/programs/celltypes/workgroups/rnaseqanalysis/Nik/Analyses_for_m1_paper/New_Human_mouse_marmoset_integration/Integration_12_09_19/Step2/Glia/sample.combined.cells")

combined_meta <- seurat_obj_exc@meta.data %>%
  as_tibble() %>%
  mutate(cluster = str_c(orig.ident, "xxx", cluster_label)) %>%
  select(sample_id, cluster)

reductions <- seurat_obj_exc@reductions

#load reductions and find species_cluster means
pc_means <- reductions$pca@cell.embeddings %>%
  as_tibble(rownames = "sample_id") %>%
  mutate(sample_id = sample_id %>% 
           str_remove_all("human_") %>%
           str_remove_all("marmoset_") %>%
           str_remove_all("mouse_")) %>%
  
  filter(sample_id %in% combined_meta$sample_id) %>%
  
  left_join(
    combined_meta,
    by = "sample_id"
  ) %>%
  
  select(-sample_id) %>%
  
  select(cluster, everything()) %>%
  
  group_by(cluster) %>%
  summarise_all(mean) %>%
  ungroup()

#calculate dist mat
dist_mat <- as.matrix(dist(pc_means[ , -1], method = "euclidean"))
colnames(dist_mat) <- rownames(dist_mat) <- pc_means$cluster

dist_tbl <- dist_mat %>% 
  
  as_tibble(rownames = "species_1") %>%
  
  gather(key = "species_2",
         value = "dist",
         -species_1) %>%
  
  separate(species_1, into = c("species_1", "cluster_1"), sep = "xxx") %>%
  separate(species_2, into = c("species_2", "cluster_2"), sep = "xxx") %>%
  
  #filter only human to other species comparisons
  filter(species_1 == "mouse",
         species_2 != "mouse") 



#plot euclid_dist heatmaps
dend <- readRDS("//allen/programs/celltypes/workgroups/rnaseqanalysis/Nik/Analyses_for_m1_paper/New_Human_mouse_marmoset_integration/Shiny_obj/Mouse_Broad_10xV3/dend.RData")
cluster_order <- labels(dend)
cluster_order <- cluster_order[cluster_order %in% dist_tbl$cluster_1]

#order factors and plot order
to_plot <- dist_tbl %>% 
  mutate(cluster_1 = cluster_1 %>% as_factor() %>% fct_relevel(cluster_order),
         species_2 = species_2 %>% as_factor() %>% fct_relevel(c("marmoset", "human"))) %>%
  arrange(species_2, cluster_2, dist) 

to_plot$cluster_1[to_plot$cluster_1 == "L5 IT Rspo1_1"] <- "L4/5 IT_1" 
to_plot$cluster_1[to_plot$cluster_1 == "L5 IT Rspo1_2"] <- "L4/5 IT_2" 
to_plot$cluster_1[to_plot$cluster_1 == "L5 IT S100b_1"] <- "L5 IT_1" 
to_plot$cluster_1[to_plot$cluster_1 == "L5 IT S100b1_2"] <- "L5 IT_2" 
to_plot$cluster_1[to_plot$cluster_1 == "L5 IT Pld5_1"] <- "L5 IT_3" 
to_plot$cluster_1[to_plot$cluster_1 == "L5 IT Pld5_2"] <- "L5 IT_4" 


species_2_plot_order <- to_plot %>%
  
  group_by(species_2, cluster_2) %>%
  summarise(cluster_2_nn = first(cluster_1)) %>%
  ungroup() %>%
  
  arrange(species_2, cluster_2_nn) %>%
  mutate(species_cluster_2 = str_c(species_2, "_", cluster_2))

#make plot


p1 <- to_plot %>%
  mutate(species_cluster_2 = str_c(species_2, "_", cluster_2) %>% 
           as_factor() %>% 
           fct_relevel(species_2_plot_order$species_cluster_2)) %>%
  
  ggplot() +
  geom_tile(aes(x = species_cluster_2, y = cluster_1, fill = log10(dist))) +
  
  labs(y = "human",
       x = "") +
  
  scale_fill_viridis_c(option = "B") +
  
  facet_grid(~species_2, scales = "free") +
  
  theme(aspect.ratio = 1,
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))

p1   

p2 <- to_plot %>%
  mutate(species_cluster_2 = str_c(species_2, "_", cluster_2) %>% 
           as_factor() %>% 
           fct_relevel(species_2_plot_order$species_cluster_2)) %>%
  
  ggplot() +
  geom_tile(aes(y = species_cluster_2, x = cluster_1, fill = log10(dist))) +
  
  labs(y = "", 
       x= "human") + 
  
  scale_fill_viridis_c(option = "B") +
  
  facet_wrap(~species_2, scales = "free_y", ncol = 1) +
  
  theme(aspect.ratio = 1,
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
        strip.background = element_blank(),
        strip.text = element_blank())
p2  

pdf(file = str_c(outs_dir, "glia_horizontal_heatmap.pdf"),
    useDingbats = FALSE,
    width = 10)
print(p1)
dev.off()

pdf(file = str_c(outs_dir, "glia_vertical_heatmap.pdf"),
    useDingbats = FALSE,
    height = 10)
print(p2)
dev.off()

write.csv(to_plot, "glia_mouse_distances.csv")
