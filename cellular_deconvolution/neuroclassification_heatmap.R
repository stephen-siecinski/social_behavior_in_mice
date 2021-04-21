# Heatmap using neuronal proportion estimate

library(scales)
library(pheatmap)

proportions_mat <- proportions[,1:14]
rownames(proportions_mat) <- proportions_mat$sample_id
proportions_mat <- proportions_mat %>%
  dplyr::select(-sample_id)
proportions_mat <- t(as.matrix(proportions_mat))

heatmap_annotations <- dplyr::select(pheno_data, sex, strain, brain_region, sociability) %>%
  mutate(sex = recode(sex, "F" = "female", "M" = "male"))

annotation_colors = list(
  strain = c(B6 = "blue", C58J = "light green"),
  brain_region = c(A = "blue", H = "light green"),
  sociability = c(high = "blue", low = "light green"),
  sex = c(female = "blue", male = "light green"))

plot <- pheatmap(proportions_mat, cluster_rows = TRUE, show_rownames = TRUE,
                 cluster_cols = TRUE, annotation_col = heatmap_annotations, 
                 annotation_colors =  annotation_colors, color = viridis(200), 
                 border_color = NA, cutree_cols = 2)

ggsave(filename = paste0(outputs, "20210219_neuro_clustering_heatmap.png"),
       plot = plot, device = "png",
       dpi = 300, height = 8, width = 10, units = "in")
