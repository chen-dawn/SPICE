library(tidyverse)
library(vroom)
library(data.table)
library(Biostrings)
library(ggpointdensity)
library(pheatmap)

output_filepath <- "/Users/dawnxi/melange/figures_outputs/fig05/"
color_palette1 <- c("#E96052", "#F08944", "#F9AB57", "#FFD16D", "#FFE7B8", "#AADDE1", "#70BDD6", "#4F8DA3", "#336695", "#1A426D")
color_palette1_custom <- colorRampPalette(color_palette1)(100)
color_palette1_custom_rev <- colorRampPalette(rev(color_palette1))(100)

# color_palette2 <- c("#4575B4", "#80AED1", "#DFF2F6", "#FFE7B8", "#F08944", "#D83629")
color_palette2 <- c(
  "#4575B4",  # deep blue
  "#85B6D6",  # slightly lighter/more even blue
  "#E2EFF2",  # less stark pastel blue
  "#FFE3B0",  # warmer, slightly less saturated yellow
  "#EF9651",  # softer orange
  "#D83629"   # red
)
color_palette2_custom <- colorRampPalette(color_palette2)(100)
color_palette2_custom_rev <- colorRampPalette(rev(color_palette2))(100)

input_path <- "/Users/dawnxi/Harvard University Dropbox/Dawn Chen/02Splicing/for_kai/no_mutagenesis_data/"
# Do the sequence scatters first.
filepaths_barcode <- list.files(path = input_path, pattern = "predictions_barcode_seed.*\\.csv", full.names = TRUE)
filepaths_gene_barcode <- list.files(path = input_path, pattern = "psi_predictions_barcode_gene_seed.*\\.csv", full.names = TRUE)
filepaths_barcode_across_celltypes <- list.files(path = input_path, pattern = "predictions_barcode_across_cell_types_seed.*\\.csv", full.names = TRUE)
filepaths_gene_barcode_across_celltypes <- list.files(path = input_path, pattern = "psi_predictions_barcode_gene_across_cell_types_seed.*\\.csv", full.names = TRUE)

######################################################
# Plot the individual plots for the barcode-only files
######################################################
for (i in 1:length(filepaths_barcode)) {
  filepath <- filepaths_barcode[i]
  df <- read_csv(filepath)
  
  # Calculate correlations
  pearson_cor <- cor.test(df$Real_PSI, df$Predicted_PSI, method = "pearson")
  spearman_cor <- cor.test(df$Real_PSI, df$Predicted_PSI, method = "spearman")
  
  # Build plot
  p1 <- ggplot(df, aes(x = Real_PSI, y = Predicted_PSI)) + 
    stat_bin_2d(bins = 100) + 
    scale_fill_gradientn(colors = color_palette2_custom, trans = "log10", name = "log10(Count)") +
    theme_classic(base_size = 14) + 
    xlim(-12, 12) + ylim(-12, 12) +
    coord_fixed() +
    annotate(
      "text",
      x = -11.5, y = 11.5, hjust = 0, vjust = 1, size = 5,
      label = paste0(
        "Pearson r = ", formatC(pearson_cor$estimate, format = "f", digits = 3),
        "\nSpearman \u03C1 = ", formatC(spearman_cor$estimate, format = "f", digits = 3)
      )
    )
  
  ggsave(
    filename = paste0(output_filepath, "fig05_sequence_scatter_barcode_", i, ".pdf"),
    plot = p1,
    width = 8,
    height = 8,
    dpi = 300
  )
}

######################################################
# Plot the individual plots for the gene-barcode files
######################################################
for (i in 1:length(filepaths_gene_barcode)) {
  filepath <- filepaths_gene_barcode[i]
  df <- read_csv(filepath)
  
  # Calculate correlations
  pearson_cor <- cor.test(df$Real_PSI, df$Predicted_PSI, method = "pearson")
  spearman_cor <- cor.test(df$Real_PSI, df$Predicted_PSI, method = "spearman")
  
  # Build plot
  p1 <- ggplot(df, aes(x = Real_PSI, y = Predicted_PSI)) + 
    stat_bin_2d(bins = 100) + 
    scale_fill_gradientn(colors = color_palette2_custom, trans = "log10", name = "log10(Count)") +
    theme_classic(base_size = 14) + 
    xlim(-12, 12) + ylim(-12, 12) +
    coord_fixed() +
    annotate(
      "text",
      x = -11.5, y = 11.5, hjust = 0, vjust = 1, size = 5,
      label = paste0(
        "Pearson r = ", formatC(pearson_cor$estimate, format = "f", digits = 3),
        "\nSpearman \u03C1 = ", formatC(spearman_cor$estimate, format = "f", digits = 3)
      )
    )
  
  ggsave(
    filename = paste0(output_filepath, "fig05_sequence_scatter_gene_barcode_", i, ".pdf"),
    plot = p1,
    width = 8,
    height = 8,
    dpi = 300
  )
}



######################################################
# Plot the individual plots for the gene-barcode-across-cell-types files
######################################################
for (i in 1:length(filepaths_barcode_across_celltypes)) {
  filepath <- filepaths_barcode_across_celltypes[i]
  df <- read_csv(filepath)
  
  # Calculate correlations
  pearson_cor <- cor.test(df$Real_PSI, df$Predicted_PSI, method = "pearson")
  spearman_cor <- cor.test(df$Real_PSI, df$Predicted_PSI, method = "spearman")
  
  # Build plot
  p1 <- ggplot(df, aes(x = Real_PSI, y = Predicted_PSI)) + 
    stat_bin_2d(bins = 100) + 
    scale_fill_gradientn(colors = color_palette2_custom, trans = "log10", name = "log10(Count)") +
    theme_classic(base_size = 14) + 
    xlim(-12, 12) + ylim(-12, 12) +
    coord_fixed() +
    annotate(
      "text",
      x = -11.5, y = 11.5, hjust = 0, vjust = 1, size = 5,
      label = paste0(
        "Pearson r = ", formatC(pearson_cor$estimate, format = "f", digits = 3),
        "\nSpearman \u03C1 = ", formatC(spearman_cor$estimate, format = "f", digits = 3)
      )
    )
  
  ggsave(
    filename = paste0(output_filepath, "fig05_sequence_scatter_celltype_barcode_", i, ".pdf"),
    plot = p1,
    width = 8,
    height = 8,
    dpi = 300
  )
}

######################################################
# Plot the individual plots for the gene-barcode-across-cell-types files
######################################################
for (i in 1:length(filepaths_gene_barcode_across_celltypes)) {
  filepath <- filepaths_gene_barcode_across_celltypes[i]
  df <- read_csv(filepath)
  
  # Calculate correlations
  pearson_cor <- cor.test(df$Real_PSI, df$Predicted_PSI, method = "pearson")
  spearman_cor <- cor.test(df$Real_PSI, df$Predicted_PSI, method = "spearman")
  
  # Build plot
  p1 <- ggplot(df, aes(x = Real_PSI, y = Predicted_PSI)) + 
    stat_bin_2d(bins = 100) + 
    scale_fill_gradientn(colors = color_palette2_custom, trans = "log10", name = "log10(Count)") +
    theme_classic(base_size = 14) + 
    xlim(-12, 12) + ylim(-12, 12) +
    coord_fixed() +
    annotate(
      "text",
      x = -11.5, y = 11.5, hjust = 0, vjust = 1, size = 5,
      label = paste0(
        "Pearson r = ", formatC(pearson_cor$estimate, format = "f", digits = 3),
        "\nSpearman \u03C1 = ", formatC(spearman_cor$estimate, format = "f", digits = 3)
      )
    )
  
  ggsave(
    filename = paste0(output_filepath, "fig05_sequence_scatter_celltype_gene_barcode_", i, ".pdf"),
    plot = p1,
    width = 8,
    height = 8,
    dpi = 300
  )
}


##### Plot the Pearson Correlations #####
# Empty results data frame
cor_results <- tibble(
  dataset = character(),
  try = integer(),
  pearson_r = numeric()
)

######################################################
# Collect correlations for the barcode-only files
######################################################
for (i in seq_along(filepaths_barcode)) {
  filepath <- filepaths_barcode[i]
  df <- read_csv(filepath)
  
  pearson_cor <- cor.test(df$Real_PSI, df$Predicted_PSI, method = "pearson")
  
  cor_results <- bind_rows(
    cor_results,
    tibble(dataset = "barcode_only", try = i, pearson_r = pearson_cor$estimate)
  )
}

######################################################
# Collect correlations for the gene-barcode files
######################################################
for (i in seq_along(filepaths_gene_barcode)) {
  filepath <- filepaths_gene_barcode[i]
  df <- read_csv(filepath)
  
  pearson_cor <- cor.test(df$Real_PSI, df$Predicted_PSI, method = "pearson")
  
  cor_results <- bind_rows(
    cor_results,
    tibble(dataset = "gene_barcode", try = i, pearson_r = pearson_cor$estimate)
  )
}

######################################################
# Collect correlations for the barcode-across-celltypes files
######################################################
for (i in seq_along(filepaths_barcode_across_celltypes)) {
  filepath <- filepaths_barcode_across_celltypes[i]
  df <- read_csv(filepath)
  
  pearson_cor <- cor.test(df$Real_PSI, df$Predicted_PSI, method = "pearson")
  
  cor_results <- bind_rows(
    cor_results,
    tibble(dataset = "celltype_barcode", try = i, pearson_r = pearson_cor$estimate)
  )
}

######################################################
# Collect correlations for the gene-barcode-across-celltypes files
######################################################
for (i in seq_along(filepaths_gene_barcode_across_celltypes)) {
  filepath <- filepaths_gene_barcode_across_celltypes[i]
  df <- read_csv(filepath)
  
  pearson_cor <- cor.test(df$Real_PSI, df$Predicted_PSI, method = "pearson")
  
  cor_results <- bind_rows(
    cor_results,
    tibble(dataset = "celltype_gene_barcode", try = i, pearson_r = pearson_cor$estimate)
  )
}
######################################################
# Plot boxplot of Pearson r across tries (non-celltype)
######################################################
cor_results_noncell <- cor_results %>%
  filter(dataset %in% c("barcode_only", "gene_barcode"))

p_box_noncell <- ggplot(cor_results_noncell, aes(x = dataset, y = pearson_r, fill = dataset)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  theme_classic(base_size = 14) +
  labs(
    title = "Pearson correlation across 10 tries (Non Cell Type)",
    x = "Dataset",
    y = "Pearson r"
  ) +
  guides(fill = "none") +
  ylim(0.6, 0.8)

ggsave(
  filename = paste0(output_filepath, "fig05_pearson_boxplot_noncell.pdf"),
  plot = p_box_noncell,
  width = 6,
  height = 6,
  dpi = 300
)


######################################################
# Plot boxplot of Pearson r across tries (celltype)
######################################################
cor_results_cell <- cor_results %>%
  filter(dataset %in% c("celltype_barcode", "celltype_gene_barcode"))

p_box_cell <- ggplot(cor_results_cell, aes(x = dataset, y = pearson_r, fill = dataset)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  theme_classic(base_size = 14) +
  labs(
    title = "Pearson correlation across 10 tries (Cell Type)",
    x = "Dataset",
    y = "Pearson r"
  ) +
  guides(fill = "none") +
  ylim(0.8, 1)

ggsave(
  filename = paste0(output_filepath, "fig05_pearson_boxplot_cell.pdf"),
  plot = p_box_cell,
  width = 6,
  height = 6,
  dpi = 300
)

########## Plot heatmap for some individual barcode sequences #########
library(pheatmap)
set.seed(42)
df_new_barcode <- read_csv(filepaths_gene_barcode[1]) %>% 
  mutate(real_included = (1001*2^(Real_PSI) - 1)/(2^(Real_PSI) + 1),
         real_excluded = 1000 - real_included,
         predicted_included = (1001*2^(Predicted_PSI) - 1)/(2^(Predicted_PSI) + 1),
         predicted_excluded = 1000 - predicted_included) %>% 
  mutate(Real_PSI_bounded = (real_included)/(real_included + predicted_excluded)) %>% 
  mutate(Predicted_PSI_bounded = (predicted_included)/(predicted_included + real_excluded)) 


# Sample 20 random barcode
sampled_barcodes <- sample(unique(df_new_barcode$Barcode), 40)

# # plot heatmap of the observed vs predicted for each cell type
# for (barcode in sampled_barcodes) {
#   df_subset <- df_new_barcode %>%
#     filter(Barcode == barcode) %>%
#     select(Real_PSI_bounded,Predicted_PSI_bounded, Celltype) 
#   
#   if (nrow(df_subset) != 44) {
#     next  # Skip to the next iteration if no data for this barcode
#   }
# 
#   df_subset_mat <- as.matrix(df_subset %>% select(-Celltype))
#   rownames(df_subset_mat) <- df_subset$Celltype
#   df_subset_mat <- t(df_subset_mat)
#   print(ncol(df_subset_mat))
# 
#   p1 <- pheatmap(df_subset_mat, 
#                  color = color_palette2_custom, 
#                  cluster_rows = F,   # Allow clustering for better visualization
#                  cluster_cols = F,   # Allow clustering for better visualization
#                  fontsize = 8,         # Increased font size for readability
#                  border_color = NA,
#                  main = paste0("PSI Heatmap for Barcode: ", barcode),
#                  treeheight_row = 10,  # Reduce tree height for better spacing
#                  treeheight_col = 10,
#                  show_rownames = TRUE,
#                  angle_col = 90)       # Tilt column labels for readability  
#   ggsave(filename = paste0(output_filepath, "/fig05_PSI_heatmap_barcode_", barcode, ".pdf"),
#          plot = p1,
#          width = 8, height = 2, dpi = 300)
# }

library(gridExtra)

plot_list <- list()
breaks <- seq(0, 1, length.out = 101)

for (barcode in sampled_barcodes) {
  df_subset <- df_new_barcode %>%
    filter(Barcode == barcode) %>%
    select(Real_PSI_bounded, Predicted_PSI_bounded, Celltype) 
  
  if (nrow(df_subset) != 44) {
    next
  }
  
  df_subset_mat <- as.matrix(df_subset %>% select(-Celltype))
  rownames(df_subset_mat) <- df_subset$Celltype
  df_subset_mat <- t(df_subset_mat)
  
  p <- pheatmap(df_subset_mat, 
                color = color_palette2_custom,
                cluster_rows = FALSE,
                cluster_cols = FALSE,
                fontsize = 8,
                border_color = NA,
                main = paste0(barcode),
                treeheight_row = 10,
                treeheight_col = 10,
                show_rownames = TRUE,
                angle_col = 90,
                breaks = breaks)
  
  plot_list[[length(plot_list) + 1]] <- p[[4]]  # Grab the gtable
}

pdf(file = paste0(output_filepath, "fig05_PSI_heatmaps_barcode_stacked.pdf"), width = 8, height = 1.5 * length(plot_list))
grid.arrange(grobs = plot_list, ncol = 1)
dev.off()

########## Plot heatmap for some individual cell types #########
library(pheatmap)
set.seed(42)
df_new_barcode <- read_csv(filepaths_gene_barcode_across_celltypes[1]) %>% 
  mutate(real_included = (1001*2^(Real_PSI) - 1)/(2^(Real_PSI) + 1),
         real_excluded = 1000 - real_included,
         predicted_included = (1001*2^(Predicted_PSI) - 1)/(2^(Predicted_PSI) + 1),
         predicted_excluded = 1000 - predicted_included) %>% 
  mutate(Real_PSI_bounded = (real_included)/(real_included + predicted_excluded)) %>% 
  mutate(Predicted_PSI_bounded = (predicted_included)/(predicted_included + real_excluded)) 


# There are only a few cell types
breaks <- seq(0, 1, length.out = 101)

for (celltype in unique(df_new_barcode$Celltype)[1:5]) {
  df_subset <- df_new_barcode %>%
    filter(Celltype == celltype) %>%
    select(Real_PSI_bounded, Predicted_PSI_bounded, Barcode) 
  
  df_subset_mat <- as.matrix(df_subset %>% select(-Barcode))
  rownames(df_subset_mat) <- df_subset$Barcode
  df_subset_mat <- t(df_subset_mat)
  
  p <- pheatmap(df_subset_mat, 
                color = color_palette2_custom,
                cluster_rows = FALSE,
                cluster_cols = TRUE,
                fontsize = 8,
                border_color = NA,
                main = paste0(celltype),
                treeheight_row = 10,
                treeheight_col = 10,
                show_rownames = FALSE,
                show_colnames = FALSE,
                angle_col = 90,
                breaks = breaks)
  
  # Save each heatmap as an individual PDF
  pdf(file = paste0(output_filepath, "fig05_PSI_heatmap_", celltype, ".pdf"), 
      width = 8, height = 3)
  grid.arrange(grobs = list(p[[4]]))
  dev.off()
  
  # Save each heatmap as an individual PNG
  png(file = paste0(output_filepath, "fig05_PSI_heatmap_", celltype, ".png"),
      width = 2400, height = 800, res = 300)
  grid.arrange(grobs = list(p[[4]]))
  dev.off()
}