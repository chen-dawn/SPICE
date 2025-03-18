library(tidyverse)
library(vroom)
library(data.table)
library(pheatmap)
library(preprocessCore)
library(purrr)
library(RColorBrewer)
library(ggpubr)

output_filepath <- "C:/Users/dawnxi/Dropbox (Harvard University)/02Splicing/SplicingManuscript/figure_outputs"
# Now I kinda want to plot the heatmap. 
final_psi_table_filtered <- fread("U:/processed_data/reprocess_250221/count_normalized_chimeric_rate_considering_included/WT_all_samples_PSI_count_table.csv")
final_psi_table_filtered <- final_psi_table_filtered %>% 
  mutate(index_offset = paste(index, offset, sep = "__")) %>% 
  filter(offset == "0:0:0") %>% 
  select(-index, -offset, -mode) %>%
  mutate(total_count = included_count + skipped_count) %>%
  filter(total_count >= 10) %>%
  mutate(PSI = included_count/(included_count + skipped_count)) 


final_psi_table_pivot <- final_psi_table_filtered %>%
  select(sample, index_offset, PSI) %>%
  pivot_wider(names_from = c(sample), values_from = PSI) 

# Plot a violin plot per cell line. 
# Define color palette
palette <- RColorBrewer::brewer.pal(8, "Set2")

# Create the plot
pastel_color <- "#A6CEE3"  # Light blue pastel from ColorBrewer

# Create the plot
p1 <- ggplot(final_psi_table_filtered, aes(x = condition, y = PSI)) + 
  geom_violin(fill = pastel_color, color = "black", trim = TRUE, alpha = 0.7) +  # Single pastel color
  # geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "black") +  # Boxplot with white fill
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "black") +  # Highlight median
  theme_classic(base_size = 14) +  
  labs(
    title = "Distribution of PSI Across Conditions",
    x = "Condition",
    y = "Percent Spliced-In (PSI)"
  ) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

ggsave(filename = file.path(output_filepath, "all_samples_violin_plot_new_data_offset_zero.pdf"), plot = p1, width = 12, height = 6, dpi = 300)

# Plot cross replicate correlation for all samples
# Calculate correlation matrix
heatmap_colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)  # Better Blue-White-Red
breaks <- seq(0.5, 1, length.out = 100)

# Compute correlation matrix
correlation_matrix <- cor(final_psi_table_pivot %>% select(-index_offset), use = "pairwise.complete.obs")

# Plot heatmap with enhanced aesthetics
p2 <- pheatmap(correlation_matrix, 
               color = heatmap_colors, 
               cluster_rows = F,   # Allow clustering for better visualization
               cluster_cols = F,   # Allow clustering for better visualization
               fontsize = 8,         # Increased font size for readability
               border_color = "grey90", # Subtle grid lines
               breaks = breaks,
               main = "Cross-Replicate Correlation Heatmap",
               treeheight_row = 10,  # Reduce tree height for better spacing
               treeheight_col = 10,
               angle_col = 45)       # Tilt column labels for readability

# Save high-resolution image
ggsave(file.path(output_filepath,"CrossReplicate_Correlation_Heatmap_new_data_offset_zero.pdf"), plot = p2, dpi = 300, width = 17, height = 16)
  
# Filter out the bad conditions.
bad_conditions <- c("DAOY", "COLO783", "JHOM1", "RVH421", "OVTOKO", "K562WT", "K562K700E")

final_psi_table_filtered_filtered <- final_psi_table_filtered %>% 
  filter(!condition %in% bad_conditions)

final_psi_table_pivot_filtered <- final_psi_table_filtered_filtered %>%
  select(sample, index_offset, PSI) %>%
  pivot_wider(names_from = c(sample), values_from = PSI)

# Compute correlation matrix
correlation_matrix_filtered <- cor(final_psi_table_pivot_filtered %>% select(-index_offset), use = "pairwise.complete.obs")

# Plot heatmap with enhanced aesthetics
p3 <- pheatmap(correlation_matrix_filtered, 
               color = heatmap_colors, 
               cluster_rows = F,   # Allow clustering for better visualization
               cluster_cols = F,   # Allow clustering for better visualization
               fontsize = 8,         # Increased font size for readability
               border_color = NA, # Subtle grid lines
               breaks = breaks,
               main = "Cross-Replicate Correlation Heatmap (Filtered Conditions)",
               treeheight_row = 10,  # Reduce tree height for better spacing
               treeheight_col = 10,
               angle_col = 45)       # Tilt column labels for readability

# Save high-resolution image
ggsave(file.path(output_filepath,"CrossReplicate_Correlation_Heatmap_Filtered.pdf"), plot = p3, dpi = 300, width = 17, height = 16)

all_samples_df <- fread("C:/Users/dawnxi/Dropbox (Harvard University)/02Splicing/latest/all_sample_reps_PSI.csv")

all_samples_df <- all_samples_df  %>% filter((count + skipped) > 10) 

all_samples_pivot <- all_samples_df %>% 
  select(index, sample, PSI) %>%
  pivot_wider(names_from = c(sample), values_from = PSI)

# Plot the violin plot for all samples.
p4 <- ggplot(all_samples_df, aes(x = condition, y = PSI)) + 
  geom_violin(fill = pastel_color, color = "black", trim = TRUE, alpha = 0.7) +  # Single pastel color
  # geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "black") +  # Boxplot with white fill
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "black") +  # Highlight median
  theme_classic(base_size = 14) +  
  labs(
    title = "Distribution of PSI Across Conditions",
    x = "Condition",
    y = "Percent Spliced-In (PSI)"
  ) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

ggsave(filename = file.path(output_filepath, "all_samples_violin_plot_old_data.pdf"), plot = p4, width = 12, height = 6, dpi = 300)

# Compute correlation matrix
correlation_matrix_all_samples <- cor(all_samples_pivot %>% select(-index), use = "pairwise.complete.obs")

# Plot heatmap with enhanced aesthetics
p5 <- pheatmap(correlation_matrix_all_samples, 
               color = heatmap_colors, 
               cluster_rows = F,   # Allow clustering for better visualization
               cluster_cols = F,   # Allow clustering for better visualization
               fontsize = 8,         # Increased font size for readability
               border_color = "grey90", # Subtle grid lines
               main = "Cross-Replicate Correlation Heatmap (All Samples)",
               treeheight_row = 10,  # Reduce tree height for better spacing
               treeheight_col = 10,
               angle_col = 45)       # Tilt column labels for readability

# Save high-resolution image
ggsave(file.path(output_filepath,"CrossReplicate_Correlation_Heatmap_old_data.pdf"), plot = p5, dpi = 300, width = 17, height = 16)