library(tidyverse)
library(vroom)
library(data.table)
library(pheatmap)
library(preprocessCore)
library(purrr)
library(RColorBrewer)
library(ggpubr)

reverse_complement <- function(dna_seq) {
  complement <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
  nucleotides <- unlist(strsplit(dna_seq, ""))
  complement_nucleotides <- complement[nucleotides]
  reverse_complement_seq <- paste(rev(complement_nucleotides), collapse = "")
  return(reverse_complement_seq)
}

############# Look at Tau and PSI #################
# Calculate Tau for each gene.
calculate_tau <- function(row){
  # Remove NA values from the row
  non_na_row <- row[!is.na(row)]
  # If the row is empty after removing NAs, return NA
  if (length(non_na_row) < 20) {
    return(NA)
  }
  # Normalize the row by dividing by the max value of the non-NA row
  norm_row <- non_na_row / max(non_na_row)
  # Calculate tau using the number of non-NA values
  tau <- sum(1 - norm_row) / (length(non_na_row) - 1)
  return(tau)
}

# upsilon is the metric for PSI. We add 1 so that values close to 0 will not be inflated.
calculate_upsilon <- function(row) {
  # Remove NA values from the row
  non_na_row <- row[!is.na(row)]
  # If the row is empty after removing NAs, return NA
  if (length(non_na_row) < 20) {
    return(NA)
  }
  # Add 1 to every value in the row.
  non_na_row <- non_na_row + 1
  # Normalize the row by dividing by the max value of the non-NA row
  norm_row <- non_na_row / max(non_na_row)
  # Calculate tau using the number of non-NA values
  tau <- sum(1 - norm_row) / (length(non_na_row) - 1)
  return(tau * 2)
}

output_filepath <- "C:/Users/dawnxi/Dropbox (Harvard University)/02Splicing/SplicingManuscript/figure_outputs"
# Now I kinda want to plot the heatmap. 
final_psi_table_filtered <- fread("U:/processed_data/reprocess_250221/count_normalized_chimeric_rate_considering_included/WT_all_samples_PSI_count_table.csv")
final_psi_table_filtered <- final_psi_table_filtered %>% 
  filter(!(condition %in% c("K562WT", "K562K700E"))) %>% 
  filter(!(condition %in% c("JHOM1", "RVH421", "KNS60", "OVTOKO"))) %>% 
  mutate(total_count = included_count + skipped_count) %>%
  filter(total_count >= 20) %>%
  mutate(index_offset = paste(index, offset, sep = "__")) %>% 
  separate(offset, into = c("upstream_offset", "downstream_offset", "const_offset"), sep = ":") %>% 
  mutate(upstream_offset = as.integer(upstream_offset)) %>% 
  mutate(downstream_offset = as.integer(downstream_offset)) %>%
  mutate(const_offset = as.integer(const_offset)) %>% 
  filter(abs(upstream_offset) != 1 & abs(downstream_offset)!= 1) %>% 
  dplyr::select(-upstream_offset, -downstream_offset, -const_offset) %>% 
  dplyr::select(-index, -mode) %>%
  mutate(PSI = included_count/(included_count + skipped_count))



dt <- as.data.table(final_psi_table_filtered)

# Compute PSI by condition efficiently
psi_by_condition <- dt[, .(PSI = mean(PSI, na.rm = TRUE), num_rep = .N), by = .(condition, index_offset)][
  num_rep >= 2, .(condition, index_offset, PSI)]  # Filter out groups with <2 replicates

psi_table_pivot <- psi_by_condition %>%
  select(condition, index_offset, PSI) %>%
  pivot_wider(names_from = c(condition), values_from = PSI)

psi_table_pivot_sample <- final_psi_table_filtered %>%
  select(sample, index_offset, PSI) %>%
  pivot_wider(names_from = c(sample), values_from = PSI) 

# Convert to matrix.
psi_table_mat <- as.matrix(psi_table_pivot_sample %>% select(-index_offset))
rownames(psi_table_mat) <- psi_table_pivot_sample$index_offset

# Calculate upsilon for each row. 
upsilon_values <- psi_table_pivot %>% 
  select(-index_offset) %>% 
  apply(1, calculate_upsilon) %>% 
  as.data.frame() %>% 
  setNames("upsilon")

# Calculate reverse upsilon which is 1-values. 
psi_table_pivot_reverse <- psi_table_pivot %>% 
  select(-index_offset) %>% 
  mutate_all(~ 1 - .)  

# Calculate upsilon for each row.
upsilon_values_reverse <- psi_table_pivot_reverse %>% 
  apply(1, calculate_upsilon) %>% 
  as.data.frame() %>% 
  setNames("upsilon_reverse")


# Calculate tau for each row.
tau_values <- psi_table_pivot %>% 
  select(-index_offset) %>% 
  apply(1, calculate_tau) %>% 
  as.data.frame() %>% 
  setNames("tau")

# Calculate reverse tau which is 1-values.
tau_values_reverse <- psi_table_pivot_reverse %>% 
  apply(1, calculate_tau) %>% 
  as.data.frame() %>% 
  setNames("tau_reverse")

# Get num NA per row. 
num_na_per_row <- psi_table_pivot %>% 
  select(-index_offset) %>% 
  apply(1, function(x) sum(is.na(x))) %>% 
  as.data.frame() %>% 
  setNames("num_na")

# Get min of each row.
min_per_row <- psi_table_pivot %>% 
  select(-index_offset) %>% 
  apply(1, function(x) min(x, na.rm = TRUE)) %>% 
  as.data.frame() %>% 
  setNames("min")

# Get max of each row.
max_per_row <- psi_table_pivot %>% 
  select(-index_offset) %>% 
  apply(1, function(x) max(x, na.rm = TRUE)) %>% 
  as.data.frame() %>% 
  setNames("max")

# Get mean of each row. 
mean_per_row <- psi_table_pivot %>% 
  select(-index_offset) %>% 
  apply(1, function(x) mean(x, na.rm = TRUE)) %>% 
  as.data.frame() %>% 
  setNames("mean")

# Get the sample where the PSI is max.
max_sample_per_row <- psi_table_pivot %>% 
  select(-index_offset) %>% 
  apply(1, function(x) names(which.max(x))) %>% 
  as.data.frame() %>% 
  setNames("max_sample")

# Get the sample where the PSI is min.
min_sample_per_row <- psi_table_pivot %>% 
  select(-index_offset) %>% 
  apply(1, function(x) names(which.min(x))) %>% 
  as.data.frame() %>% 
  setNames("min_sample")

# Calculate matrix residuals of the PSI table.
psi_table_mat <- as.matrix(psi_table_pivot %>% select(-index_offset))
psi_table_mat_residuals <- psi_table_mat - rowMeans(psi_table_mat, na.rm = TRUE)
# Get residual sum across all rows. 
# Compute the number of non-NA values per row
sample_size <- rowSums(!is.na(psi_table_mat_residuals))
# Compute the residual sum and normalize by sample size
residual_sum <- rowSums(psi_table_mat_residuals^2, na.rm = TRUE) / sample_size
residual_sum <- as.data.frame(residual_sum) %>%
  setNames("residual_sum")

cell_specific_score <- data.frame(index_offset = psi_table_pivot$index_offset, 
                                  upsilon = upsilon_values$upsilon, 
                                  upsilon_reverse = upsilon_values_reverse$upsilon_reverse,
                                  tau = tau_values$tau,
                                  tau_reverse = tau_values_reverse$tau_reverse,
                                  num_na_per_row = num_na_per_row$num_na,
                                  row_max = max_per_row$max,
                                  row_min = min_per_row$min,
                                  mean_per_row = mean_per_row$mean,
                                  max_sample = max_sample_per_row$max_sample,
                                  min_sample = min_sample_per_row$min_sample,
                                  residual_sum = residual_sum$residual_sum) %>%
  mutate(max_min_diff = row_max - row_min) 

# Convert all numeric-like columns to numeric
cell_specific_score_numeric <- cell_specific_score %>%
  filter(num_na_per_row <= 10) %>%  
  select(-index_offset, - max_sample, - min_sample, -num_na_per_row) %>% 
  select(where(is.numeric))  # Keep only numeric columns

# Convert to long format for plotting
cell_specific_score_long <- cell_specific_score_numeric %>%

  pivot_longer(cols = everything(), names_to = "metric", values_to = "value")

##############################################################
##### General Histogram Plot #################################
##############################################################
histogram_plot <- ggplot(cell_specific_score_long, aes(x = value)) +
  geom_histogram(bins = 50, fill = "#1f78b4", color = "black", alpha = 0.8) +  # High contrast colors
  facet_wrap(~metric, scales = "free", ncol = 3) +  # Arrange in 3 columns for better layout
  labs(title = "Distribution of Metrics in Cell-Specific Score Table",
       x = "Value", y = "Frequency") +
  theme_classic(base_size = 16) +  # Larger base font size for readability
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18),  # Center title
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 14),
    strip.text = element_text(face = "bold", size = 14),  # Bold facet labels
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()
  )

ggsave(filename = paste0(output_filepath, "/all_general_metrics_histogram.pdf"), 
       plot = histogram_plot, 
       width = 12, height = 8, dpi = 300)


##############################################################
##### Plot Splice Site motifs ################################
##############################################################
library(ggseqlogo)
index_offset_table <- psi_table_pivot %>% 
  select(index_offset) %>% 
  separate(index_offset, into = c("ID", "offset"), sep = "__", remove = F) %>%
  separate(offset, into = c("upstream_offset", "downstream_offset", "const_offset"), sep = ":", remove = F) %>% 
  mutate(upstream_offset = as.integer(upstream_offset)) %>%
  mutate(downstream_offset = as.integer(downstream_offset)) %>%
  mutate(const_offset = as.integer(const_offset))


twist_barcodes <- read_csv("U:/melange/data/guide_library/20230130_twist_library_v3.csv") %>%
  mutate(barcodeRevcomp = sapply(barcode, reverse_complement))

index_offset_table <- index_offset_table %>% 
  left_join(twist_barcodes, by = c("ID" ))

# Adjust the sequences upstreamIntronSeq, skippedExonSeq, downstreamIntronSeq based on offset.
index_offset_table <- index_offset_table %>%
  mutate(upstreamIntron_len = nchar(upstreamIntronSeq)) %>%
  mutate(downstreamIntron_len = nchar(downstreamIntronSeq)) %>%
  mutate(skippedExon_len = nchar(skippedExonSeq)) %>%
  mutate(upstreamIntron_len_adj = upstreamIntron_len + upstream_offset) %>%
  mutate(downstreamIntron_len_adj = downstreamIntron_len - downstream_offset) %>%
  mutate(skippedExon_len_adj = skippedExon_len - upstream_offset + downstream_offset) %>%
  mutate(upstreamIntronSeq_adj = substr(librarySequence, 1, upstreamIntron_len_adj)) %>%
  mutate(skippedExonSeq_adj = substr(librarySequence, upstreamIntron_len_adj + 1, upstreamIntron_len_adj + skippedExon_len_adj)) %>%
  mutate(downstreamIntronSeq_adj = substr(librarySequence, upstreamIntron_len_adj + skippedExon_len_adj + 1, upstreamIntron_len_adj + skippedExon_len_adj + downstreamIntron_len_adj))


# Reference offset only.
ref_offset <- index_offset_table %>% 
  filter(offset == "0:0:0")
# Function to safely extract substrings
safe_substr <- function(seq, start, stop) {
  ifelse(nchar(seq) >= stop, substr(seq, start, stop), NA)
}

# Extract Upstream Splice Site Motif (Last 20 bp of Upstream Intron + First 10 bp of Skipped Exon)
ref_offset <- ref_offset %>%
  mutate(
    upstream_motif = paste0(
      safe_substr(upstreamIntronSeq_adj, nchar(upstreamIntronSeq_adj) - 19, nchar(upstreamIntronSeq_adj)),
      safe_substr(skippedExonSeq_adj, 1, 3)
    )
  )

# Extract Downstream Splice Site Motif (Last 10 bp of Skipped Exon + First 20 bp of Downstream Intron)
ref_offset <- ref_offset %>%
  mutate(
    downstream_motif = paste0(
      safe_substr(skippedExonSeq_adj, nchar(skippedExonSeq_adj) - 2, nchar(skippedExonSeq_adj)),
      safe_substr(downstreamIntronSeq_adj, 1, 6)
    )
  )

# Remove any NA values
ref_offset <- ref_offset %>%
  filter(!is.na(upstream_motif) & !is.na(downstream_motif))

# Plot Upstream Splice Site Motif Logo
p_upstream <- ggseqlogo(ref_offset$upstream_motif, method = "bits") 
p_downstream <- ggseqlogo(ref_offset$downstream_motif, method = "bits")

# Save the plots
p <- gridExtra::grid.arrange(p_upstream, p_downstream, ncol = 2)
ggsave(filename = paste0(output_filepath, "/splice_site_motif_ref_logo.pdf"), 
       plot = p, 
       width = 12, height = 3, dpi = 300)


# Plot for non-ref
non_ref_offset <- index_offset_table %>% 
  filter(offset != "0:0:0")

# Extract Upstream Splice Site Motif (Last 20 bp of Upstream Intron + First 10 bp of Skipped Exon)
non_ref_offset <- non_ref_offset %>%
  mutate(
    upstream_motif = paste0(
      safe_substr(upstreamIntronSeq_adj, nchar(upstreamIntronSeq_adj) - 19, nchar(upstreamIntronSeq_adj)),
      safe_substr(skippedExonSeq_adj, 1, 3)
    )
  )

# Extract Downstream Splice Site Motif (Last 10 bp of Skipped Exon + First 20 bp of Downstream Intron)
non_ref_offset <- non_ref_offset %>%
  mutate(
    downstream_motif = paste0(
      safe_substr(skippedExonSeq_adj, nchar(skippedExonSeq_adj) - 2, nchar(skippedExonSeq_adj)),
      safe_substr(downstreamIntronSeq_adj, 1, 6)
    )
  )

# Remove any NA values
non_ref_offset <- non_ref_offset %>%
  filter(!is.na(upstream_motif) & !is.na(downstream_motif)) %>% 
  filter(nchar(upstream_motif) == 23 & nchar(downstream_motif) == 9)

# Plot Upstream Splice Site Motif Logo
p_upstream <- ggseqlogo(non_ref_offset$upstream_motif, method = "bits") + ylim(0,2)
p_downstream <- ggseqlogo(non_ref_offset$downstream_motif, method = "bits") + ylim(0,2)
p <- gridExtra::grid.arrange(p_upstream, p_downstream, ncol = 2)
ggsave(filename = paste0(output_filepath, "/splice_site_motif_non_ref_logo.pdf"), 
       plot = p, 
       width = 12, height = 3, dpi = 300)


##############################################################
##### Plot Clustering by tissue type #########################
##############################################################
library(ggplot2)
library(ggrepel)
library(viridis)
library(umap)
library(dplyr)
library(readr)
library(tidyr)
library(FactoMineR)

# Load metadata
cellline_metadata <- read_csv("U:/melange/data/cellline_data_full_metadata.csv") %>% 
  select(StrippedName, Disease, `Disease Subtype`, lineage)

# Preprocess PSI matrix
psi_table_mat_clean <- psi_table_pivot %>%
  select(-index_offset) %>%
  mutate(across(everything(), as.numeric)) %>%  # Convert to numeric
  select(-HEK)

# Change the Kelly column to KELLY
colnames(psi_table_mat_clean) <- gsub("Kelly", "KELLY", colnames(psi_table_mat_clean))

# Remove rows that have > 10% NA
psi_table_mat_clean <- psi_table_mat_clean[rowMeans(is.na(psi_table_mat_clean)) <= 0.01, ]

# Replace NA values with row medians
psi_table_mat_clean <- psi_table_mat_clean %>%
  mutate(across(everything(), ~ replace_na(., median(., na.rm = TRUE))))

# Convert to matrix
psi_matrix <- as.matrix(psi_table_mat_clean)

# Ensure row names are properly assigned
rownames(psi_matrix) <- psi_table_pivot$index_offset[1:nrow(psi_matrix)]

# Remove rows with zero variance (needed for PCA & UMAP)
psi_matrix <- psi_matrix[apply(psi_matrix, 1, var, na.rm = TRUE) > 0, ]

# **Transpose the matrix** (samples = rows, features = columns)
psi_matrix_t <- t(psi_matrix)

### --- PCA ANALYSIS --- ###

# Perform PCA
pca_result <- prcomp(psi_matrix_t, center = TRUE, scale. = TRUE)

# Extract PCA scores
pca_scores <- as.data.frame(pca_result$x) %>%
  rownames_to_column(var = "StrippedName") %>%
  left_join(cellline_metadata, by = "StrippedName") %>% 
  distinct()

# Extract variance explained for axis labels
pca_var_explained <- summary(pca_result)$importance[2, ] * 100  # Percentage variance explained

# Plot PCA
pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Disease, label = StrippedName)) +
  geom_point(alpha = 0.8, size = 3) +  # Solid dots
  scale_color_brewer(palette = "Set1") +  
  geom_text_repel(aes(label = StrippedName), size = 3, box.padding = 0.5) +  # Add sample labels
  labs(
    title = "PCA of PSI Values",
    x = paste0("PC1 (", round(pca_var_explained[1], 1), "% Variance)"),
    y = paste0("PC2 (", round(pca_var_explained[2], 1), "% Variance)"),
    color = "Disease Type"
  ) +
  theme_classic(base_size = 16) +  # Classic theme without grid lines
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(size = 1),  # Keep x/y axis lines
    legend.position = "right",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(filename = paste0(output_filepath, "/pca_plot.pdf"), 
       plot = pca_plot, 
       width = 8, height = 6, dpi = 600)  # High DPI for publication quality

 ### --- UMAP ANALYSIS --- ###

# Perform UMAP
umap_result <- umap(psi_matrix_t)

# Extract UMAP scores
umap_scores <- as.data.frame(umap_result$layout) %>%
  rownames_to_column(var = "StrippedName") %>%
  rename(UMAP1 = V1, UMAP2 = V2) %>%
  left_join(cellline_metadata, by = "StrippedName") %>% 
  distinct()

# Plot UMAP
umap_plot <- ggplot(umap_scores, aes(x = UMAP1, y = UMAP2, color = Disease, label = StrippedName)) +
  geom_point(alpha = 0.8, size = 3) +  # Solid dots
  scale_color_brewer(palette = "Set1") +  
  geom_text_repel(aes(label = StrippedName), size = 3, box.padding = 0.5) +  # Add sample labels
  labs(
    title = "UMAP of PSI Values",
    x = "UMAP1",
    y = "UMAP2",
    color = "Disease Type"
  ) +
  theme_classic(base_size = 16) +  # Classic theme without grid lines
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(size = 1),  # Keep x/y axis lines
    legend.position = "right",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(filename = paste0(output_filepath, "/umap_plot.pdf"), 
       plot = umap_plot, 
       width = 8, height = 6, dpi = 600)  # High DPI for publication quality

# Plot UMAP
umap_plot <- ggplot(umap_scores, aes(x = UMAP1, y = UMAP2, color = `Disease Subtype`, label = StrippedName)) +
  geom_point(alpha = 0.8, size = 3) +  # Solid dots
  scale_color_brewer(palette = "Set1") +  
  geom_text_repel(aes(label = StrippedName), size = 3, box.padding = 0.5) +  # Add sample labels
  labs(
    title = "UMAP of PSI Values",
    x = "UMAP1",
    y = "UMAP2",
    color = "Disease Type"
  ) +
  theme_classic(base_size = 16) +  # Classic theme without grid lines
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(size = 1),  # Keep x/y axis lines
    legend.position = "right",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(filename = paste0(output_filepath, "/umap_plot_subtype.pdf"), 
       plot = umap_plot, 
       width = 8, height = 6, dpi = 600)  # High DPI for publication quality
 
# What are the sequences that contribute most to PC2? 
# Get the loadings for PC2
pca_loadings <- as.data.frame(pca_result$rotation) %>%
  rownames_to_column(var = "index_offset") %>%
  select(index_offset, PC2) %>%
  arrange(desc(abs(PC2))) %>% 
  head(20)
# Plot PC2 loadings
pca_loadings_plot <- ggplot(pca_loadings, aes(x = reorder(index_offset, PC2), y = PC2)) +
  geom_bar(stat = "identity", fill = "#1f78b4", color = "black") +  # High contrast colors
  coord_flip() +  # Horizontal bars
  labs(
    title = "Top Loadings for PC2",
    x = "Index Offset",
    y = "PC2 Loading"
  ) +
  theme_classic(base_size = 16) +  # Classic theme without grid lines
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(size = 1),  # Keep x/y axis lines
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
# pheatmap of these top loadings
top_loadings_mat <- as.matrix(psi_table_pivot_sample %>% 
                               filter(index_offset %in% pca_loadings$index_offset) %>% 
                               select(-index_offset))

# Plot pheatmap
pheatmap(
  top_loadings_mat, 
  cluster_rows = T, 
  cluster_cols = F, 
  color = viridis(100), 
  breaks = seq(0, 1, length.out = 100), 
  main = "Top Loadings for PC2",
  show_rownames = T,
  show_colnames = T 
)
