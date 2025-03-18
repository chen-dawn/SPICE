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

cell_specific_score <- data.frame(index_offset = psi_table_pivot$index_offset, 
                                  upsilon = upsilon_values$upsilon, 
                                  upsilon_reverse = upsilon_values_reverse$upsilon_reverse,
                                  tau = tau_values$tau,
                                  tau_reverse = tau_values_reverse$tau_reverse,
                                  num_na_per_row = num_na_per_row$num_na,
                                  row_max = max_per_row$max,
                                  row_min = min_per_row$min,
                                  max_sample = max_sample_per_row$max_sample,
                                  min_sample = min_sample_per_row$min_sample) %>%
  filter(num_na_per_row < 10)

# ggplot histogram
ggplot(cell_specific_score, aes(x = upsilon)) + 
  geom_histogram(bins = 50) + 
  labs(title = "Distribution of Upsilon Values", x = "Upsilon", y = "Frequency") + 
  theme_classic(base_size = 14) + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggplot(cell_specific_score, aes(x = tau)) + 
  geom_histogram(bins = 50) + 
  labs(title = "Distribution of Tau Values", x = "Tau", y = "Frequency") + 
  theme_classic(base_size = 14) + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5))


# Shortlist to upsilon > 0.75. 
cell_specific_score_shortlist <- cell_specific_score %>% 
  filter(upsilon > 0.7 | upsilon_reverse > 0.7) %>% 
  # filter(upsilon >0.5 & upsilon < 0.51) %>% 
  arrange(desc(upsilon))


# shortlist psi_table_mat
psi_table_mat_high_upsilon <- psi_table_mat[rownames(psi_table_mat) %in% cell_specific_score_shortlist$index_offset, ]
# Heatmap
heatmap_colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)  # Better Blue-White-Red

p1 <- pheatmap(psi_table_mat_high_upsilon, 
         color = heatmap_colors, 
         cluster_rows = T,   # Allow clustering for better visualization
         cluster_cols = F,   # Allow clustering for better visualization
         fontsize = 8,         # Increased font size for readability
         border_color = "grey90", # Subtle grid lines
         main = "PSI Heatmap for High Upsilon Values",
         treeheight_row = 10,  # Reduce tree height for better spacing
         treeheight_col = 10,
         show_rownames = FALSE,
         angle_col = 45)       # Tilt column labels for readability
ggsave(filename = paste0(output_filepath, "/PSI_heatmap_high_upsilon.pdf"), 
       plot = p1, 
       width = 12, height = 8, dpi = 300)

## Get number of cell types in ecah group
num_cell_type_high_tau <- cell_specific_score_shortlist %>% 
  mutate(target_cell_type = ifelse(upsilon > 0.5, max_sample, min_sample)) %>% 
  group_by(target_cell_type) %>%
  summarise(num_cell_type = n()) 

theme_publication <- theme_classic(base_size = 9) + 
  theme(
    strip.text = element_text(face = "bold", size = 9),  # Bold facet labels
    axis.text = element_text(size = 9),  # Readable axis labels
    axis.title = element_text(size = 9, face = "bold"),  # Emphasized axis titles
    legend.position = "right",  # Show legend
    panel.spacing = unit(0.8, "lines")  # Reduce facet spacing
  )

ggplot(num_cell_type_high_tau, aes(x = target_cell_type, y = num_cell_type)) + 
  geom_bar(stat = "identity", fill = "#4CAF50") + 
  labs(title = "Number of Cell Types with High Upsilon Values", x = "Cell Type", y = "Count") + 
  theme_publication + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) # Add some space at the top
ggsave(filename = paste0(output_filepath, "/num_cell_type_high_tau.pdf"), 
       width = 12, height = 8, dpi = 300)

# Plot only the sequences high in each cell type.
for (i in 1:nrow(num_cell_type_high_tau)) {
  cell_type <- num_cell_type_high_tau$target_cell_type[i]
  num_cell_type <- num_cell_type_high_tau$num_cell_type[i]
  print(paste("Processing cell type:", cell_type, "with", num_cell_type, "sequences"))
  
  # Get the indices of the top num_cell_type sequences for this cell type
  top_indices <- cell_specific_score_shortlist %>% 
    mutate(target_cell_type = ifelse(upsilon > 0.5, max_sample, min_sample)) %>% 
    filter(target_cell_type == cell_type) %>% 
    arrange(desc(upsilon)) %>% 
    head(num_cell_type) %>% 
    pull(index_offset)
 
  # Create a heatmap for these sequences
  p2 <- pheatmap(psi_table_mat_high_upsilon[top_indices, ], 
           color = heatmap_colors, 
           cluster_rows = F,   # Allow clustering for better visualization
           cluster_cols = F,   # Allow clustering for better visualization
           fontsize = 8,         # Increased font size for readability
           border_color = "grey90", # Subtle grid lines
           main = paste("PSI Heatmap for High Upsilon Values in", cell_type),
           treeheight_row = 10,  # Reduce tree height for better spacing
           treeheight_col = 10,
           show_rownames = T,
           angle_col = 45)       # Tilt column labels for readability
  ggsave(filename = paste0(output_filepath, "/PSI_heatmap_high_upsilon_", cell_type, ".pdf"), p2,
         width = 20, height = 8, dpi = 300)
}

# ###### Look at sequences that are high in Kelly #####
# twist_barcodes <- read_csv("U:/melange/data/guide_library/20230130_twist_library_v3.csv") %>%
#   mutate(barcodeRevcomp = sapply(barcode, reverse_complement))
# 
# high_in_kelly_index <- cell_specific_score_shortlist %>% 
#   mutate(target_cell_type = ifelse(upsilon > 0.5, max_sample, min_sample)) %>% 
#   filter(target_cell_type == "Kelly") %>% 
#   # Split index offset by __
#   separate(index_offset, into = c("ID", "offset"), sep = "__", remove = F) %>% 
#   # Separate offset.
#   separate(offset, into = c("upstream_offset", "downstream_offset", "const_offset"), sep = ":") %>% 
#   mutate(upstream_offset = as.integer(upstream_offset)) %>% 
#   mutate(downstream_offset = as.integer(downstream_offset)) %>%
#   mutate(const_offset = as.integer(const_offset)) 
# 
# # > colnames(high_in_kelly_sequences)
# # [1] "ID"                  "offset"              "upsilon"            
# # [4] "upsilon_reverse"     "tau"                 "tau_reverse"        
# # [7] "num_na_per_row"      "row_max"             "row_min"            
# # [10] "max_sample"          "min_sample"          "target_cell_type"   
# # [13] "barcode"             "upstreamIntronSeq"   "skippedExonSeq"     
# # [16] "downstreamIntronSeq" "librarySequence"     "twistSequence"      
# # [19] "barcodeRevcomp"     
# # Pull the sequences from the twist barcodes
# high_in_kelly_sequences <- high_in_kelly_index %>% 
#   left_join(twist_barcodes, by = c("ID" )) 
# 
# # Adjust the sequences upstreamIntronSeq, skippedExonSeq, downstreamIntronSeq based on offset.
# high_in_kelly_sequences <- high_in_kelly_sequences %>% 
#   mutate(upstreamIntron_len = nchar(upstreamIntronSeq)) %>%
#   mutate(downstreamIntron_len = nchar(downstreamIntronSeq)) %>%
#   mutate(skippedExon_len = nchar(skippedExonSeq)) %>% 
#   mutate(upstreamIntron_len_adj = upstreamIntron_len + upstream_offset) %>% 
#   mutate(downstreamIntron_len_adj = downstreamIntron_len - downstream_offset) %>%
#   mutate(skippedExon_len_adj = skippedExon_len - upstream_offset + downstream_offset) %>%
#   mutate(upstreamIntronSeq_adj = substr(librarySequence, 1, upstreamIntron_len_adj)) %>%
#   mutate(skippedExonSeq_adj = substr(librarySequence, upstreamIntron_len_adj + 1, upstreamIntron_len_adj + skippedExon_len_adj)) %>%
#   mutate(downstreamIntronSeq_adj = substr(librarySequence, upstreamIntron_len_adj + skippedExon_len_adj + 1, upstreamIntron_len_adj + skippedExon_len_adj + downstreamIntron_len_adj))
#   
# 
# RBFOX_MOTIF <- "GCATG"
# # Look for binding sequence in the upstreamIntronSeq, skippedExonSeq, downstreamIntronSeq. Count how many per seq.
# high_in_kelly_sequences <- high_in_kelly_sequences %>%
#   mutate(
#     upstream_motif_count = str_count(upstreamIntronSeq_adj, RBFOX_MOTIF),
#     exon_motif_count = str_count(skippedExonSeq_adj, RBFOX_MOTIF),
#     downstream_motif_count = str_count(downstreamIntronSeq_adj, RBFOX_MOTIF),
#     total_motif_count = upstream_motif_count + exon_motif_count + downstream_motif_count
#   )
# 
# # Get the counts from final_psi_table_filtered.
# high_in_kelly_og <- final_psi_table_filtered %>% filter(condition == "Kelly") %>% 
#   filter(index_offset %in% high_in_kelly_sequences$index_offset) 
# 
# test <- final_psi_table_filtered %>% 
#   filter(grepl("ENSG00000113100.10;CDH9;chr5-27028253-27028370-26988105-26988382-27038462-27038583", index_offset))
