library(tidyverse)
library(vroom)
library(data.table)
library(pheatmap)
library(preprocessCore)
library(purrr)
library(RColorBrewer)
library(ggpubr)
library(viridis)

reverse_complement <- function(dna_seq) {
  complement <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
  nucleotides <- unlist(strsplit(dna_seq, ""))
  complement_nucleotides <- complement[nucleotides]
  reverse_complement_seq <- paste(rev(complement_nucleotides), collapse = "")
  return(reverse_complement_seq)
}

output_filepath <- "C:/Users/dawnxi/Dropbox (Harvard University)/02Splicing/SplicingManuscript/figure_outputs"

all_files_df_path <- "U:/processed_data/reprocess_250221/count_normalized_chimeric_rate_considering_included/WT_all_samples_raw_counts.csv"
all_files_df <- fread(all_files_df_path)

# > head(all_files_df)
# index
# <char>
#   1: ENSG00000000003.15;TSPAN6;chrX-100632484-100632568-100630758-100630866-100633404-100633539
# 2: ENSG00000000003.15;TSPAN6;chrX-100632484-100632568-100630758-100630866-100633404-100633539
# 3: ENSG00000000003.15;TSPAN6;chrX-100632484-100632568-100630758-100630866-100633404-100633539
# 4: ENSG00000000003.15;TSPAN6;chrX-100632484-100632568-100630758-100630866-100633404-100633539
# 5: ENSG00000000003.15;TSPAN6;chrX-100633930-100634029-100632484-100632568-100635177-100635252
# 6: ENSG00000000003.15;TSPAN6;chrX-100633930-100634029-100632484-100632568-100635177-100635252
# mode offset count count_scaled    sample condition
# <char> <char> <int>        <int>    <char>    <char>
#   1:  INCLUDED 0:-1:0     1            1 769P-rep1      769P
# 2:  INCLUDED  0:0:0    29           29 769P-rep1      769P
# 3:   SKIPPED      0     3            0 769P-rep1      769P
# 4: UNSPLICED    113     3            0 769P-rep1      769P
# 5:  INCLUDED  0:0:0    88           88 769P-rep1      769P
# 6:   SKIPPED      0    15            1 769P-rep1      769P

all_files_df_filtered <- all_files_df %>% 
  filter(!(condition %in% c("K562WT", "K562K700E"))) %>% 
  filter(!(condition %in% c("JHOM1", "RVH421", "KNS60", "OVTOKO"))) %>% 
  group_by(sample, index) %>% 
  mutate(total_count = sum(count_scaled)) %>% 
  filter(total_count >= 20) 

# Compute unspliced percentage
included_only <- all_files_df_filtered %>%
  filter(mode == "INCLUDED") %>%
  group_by(sample, index) %>% 
  mutate(total_included  = sum(count_scaled)) %>%
  filter(total_included >= 30) %>% 
  ungroup() %>%
  mutate(included_perc = count_scaled / total_included) 

included_only <- included_only %>% 
  separate(offset, into = c("upstream_offset", "downstream_offset", "const_offset"), sep = ":", remove = F) %>% 
  mutate(upstream_offset = as.integer(upstream_offset),
         downstream_offset = as.integer(downstream_offset)) %>%
  filter(abs(upstream_offset) != 1) %>% 
  filter(abs(downstream_offset) != 1) 

# Get how many included reads are there per index per sample. 
num_included_events_per_sample <- included_only %>% 
  group_by(sample, index) %>% 
  summarise(num_included_events = n_distinct(offset))

num_included_events_per_sample_wide <- num_included_events_per_sample %>% 
  pivot_wider(names_from = sample, values_from = num_included_events)

included_by_condition <- included_only %>% 
  group_by(condition, index, offset, upstream_offset, downstream_offset, const_offset) %>% 
  summarise(included_perc = mean(included_perc, na.rm = T), n= n()) %>%
  filter(n >= 2)


# Get number of splice sites per index
num_splice_sites_per_index <- included_by_condition %>% 
  group_by(index, offset) %>%
  mutate(num_occurances = n()) %>%
  filter(num_occurances >= 21) %>%
  ungroup() %>% 
  group_by(index) %>%
  summarise(num_splice_sites = n_distinct(offset))

p1 <- ggplot(num_splice_sites_per_index, aes(x = num_splice_sites)) +
  geom_histogram(binwidth = 1, fill = "gray70", color = "black") +  # High contrast bars
  scale_x_continuous(breaks = seq(min(num_splice_sites_per_index$num_splice_sites), 
                                  max(num_splice_sites_per_index$num_splice_sites), by = 1)) +  # Ensure integer ticks
  labs(title = "Distribution of Number of Splice Sites per Index",
       x = "Number of Alternative 3'ss",
       y = "Frequency") +
  theme_classic(base_size = 16) +  # Remove grid lines, improve readability
  theme(
    axis.text.x = element_text(size = 14),  # Readable x-axis
    axis.text.y = element_text(size = 14),  # Readable y-axis
    axis.title = element_text(size = 16),   # Readable axis titles
    plot.title = element_text(size = 18, hjust = 0.5),  # Centered title
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )

# Save as a high-resolution PDF
ggsave(file.path(output_filepath, "splice_sites_distribution.pdf"), 
       plot = p1, width = 10, height = 6, dpi = 300)  # Publication quality

# Select top 2 sequences with most alternative 3'ss
top_5_sequences <- num_splice_sites_per_index %>%
  top_n(2, num_splice_sites) %>%
  pull(index)

# Prepare data
included_only_filtered <- included_by_condition %>%
  ungroup() %>% 
  filter(index %in% top_5_sequences) %>%  
  mutate(index_offset = paste(index, offset, sep = "__")) %>%
  select(index_offset, condition, included_perc) %>% 
  pivot_wider(names_from = condition, values_from = included_perc) 

# Order rows alphabetically by index_offset
included_only_filtered <- included_only_filtered %>%
  arrange(index_offset)

# Filter rows that have < 50% missing values
included_only_filtered <- included_only_filtered %>%
  filter(rowMeans(is.na(included_only_filtered[, -1])) < 0.5)

# Convert to matrix
mat <- as.matrix(included_only_filtered[, -1])
rownames(mat) <- included_only_filtered$index_offset

# Identify row indices for adding a gap
index_labels <- included_only_filtered$index_offset
sequence_labels <- sapply(strsplit(index_labels, "__"), `[`, 1)  # Extract sequence ID
sequence_change_points <- which(diff(as.numeric(factor(sequence_labels))) != 0)  # Find index where sequence changes

# Save heatmap as high-resolution PDF
pdf(file.path(output_filepath, "3ss_top_2_sequences_heatmap.pdf"), width = 24, height = 8)  # Adjust dimensions for readability
pheatmap(mat,
         cluster_rows = FALSE,  # Maintain alphabetical ordering
         cluster_cols = F,   # Cluster samples (columns)
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = viridis(100),
         main = "Top 2 Sequences with Most Alternative 3'ss",
         fontsize_row = 10,  # Readable row labels
         fontsize_col = 10,
         gaps_row = sequence_change_points)  # Add divider between the two sequences
dev.off()  # Close the PDF device


###################################################
###### Calculate Tau for these splice sites. ######
###################################################
library(matrixStats) 
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


included_by_condition_pivot <- included_by_condition %>%
  ungroup() %>% 
  mutate(index_offset = paste(index, offset, sep = "__")) %>%
  select(index_offset, condition, included_perc) %>% 
  pivot_wider(names_from = condition, values_from = included_perc) 

# Ensure all numeric columns for calculations
numeric_cols <- included_by_condition_pivot %>%
  select(-index_offset) %>%  # Exclude non-numeric columns
  mutate(across(everything(), as.numeric))

# Calculate tau, mean, and median per row
tau_per_row <- apply(numeric_cols, 1, calculate_tau)
mean_per_row <- rowMeans(numeric_cols, na.rm = TRUE)
median_per_row <- rowMedians(as.matrix(numeric_cols), na.rm = TRUE)

# Create dataframe for plotting
tau_per_row_df <- data.frame(
  index_offset = included_by_condition_pivot$index_offset, 
  tau = tau_per_row, 
  mean_per_row = mean_per_row,
  median_per_row = median_per_row
)

tau_per_row_df_filtered <- tau_per_row_df %>% filter(median_per_row > 0.05)

p1 <- ggplot(tau_per_row_df_filtered, aes(x = tau)) +
  geom_histogram(bins = 30, fill = "gray70", color = "black") +  # High contrast bars
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +  # Ensure integer ticks
  labs(title = "Distribution of Tau Values for 3'ss",
       x = "Tau",
       y = "Frequency") +
  theme_classic(base_size = 16) +  # Remove grid lines, improve readability
  theme(
    axis.text.x = element_text(size = 14),  # Readable x-axis
    axis.text.y = element_text(size = 14),  # Readable y-axis
    axis.title = element_text(size = 16),   # Readable axis titles
    plot.title = element_text(size = 18, hjust = 0.5),  # Centered title
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )
ggsave(file.path(output_filepath, "3ss_tau_distribution.pdf"), 
       plot = p1, width = 10, height = 6, dpi = 300)  # Publication quality

# Plot heatmap of sequences that are tau > 0.9
high_tau_seq <- tau_per_row_df %>% filter(tau > 0.8) %>% pull(index_offset)
high_tau_pivot <- included_by_condition_pivot %>% 
  filter(index_offset %in% high_tau_seq) %>% 
  select(-index_offset) %>% 
  mutate(across(everything(), as.numeric))

# Convert to matrix
mat <- as.matrix(high_tau_pivot)
rownames(mat) <- high_tau_seq
# Filter out rows with > 50% NA
mat <- mat[rowMeans(is.na(mat)) < 0.1, ]

pdf(file.path(output_filepath, "3ss_high_tau_heatmap.pdf"), width = 24, height = 16) 
pheatmap(mat,
         cluster_rows = T,  # Maintain alphabetical ordering
         cluster_cols = F,   # Cluster samples (columns)
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = viridis(100),
         main = "Sequences with 3'ss Tau > 0.8",
         fontsize_row = 10,  # Readable row labels
         fontsize_col = 10)
dev.off()
 
 