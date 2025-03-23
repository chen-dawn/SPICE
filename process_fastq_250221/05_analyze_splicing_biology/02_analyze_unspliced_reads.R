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

# Extract all possible (index, sample, offset) combinations
all_distinct_sample <- unique(all_files_df_filtered$sample)

all_distinct_offsets <- unique(all_files_df_filtered$offset)

all_distinct_index <- unique(all_files_df_filtered$index)

num_unique_index_per_sample <- all_files_df_filtered %>% 
  group_by(sample) %>% 
  summarise(num_unique_index = n_distinct(index))

# Compute unspliced percentage
unspliced_only <- all_files_df_filtered %>%
  filter(mode == "UNSPLICED") %>%
  mutate(unspliced_perc = count_scaled / total_count) %>% 
  left_join(num_unique_index_per_sample, by = "sample") 

unspliced_by_sample <- unspliced_only %>% 
  ungroup() %>% 
  group_by(sample, offset) %>% 
  summarise(mean_unspliced_perc = sum(unspliced_perc)/ num_unique_index[1]) 
  

##############################################################
##### General Heatmap Plot ###################################
##############################################################
unspliced_by_sample_pivot <- unspliced_by_sample %>% 
  pivot_wider(names_from = offset, values_from = mean_unspliced_perc, values_fill = 0) %>%
  relocate(order(as.integer(names(.)[-1])))  # Reorder offset columns numerically

# Reorder columns numerically
offsets_numeric <- as.integer(names(unspliced_by_sample_pivot)[-1])
sorted_order <- order(offsets_numeric) + 1  # +1 to account for "sample" column
unspliced_by_sample_pivot <- unspliced_by_sample_pivot[, c(1, sorted_order)]  

# Convert to matrix
unspliced_mat <- as.matrix(unspliced_by_sample_pivot[, -1])  
rownames(unspliced_mat) <- unspliced_by_sample_pivot$sample  

# Save heatmap to high-resolution PDF
pdf(file.path(output_filepath, "unspliced_reads_heatmap.pdf"), width = 20, height = 10)  # Adjust width & height for readability
pheatmap(unspliced_mat, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         show_rownames = TRUE, 
         show_colnames = TRUE,
         main = "Unspliced Reads Percentage",
         fontsize_row = 8,   # Increase for readability
         fontsize_col = 8,   # Increase for clarity
         color = viridis(100))
dev.off()  # Close the PDF device

# Filter out offset == 113. 
unspliced_by_sample_pivot_filtered <- unspliced_by_sample_pivot %>% 
  select(-`113`, -`0`, -`144`)  # Remove the column for offset == 113

# Convert to matrix
unspliced_mat_filtered <- as.matrix(unspliced_by_sample_pivot_filtered[, -1])
rownames(unspliced_mat_filtered) <- unspliced_by_sample_pivot_filtered$sample

# Save heatmap to high-resolution PDF
pdf(file.path(output_filepath, "unspliced_reads_heatmap_filtered.pdf"), width = 20, height = 10)  # Adjust width & height for readability
pheatmap(unspliced_mat_filtered, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         show_rownames = TRUE, 
         show_colnames = TRUE,
         main = "Unspliced Reads Percentage (Filtered)",
         fontsize_row = 8,   # Increase for readability
         fontsize_col = 8,   # Increase for clarity
         color = viridis(100))
dev.off()  # Close the PDF device

# Plot a bar plot.
seq <- "GTGAAGTGATGATGATGGCTGACCAGGCGTTACAGTGTCTCTAGGCAGTTGCTGGGAACTGGCTAGAGACATAAGGTTAAGATGTGAGGAGATGGGTTTTGATTTCTGGACAGGGGAAAGGAAGTAATCTGAGATTGAATCCAGGAAATGA"
# Convert sequence into individual characters
seq_bases <- unlist(strsplit(seq, ""))
unspliced_grouped <- unspliced_by_sample %>%
  group_by(offset) %>%
  summarise(mean_unspliced_perc = mean(mean_unspliced_perc, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    offset = as.integer(offset),  # Ensure numeric ordering
    base = ifelse(offset > 0 & offset <= length(seq_bases), seq_bases[offset], NA)  # Map bases
  ) %>% 
  filter(mean_unspliced_perc > 0) %>%   # Filter out rows with mean_unspliced_perc == 0
  filter(offset >= 30)


num_per_base <- unspliced_grouped %>% group_by(base) %>% 
  summarise(n = n()) %>% 
  # Add a row for each base with n = 0
  complete(base = c("A", "C", "G", "T"), fill = list(n = 0)) 
p1 <- ggplot(num_per_base, aes(x = base, y = n)) +
  geom_bar(stat = "identity", fill = "gray70", color = "black", width = 0.8) +  # Adjust bar width
  labs(title = "Number of Bases per Position",
       x = "Base",
       y = "Count") +
  theme_classic(base_size = 16) +  # Remove grid lines, increase font size
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Rotate x-axis labels
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18, hjust = 0.5),  # Center title
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )

ggsave(file.path(output_filepath, "number_of_bases_per_position.pdf"), plot = p1,
       width = 6, height = 4, dpi = 300)  

##############################################################
##### Are there enriched indexes? ############################
##############################################################
unspliced_by_index <- 
  unspliced_only %>% 
  group_by(index, offset) %>% 
  summarise(mean_unspliced_perc = sum(unspliced_perc)/ n())

high_unspliced <- unspliced_by_index %>% 
  filter(mean_unspliced_perc > 0.5) %>% 
  mutate(offset = as.integer(offset)) %>% 
  filter(offset == 113) %>%
  filter(offset >= 30) %>% 
  separate(index, into = c("ENSG", "refseq_gene", "coordinate"), sep = ";", remove = FALSE) 

num_per_gene <- high_unspliced %>% group_by(refseq_gene) %>% summarise(n = n()) %>% arrange(desc(n))

