library(tidyverse)
library(vroom)
library(data.table)
library(Biostrings)
library(ggpointdensity)
library(pheatmap)

reverse_complement <- function(dna_seq) {
  complement <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
  nucleotides <- unlist(strsplit(dna_seq, ""))
  complement_nucleotides <- complement[nucleotides]
  reverse_complement_seq <- paste(rev(complement_nucleotides), collapse = "")
  return(reverse_complement_seq)
}

out_dir <- "~/Dropbox (Harvard University)/02Splicing/library_47k_missplicing/data_analysis/"
dir.create(out_dir, showWarnings = FALSE)
all_files_df <- fread("~/Dropbox (Harvard University)/02Splicing/latest/umi_count_merged_to_ref_normalized.csv")
K700E_df <- fread("~/Dropbox (Harvard University)/02Splicing/latest/K700E_umi_count_merged_to_ref_normalized.csv")
# nova240826_df <- fread("~/Dropbox (Harvard University)/02Splicing/latest/Nova240826_umi_count_merged_to_ref_normalized.csv")
alt_ref_file <- read_tsv("~/melange/data/guide_library_cleaned/ref_test_alt_ref_dict.tsv")
barcodes <- read_csv("~/melange/data/guide_library/20230130_twist_library_v3.csv") %>%
  mutate(barcodeRevcomp = sapply(barcode, reverse_complement))


##### Let's just check out the K700E samples first #####
# unspliced <- K700E_df %>% filter(mode == "UNSPLICED") %>% 
#   mutate(count = as.integer(count)) %>%
#   group_by(sample, condition, index) %>%
#   mutate(total_count = sum(count)) %>%
#   filter(total_count > 20) %>%
#   ungroup() %>%
#   mutate(fraction = count / total_count) %>%
#   mutate(offset = as.integer(offset))
# ggplot(unspliced, aes(offset, fraction)) + geom_point() + 
#   facet_wrap(~sample, scales = "free_y") + 
#   ggtitle("Unspliced counts for K700E samples")
# 
# # For every index take the average of the unsplice counts for condition.
# unspliced_avg <- unspliced %>% 
#   group_by(index, condition, offset) %>% 
#   summarise(frac_sum = sum(fraction)) %>% 
#   ungroup() %>% 
#   pivot_wider(names_from = condition, values_from = frac_sum, values_fill = 0) 
# 
# # This looks kinda weird. 
# ggplot(unspliced_avg %>% filter(offset == 0), aes(K562_WT, K562_K700E)) + 
#   geom_pointdensity() + 
#   ggtitle("Unspliced counts for K562_WT vs K562_K700E")
# 
# unspliced_avg <- unspliced %>% group_by(condition, offset) %>% 
#   summarise(fraction_sum = sum(fraction)) %>% ungroup()
# 
# ggplot(unspliced_avg, aes(offset, fraction_sum)) + geom_bar(stat = "identity") + 
#   ggtitle("Unspliced counts for K562_WT vs K562_K700E") + facet_wrap(~condition)

# Get included reads only.
included_reads <- K700E_df %>% filter(mode == "INCLUDED") %>% 
  # Filter out read with offset that contain 9999
  # filter(!grepl("9999", offset)) %>% 
  # Separate the offset column into 3 columns by :
  separate(offset, c("offset_mid_start", "offset_mid_end", "offset_down_start"), sep = ":", convert = T, remove = F) %>% 
  filter(offset_down_start == 0) %>%
  group_by(sample, index) %>% 
  mutate(total_sum = sum(count)) %>% 
  ungroup() %>%
  filter(total_sum > 30) %>%
  mutate(fraction = count / total_sum) %>%
  group_by(condition, index, offset_mid_start, offset_mid_end, offset_down_start, offset) %>%
  summarise(fraction = mean(fraction)) 


# Pivot wider by count. 
included_reads_wide <- included_reads %>% 
  pivot_wider(names_from = condition, values_from = fraction, values_fill = NA) %>% 
  mutate(diff = K562_K700E - K562_WT) %>% 
  mutate(log_ratio = log2(K562_K700E / K562_WT))

high_in_K700E <- included_reads_wide %>% 
  # filter(offset_mid_start ==0 & offset_mid_end == 0 & offset_down_start == 0) %>% 
  # filter(diff < -0.3 & K562_WT > 0.9) %>% 
  # arrange(diff)
  filter(offset != "0:0:0") %>% 
  # filter(offset_mid_start == 0 & offset_mid_end == 0 & offset_down_start == 0) %>% 
  filter(abs(diff) > 0.1 & abs(log_ratio) > 2) %>% 
  arrange(desc(log_ratio))

# Get the ones that are high in K700E and check the entire counts.
shortlisted_elements <- included_reads_wide %>% 
  filter(index %in% high_in_K700E$index) %>%
  pivot_longer(cols = c(K562_WT, K562_K700E), names_to = "condition", values_to = "fraction") %>% 
  mutate(gene_name = str_split(index, ";")[[1]][2])

ggplot(shortlisted_elements, aes(offset, fraction, fill = condition)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  facet_wrap(~gene_name, scales = "free_x") + 
  ggtitle("High in K700E, low in WT") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# test <-shortlisted_elements %>% group_by(index, condition) %>% summarise(total = sum(fraction, na.rm = T))
# # Get the barcodes.
# shortlisted_barcodes <- merge(shortlisted_elements, barcodes, by.x = "index", by.y = "ID")
# write_csv(shortlisted_barcodes, file.path(out_dir, "high_WT_shortlisted_barcodes.csv"))
# write_csv(shortlisted_barcodes, "/Volumes/broad_dawnccle/processed_data/K700E_shortlist/high_WT_shortlisted_elements.csv")

# # Check out 2 random samples. K562_1ugNuc and K562_2ugNuc conditions.
# set1 <- all_files_df %>% filter(condition %in% c("K562_1ugNuc", "K562_2ugNuc")) %>% 
#   filter(mode == "INCLUDED") %>% 
#   filter(!grepl("9999", offset)) %>% 
#   separate(offset, c("offset_mid_start", "offset_mid_end", "offset_down_start"), sep = ":", convert = T, remove = F) %>% 
#   filter(offset_down_start == 0) %>%
#   group_by(sample, index) %>% 
#   mutate(total_sum = sum(count)) %>% 
#   ungroup() %>%
#   filter(total_sum > 20) %>%
#   mutate(fraction = count / total_sum) %>%
#   group_by(condition, index, offset_mid_start, offset_mid_end, offset_down_start, offset) %>%
#   summarise(fraction = mean(fraction))
# 
# set1_wide <- set1 %>%
#   pivot_wider(names_from = condition, values_from = fraction, values_fill = NA) %>% 
#   mutate(diff = K562_2ugNuc - K562_1ugNuc)
# 
# high_in_2ugNuc <- set1_wide %>%
#   filter(offset_mid_start ==0 & offset_mid_end == 0 & offset_down_start == 0) %>% 
#   filter(abs(diff)> 0.3) %>% 
#   arrange(diff)
# 
# # ok there's nothing. Good.
# shortlisted_elements_2ugNuc <- set1_wide %>%
#   filter(index %in% high_in_2ugNuc$index) %>%
#   pivot_longer(cols = c(K562_1ugNuc, K562_2ugNuc), names_to = "condition", values_to = "fraction") %>% 
#   mutate(gene_name = str_split(index, ";")[[1]][2])
# 
# ggplot(shortlisted_elements_2ugNuc, aes(offset, fraction, fill = condition)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   facet_wrap(~gene_name, scales = "free_x") +
#   ggtitle("K562 2ugNuc, 1ugNuc, abs diff >0.3") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))
# 
# # check out 8MGBA_Nuc and 8MGBA_tfx
# set2 <- all_files_df %>% filter(condition %in% c("8MGBA_Nuc", "8MGBA_tfx")) %>% 
#   filter(mode == "INCLUDED") %>% 
#   filter(!grepl("9999", offset)) %>% 
#   separate(offset, c("offset_mid_start", "offset_mid_end", "offset_down_start"), sep = ":", convert = T, remove = F) %>% 
#   filter(offset_down_start == 0) %>%
#   group_by(sample, index) %>% 
#   mutate(total_sum = sum(count)) %>% 
#   ungroup() %>%
#   filter(total_sum > 20) %>%
#   mutate(fraction = count / total_sum) %>%
#   group_by(condition, index, offset_mid_start, offset_mid_end, offset_down_start, offset) %>%
#   summarise(fraction = mean(fraction))
# 
# set2_wide <- set2 %>%
#   pivot_wider(names_from = condition, values_from = fraction, values_fill = NA) %>% 
#   mutate(diff = `8MGBA_tfx` - `8MGBA_Nuc`)
# 
# high_in_tfx <- set2_wide %>%
#   filter(offset_mid_start ==0 & offset_mid_end == 0 & offset_down_start == 0) %>% 
#   filter(abs(diff) > 0.5) %>% 
#   arrange(diff)
# 
# shortlisted_elements_tfx <- set2_wide %>%
#   filter(index %in% high_in_tfx$index) %>%
#   pivot_longer(cols = c(`8MGBA_Nuc`, `8MGBA_tfx`), names_to = "condition", values_to = "fraction") %>% 
#   mutate(gene_name = str_split(index, ";")[[1]][2])
# 
# ggplot(shortlisted_elements_tfx, aes(offset, fraction, fill = condition)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   facet_wrap(~gene_name, scales = "free_x") +
#   ggtitle("High in 8MGBA_tfx, low in 8MGBA_Nuc") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))
# 
# # Check out HEK vs K562_2ugNuc
# set3 <- all_files_df %>% filter(condition %in% c("HEK", "K562_2ugNuc")) %>% 
#   filter(mode == "INCLUDED") %>% 
#   filter(!grepl("9999", offset)) %>% 
#   separate(offset, c("offset_mid_start", "offset_mid_end", "offset_down_start"), sep = ":", convert = T, remove = F) %>% 
#   filter(offset_down_start == 0) %>%
#   group_by(sample, index) %>% 
#   mutate(total_sum = sum(count)) %>% 
#   ungroup() %>%
#   filter(total_sum > 20) %>%
#   mutate(fraction = count / total_sum) %>%
#   group_by(condition, index, offset_mid_start, offset_mid_end, offset_down_start, offset) %>%
#   summarise(fraction = mean(fraction))
# 
# set3_wide <- set3 %>%
#   pivot_wider(names_from = condition, values_from = fraction, values_fill = NA) %>% 
#   mutate(diff = K562_2ugNuc - HEK)
# 
# high_in_K562 <- set3_wide %>%
#   filter(offset_mid_start ==0 & offset_mid_end == 0 & offset_down_start == 0) %>% 
#   filter(abs(diff) > 0.4) %>% 
#   arrange(diff)
# 
# shortlisted_elements_K562 <- set3_wide %>%
#   filter(index %in% high_in_K562$index) %>%
#   pivot_longer(cols = c(HEK, K562_2ugNuc), names_to = "condition", values_to = "fraction") %>% 
#   mutate(gene_name = str_split(index, ";")[[1]][2])
# 
# ggplot(shortlisted_elements_K562, aes(offset, fraction, fill = condition)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   facet_wrap(~gene_name, scales = "free_x") +
#   ggtitle("K562_2ugNuc/ HEK diff > 0.4") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))
# 
# # Check out T47D vs Kelly
# set4 <- all_files_df %>% filter(condition %in% c("T47D", "Kelly")) %>% 
#   filter(mode == "INCLUDED") %>% 
#   filter(!grepl("9999", offset)) %>% 
#   separate(offset, c("offset_mid_start", "offset_mid_end", "offset_down_start"), sep = ":", convert = T, remove = F) %>% 
#   filter(offset_down_start == 0) %>%
#   group_by(sample, index) %>% 
#   mutate(total_sum = sum(count)) %>% 
#   ungroup() %>%
#   filter(total_sum > 20) %>%
#   mutate(fraction = count / total_sum) %>%
#   group_by(condition, index, offset_mid_start, offset_mid_end, offset_down_start, offset) %>%
#   summarise(fraction = mean(fraction))
# 
# set4_wide <- set4 %>%
#   pivot_wider(names_from = condition, values_from = fraction, values_fill = NA) %>% 
#   mutate(diff = Kelly - T47D)
# 
# high_in_Kelly <- set4_wide %>%
#   filter(offset_mid_start ==0 & offset_mid_end == 0 & offset_down_start == 0) %>% 
#   filter(abs(diff) > 0.5) %>% 
#   arrange(diff)
# 
# shortlisted_elements_Kelly <- set4_wide %>%
#   filter(index %in% high_in_Kelly$index) %>%
#   pivot_longer(cols = c(T47D, Kelly), names_to = "condition", values_to = "fraction") %>% 
#   mutate(gene_name = str_split(index, ";")[[1]][2])
# 
# ggplot(shortlisted_elements_Kelly, aes(offset, fraction, fill = condition)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   facet_wrap(~gene_name, scales = "free_x") +
#   ggtitle("Kelly/ T47D diff > 0.4") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Do for all samples.
all_samples <- all_files_df %>% 
  filter(mode == "INCLUDED") %>% 
  # filter(!grepl("9999", offset)) %>% 
  separate(offset, c("offset_mid_start", "offset_mid_end", "offset_down_start"), sep = ":", convert = T, remove = F) %>% 
  filter(offset_down_start == 0) %>%
  group_by(sample, index) %>% 
  mutate(total_sum = sum(count)) %>% 
  ungroup() %>%
  filter(total_sum > 20) %>%
  mutate(fraction = count / total_sum) %>%
  group_by(condition, index, offset_mid_start, offset_mid_end, offset_down_start, offset) %>%
  summarise(fraction = mean(fraction))

fwrite(all_samples, file.path(out_dir, "all_samples_cleaned_with_good_upstream_downstream_const.csv"))
all_samples <- fread(file.path(out_dir, "all_samples_cleaned_with_good_upstream_downstream_const.csv"))
all_samples_wide <- all_samples %>%
  filter(offset_mid_start == 0 & offset_mid_end == 0 & offset_down_start == 0) %>% 
  filter(!condition %in% c("OSRC2", "Kelly_old", "SKNAS_Nuc", "A172")) %>% 
  # Set NA values to -1.
  pivot_wider(names_from = condition, values_from = fraction, values_fill = -1) 

all_samples_wide_matrix <- all_samples_wide %>%
  ungroup() %>%
  select(-index, -offset_mid_start, -offset_mid_end, -offset_down_start, -offset) %>%
  as.matrix()

# Set the row names.
rownames(all_samples_wide_matrix) <- all_samples %>%
  filter(offset_mid_start == 0 & offset_mid_end == 0 & offset_down_start == 0) %>% 
  filter(!condition %in% c("OSRC2", "Kelly_old", "SKNAS_Nuc", "A172")) %>% 
  # Set NA values to -1.
  pivot_wider(names_from = condition, values_from = fraction, values_fill = -1) %>%
  pull(index)


pheatmap(
  all_samples_wide_matrix,
  cluster_rows = T,
  cluster_cols = T,
  color = colorRampPalette(c("black", "black", "#45abd7", "white", "#cf1c1a"))(100),
  show_rownames = F
)


#### Let's label the NA as NA.
all_samples_wide_matrix_NA <- all_samples %>%
  filter(offset_mid_start == 0 & offset_mid_end == 0 & offset_down_start == 0) %>% 
  filter(!condition %in% c("OSRC2", "Kelly_old", "SKNAS_Nuc", "A172")) %>% 
  pivot_wider(names_from = condition, values_from = fraction, values_fill = NA) %>%
  ungroup() %>%
  select(-index, -offset_mid_start, -offset_mid_end, -offset_down_start, -offset) %>%
  as.matrix()

# Set the row names.
rownames(all_samples_wide_matrix_NA) <- all_samples %>%
  filter(offset_mid_start == 0 & offset_mid_end == 0 & offset_down_start == 0) %>% 
  filter(!condition %in% c("OSRC2", "Kelly_old", "SKNAS_Nuc", "A172")) %>% 
  # Set NA values to -1.
  pivot_wider(names_from = condition, values_from = fraction, values_fill = NA) %>%
  pull(index)

# Get the average of each row.
average_each_row <- rowMeans(all_samples_wide_matrix_NA, na.rm = T)
# Get the min of each row.
min_each_row <- apply(all_samples_wide_matrix_NA, 1, min, na.rm = T)

# Make a dataframe.
average_each_row_df <- data.frame(index = rownames(all_samples_wide_matrix_NA), average = average_each_row, min = min_each_row)
seq_shortlist <- average_each_row_df %>% filter(average >0.9 & min < 0.5) 

# filter all_samples_wide_matrix
all_samples_wide_matrix_shortlist <- all_samples_wide_matrix[seq_shortlist$index, ]
# Heatmap
pheatmap(
  all_samples_wide_matrix_shortlist,
  cluster_rows = T,
  cluster_cols = T,
  color = colorRampPalette(c("black", "black", "#45abd7", "white", "#cf1c1a"))(100),
  show_rownames = F
)

# Plot maps for all of them individually. 
shortlisted_sequences <- seq_shortlist$index
shortlisted_sequences <- unique(top_seq$ExonID)
for(seq in shortlisted_sequences){
  print(seq)
  gene_name <- str_split(seq, ";")[[1]][2]
  chr_pos <- str_split(seq, ";")[[1]][3]
  tmp_filename <- paste0("V5_", gene_name, "_", chr_pos, ".pdf")
  tmp_sample <- all_samples %>% filter(index == seq) 
  plot <- ggplot(tmp_sample, aes(offset, fraction, fill = offset == "0:0:0")) + 
    geom_bar(stat = "identity", position = "dodge") + 
    theme_bw() +
    scale_fill_manual(values = c("FALSE" = "red", "TRUE" = "grey60")) + 
    facet_wrap(~condition, scales = "free_x") + 
    ggtitle(seq) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none")
  
  ggsave(filename = file.path(out_dir, tmp_filename), plot = plot, dpi = 300, width = 4000, height = 3000, units = "px")
}

# Find the shortlisted sequences in the twist library.
twist_library <- read_csv("~/melange/data/guide_library/20230130_twist_library_v3.csv") %>%
  mutate(barcodeRevcomp = sapply(barcode, reverse_complement)) %>% 
  mutate(index  = ID)
missplicing_barcodes <- twist_library %>% filter(ID %in% shortlisted_sequences)
# Create dir
dir.create("/Volumes/broad_dawnccle/processed_data/KMRC1_shortlist", showWarnings = FALSE)
write_csv(missplicing_barcodes, "/Volumes/broad_dawnccle/processed_data/KMRC1_shortlist/high_KMRC1_shortlisted_elements.csv")

##### Take the all_samples_wide_matrix_NA and calculate tau score.
# Take 1 minus all value.
all_samples_wide_matrix_NA_1_minus <- 1 - all_samples_wide_matrix_NA
calculate_tau <- function(row) {
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

# Apply the function to each row of the matrix
tau_score <- apply(all_samples_wide_matrix_NA_1_minus, 1, calculate_tau)

# Plot a histogram of tau score.
ggplot(data.frame(tau_score = tau_score), aes(tau_score)) + geom_histogram(bins = 100) + ggtitle("Tau score histogram")


# Get the elements with tau score = 1 and plot heatmap.
tau_score_df <- data.frame(index = rownames(all_samples_wide_matrix_NA_1_minus), tau_score = tau_score)
tau_score_1 <- tau_score_df %>% filter(tau_score >0.6)
all_samples_wide_matrix_tau_1 <- all_samples_wide_matrix[tau_score_1$index, ]

pheatmap(
  all_samples_wide_matrix_tau_1,
  cluster_rows = T,
  cluster_cols = F,
  color = colorRampPalette(c("black", "black", "#45abd7", "white", "#cf1c1a"))(100),
  show_rownames = T, 
  # Add main title
  main = "Upsilon score > 0.6, filtered V5 matrix"
)

# Plot maps for all of them individually.
shortlisted_sequences_tau <- tau_score_1$index

for(seq in shortlisted_sequences_tau){
  print(seq)
  gene_name <- str_split(seq, ";")[[1]][2]
  chr_pos <- str_split(seq, ";")[[1]][3]
  tmp_filename <- paste0("V5_", gene_name, "_", chr_pos, "_tau.pdf")
  tmp_sample <- all_samples %>% filter(index == seq) 
  plot <- ggplot(tmp_sample, aes(offset, fraction, fill = offset == "0:0:0")) + 
    geom_bar(stat = "identity", position = "dodge") + 
    theme_bw() +
    scale_fill_manual(values = c("FALSE" = "red", "TRUE" = "grey60")) + 
    facet_wrap(~condition, scales = "free_x") + 
    ggtitle(seq) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none")
  
  ggsave(filename = file.path(out_dir, tmp_filename), plot = plot, dpi = 300, width = 4000, height = 3000, units = "px")
}

# Find the shortlisted sequences in the twist library.
missplicing_barcodes_tau <- twist_library %>% filter(ID %in% shortlisted_sequences_tau)
write_csv(missplicing_barcodes_tau, "/Volumes/broad_dawnccle/processed_data/KMRC1_shortlist/high_KMRC1_shortlisted_elements_tau.csv")

###### Look at sample replicates ######
all_sample_reps <- all_files_df %>% 
  filter(mode == "INCLUDED") %>% 
  # filter(!grepl("9999", offset)) %>% 
  separate(offset, c("offset_mid_start", "offset_mid_end", "offset_down_start"), sep = ":", convert = T, remove = F) %>% 
  filter(offset_down_start == 0) %>%
  group_by(sample, index) %>% 
  mutate(total_sum = sum(count)) %>% 
  ungroup() %>%
  filter(total_sum > 20) %>%
  mutate(fraction = count / total_sum) %>%
  group_by(sample, condition, index, offset_mid_start, offset_mid_end, offset_down_start, offset) %>%
  summarise(fraction = mean(fraction))

fwrite(all_sample_reps, file.path(out_dir, "all_samples_cleaned_with_good_upstream_downstream_const_sample_granularity.csv"))

all_samples_wide <- all_sample_reps %>%
  filter(offset_mid_start == 0 & offset_mid_end == 0 & offset_down_start == 0) %>% 
  filter(!condition %in% c("OSRC2", "Kelly_old", "SKNAS_Nuc", "A172")) %>% 
  # Set NA values to -1.
  pivot_wider(names_from = sample, values_from = fraction, values_fill = -1) 

all_samples_wide_matrix <- all_samples_wide %>%
  ungroup() %>%
  select(-index, -offset_mid_start, -offset_mid_end, -offset_down_start, -offset, -condition) %>%
  as.matrix()

# Set the row names.
rownames(all_samples_wide_matrix) <- all_sample_reps %>%
  filter(offset_mid_start == 0 & offset_mid_end == 0 & offset_down_start == 0) %>% 
  filter(!condition %in% c("OSRC2", "Kelly_old", "SKNAS_Nuc", "A172")) %>% 
  # Set NA values to -1.
  pivot_wider(names_from = sample, values_from = fraction, values_fill = -1) %>%
  pull(index)


#### Let's label the NA as NA.
all_samples_wide_matrix_NA <- all_sample_reps %>%
  filter(offset_mid_start == 0 & offset_mid_end == 0 & offset_down_start == 0) %>% 
  filter(!condition %in% c("OSRC2", "Kelly_old", "SKNAS_Nuc", "A172")) %>% 
  pivot_wider(names_from = sample, values_from = fraction, values_fill = NA) %>%
  ungroup() %>%
  select(-index, -offset_mid_start, -offset_mid_end, -offset_down_start, -offset, -condition) %>%
  as.matrix()

# Set the row names.
rownames(all_samples_wide_matrix_NA) <- all_sample_reps %>%
  filter(offset_mid_start == 0 & offset_mid_end == 0 & offset_down_start == 0) %>% 
  filter(!condition %in% c("OSRC2", "Kelly_old", "SKNAS_Nuc", "A172")) %>% 
  # Set NA values to -1.
  pivot_wider(names_from = sample, values_from = fraction, values_fill = NA) %>%
  pull(index)

# Get the average of each row.
average_each_row <- rowMeans(all_samples_wide_matrix_NA, na.rm = T)
# Get the min of each row.
min_each_row <- apply(all_samples_wide_matrix_NA, 1, min, na.rm = T)

# Make a dataframe.
average_each_row_df <- data.frame(index = rownames(all_samples_wide_matrix_NA), average = average_each_row, min = min_each_row)
seq_shortlist <- average_each_row_df %>% filter(average >0.9 & min < 0.5) 

# filter all_samples_wide_matrix
all_samples_wide_matrix_shortlist <- all_samples_wide_matrix[seq_shortlist$index, ]
# Heatmap
pheatmap(
  all_samples_wide_matrix_shortlist,
  cluster_rows = T,
  cluster_cols = T,
  color = colorRampPalette(c("black", "black", "#45abd7", "white", "#cf1c1a"))(100),
  show_rownames = F
)

# Plot maps for all of them individually. 
shortlisted_sequences <- seq_shortlist$index

for(seq in shortlisted_sequences){
  print(seq)
  gene_name <- str_split(seq, ";")[[1]][2]
  chr_pos <- str_split(seq, ";")[[1]][3]
  tmp_filename <- paste0("V5_", gene_name, "_", chr_pos, ".pdf")
  tmp_sample <- all_samples %>% filter(index == seq) 
  plot <- ggplot(tmp_sample, aes(offset, fraction, fill = offset == "0:0:0")) + 
    geom_bar(stat = "identity", position = "dodge") + 
    theme_bw() +
    scale_fill_manual(values = c("FALSE" = "red", "TRUE" = "grey60")) + 
    facet_wrap(~condition, scales = "free_x") + 
    ggtitle(seq) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none")
  
  ggsave(filename = file.path(out_dir, tmp_filename), plot = plot, dpi = 300, width = 4000, height = 3000, units = "px")
}

# Find the shortlisted sequences in the twist library.
twist_library <- read_csv("~/melange/data/guide_library/20230130_twist_library_v3.csv") %>%
  mutate(barcodeRevcomp = sapply(barcode, reverse_complement)) %>% 
  mutate(index  = ID)
missplicing_barcodes <- twist_library %>% filter(ID %in% shortlisted_sequences)
# Create dir
dir.create("/Volumes/broad_dawnccle/processed_data/KMRC1_shortlist", showWarnings = FALSE)
write_csv(missplicing_barcodes, "/Volumes/broad_dawnccle/processed_data/KMRC1_shortlist/high_KMRC1_shortlisted_elements.csv")

##### Take the all_samples_wide_matrix_NA and calculate tau score.
# Take 1 minus all value.
all_samples_wide_matrix_NA_1_minus <- 1 - all_samples_wide_matrix_NA
calculate_tau <- function(row) {
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

# Apply the function to each row of the matrix
tau_score <- apply(all_samples_wide_matrix_NA_1_minus, 1, calculate_tau)

# Plot a histogram of tau score.
ggplot(data.frame(tau_score = tau_score), aes(tau_score)) + geom_histogram(bins = 100) + ggtitle("Tau score histogram")


# Get the elements with tau score = 1 and plot heatmap.
tau_score_df <- data.frame(index = rownames(all_samples_wide_matrix_NA_1_minus), tau_score = tau_score)
tau_score_1 <- tau_score_df %>% filter(tau_score >0.6)
all_samples_wide_matrix_tau_1 <- all_samples_wide_matrix[tau_score_1$index, ]

pheatmap(
  all_samples_wide_matrix_tau_1,
  cluster_rows = T,
  cluster_cols = F,
  color = colorRampPalette(c("black", "black", "#45abd7", "white", "#cf1c1a"))(100),
  show_rownames = T, 
  # Add main title
  main = "Upsilon score > 0.6, filtered V5 matrix"
)

# Plot maps for all of them individually.
shortlisted_sequences_tau <- tau_score_1$index

for(seq in shortlisted_sequences_tau){
  print(seq)
  gene_name <- str_split(seq, ";")[[1]][2]
  chr_pos <- str_split(seq, ";")[[1]][3]
  tmp_filename <- paste0("V5_", gene_name, "_", chr_pos, "_tau.pdf")
  tmp_sample <- all_samples %>% filter(index == seq) 
  plot <- ggplot(tmp_sample, aes(offset, fraction, fill = offset == "0:0:0")) + 
    geom_bar(stat = "identity", position = "dodge") + 
    theme_bw() +
    scale_fill_manual(values = c("FALSE" = "red", "TRUE" = "grey60")) + 
    facet_wrap(~condition, scales = "free_x") + 
    ggtitle(seq) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none")
  
  ggsave(filename = file.path(out_dir, tmp_filename), plot = plot, dpi = 300, width = 4000, height = 3000, units = "px")
}

