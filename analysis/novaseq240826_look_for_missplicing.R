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

out_dir <- "~/Dropbox (Harvard University)/02Splicing/SplicingManuscript/figure_outputs/"
dir.create(out_dir, showWarnings = FALSE)
# all_files_df <- fread("~/Dropbox (Harvard University)/02Splicing/latest/umi_count_merged_to_ref_normalized.csv")
# K700E_df <- fread("~/Dropbox (Harvard University)/02Splicing/latest/K700E_umi_count_merged_to_ref_normalized.csv")
nova240826_df <- fread("~/Dropbox (Harvard University)/02Splicing/latest/Nova240826_umi_count_merged_to_ref_normalized.csv")
alt_ref_file <- read_tsv("~/melange/data/guide_library_cleaned/ref_test_alt_ref_dict.tsv")
barcodes <- read_csv("~/melange/data/guide_library/20230130_twist_library_v3.csv") %>%
  mutate(barcodeRevcomp = sapply(barcode, reverse_complement))


# Normalize every sample so the total count is 50 million.
# EDIT I find that maybe i shouldn't do this because this is UMI counts? Not clear. 
nova240826_df <- nova240826_df %>% 
  group_by(sample) %>% 
  mutate(total_count = sum(count)) %>% 
  ungroup() %>% 
  mutate(count = count * 5e7 / total_count) %>% 
  mutate(count = as.integer(count)) %>%
  select(-total_count)


# Get included reads only.
included_reads <- nova240826_df %>% filter(mode == "INCLUDED") %>% 
  filter(count > 5) %>% 
  filter(condition %in% c("MUT-sgCh3-1", "MUT-sgRUBP1")) %>% 
  # filter(condition %in% c("MUT-plx-LacZ", "MUT-sgRUBP1")) %>%
  # Filter out read with offset that contain 9999
  # filter(!grepl("9999", offset)) %>% 
  # Separate the offset column into 3 columns by :
  separate(offset, c("offset_mid_start", "offset_mid_end", "offset_down_start"), sep = ":", convert = T, remove = F) %>% 
  filter(offset_down_start == 0) %>%
  group_by(sample, index) %>% 
  mutate(total_sum = sum(count)) %>% 
  ungroup() %>%
  filter(total_sum > 20) %>%
  mutate(fraction = count / total_sum) %>%
  group_by(condition, index, offset_mid_start, offset_mid_end, offset_down_start, offset) %>%
  summarise(fraction = mean(fraction)) 


# Pivot wider by count. 
included_reads_wide <- included_reads %>% 
  pivot_wider(names_from = condition, values_from = fraction, values_fill = NA) %>% 
  mutate(diff = `MUT-sgRUBP1` - `MUT-sgCh3-1`) %>% 
  mutate(ratio = `MUT-sgRUBP1` / `MUT-sgCh3-1`) %>% 
  mutate(log_ratio = log2(ratio))

high_in_RUPB1 <- included_reads_wide %>% 
  filter(offset != "0:0:0") %>% 
  # filter(offset_mid_start == 0 & offset_mid_end == 0 & offset_down_start == 0) %>% 
  filter(abs(diff) > 0.1 & abs(log_ratio) > 2) %>% 
  arrange(desc(ratio)) 

# Get the ones that are high in RUBP1 mut and check the entire counts.
shortlisted_elements <- included_reads_wide %>% 
  filter(index %in% high_in_RUPB1$index) %>%
  pivot_longer(cols = c(`MUT-sgCh3-1`, `MUT-sgRUBP1`), names_to = "condition", values_to = "fraction") %>% 
  mutate(gene_name = str_split(index, ";")[[1]][2])

ggplot(shortlisted_elements, aes(offset, fraction, fill = condition)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  facet_wrap(~gene_name, scales = "free_x") + 
  ggtitle("High in MUT, low in WT") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Get the sequences in all cells.
shortlisted_elements_all_modes <- nova240826_df %>% 
  filter(index %in% high_in_RUPB1$index) %>% 
  filter(condition %in% c("MUT-sgCh3-1", "MUT-sgRUBP1")) 


# Calculate the unspliced rate for each sequence.
unspliced_rate <- shortlisted_elements_all_modes %>% 
  group_by(sample, index, condition) %>% 
  mutate(total_count = sum(count)) %>% 
  ungroup() %>% 
  filter(total_count > 20) %>% 
  mutate(fraction = count / total_count) %>% 
  filter(mode == "UNSPLICED") %>%
  group_by(condition, index, offset) %>%
  summarise(fraction = mean(fraction)) %>%
  pivot_wider(names_from = condition, values_from = fraction, values_fill = NA) %>%
  mutate(diff = `MUT-sgRUBP1` - `MUT-sgCh3-1`) %>%
  mutate(ratio = `MUT-sgRUBP1` / `MUT-sgCh3-1`) %>% 
  filter(abs(diff) > 0.1) 
  
###### Do the same plot for MUT-plx317_U2AF1_WT vs MUT-plx317_U2AF1_Q157A #####
included_reads <- nova240826_df %>% filter(mode == "INCLUDED") %>% 
  # filter(condition %in% c("MUT-plx317_U2AF1_WT", "MUT-plx317_U2AF1_Q157A")) %>%
  filter(count > 5) %>% 
  filter(grepl("MUT-", condition)) %>%
  separate(offset, c("offset_mid_start", "offset_mid_end", "offset_down_start"), sep = ":", convert = T, remove = F) %>% 
  filter(offset_down_start == 0) %>%
  group_by(sample, index) %>% 
  mutate(total_sum = sum(count)) %>% 
  ungroup() %>%
  filter(total_sum > 30) %>%
  mutate(fraction = count / total_sum) %>%
  group_by(condition, index, offset_mid_start, offset_mid_end, offset_down_start, offset) %>%
  summarise(fraction = mean(fraction)) 

###### Look at U2AF1_Q157A #####
included_reads_wide <- included_reads %>% 
  filter(condition %in% c("MUT-plx317_U2AF1_WT", "MUT-plx317_U2AF1_Q157A")) %>%
  pivot_wider(names_from = condition, values_from = fraction, values_fill = NA) %>% 
  mutate(diff = `MUT-plx317_U2AF1_Q157A` - `MUT-plx317_U2AF1_WT`) %>% 
  mutate(ratio = `MUT-plx317_U2AF1_Q157A` / `MUT-plx317_U2AF1_WT`) %>% 
  mutate(log_ratio = log2(ratio))

high_in_U2AF1 <- included_reads_wide %>%
  filter(offset != "0:0:0") %>% 
  # filter(offset_mid_start == 0 & offset_mid_end == 0 & offset_down_start == 0) %>% 
  filter(abs(diff) > 0.1 & abs(log_ratio) > 2) %>% 
  arrange(desc(ratio))

shortlisted_elements <- included_reads_wide %>%
  filter(index %in% high_in_U2AF1$index) %>%
  pivot_longer(cols = c(`MUT-plx317_U2AF1_WT`, `MUT-plx317_U2AF1_Q157A`), names_to = "condition", values_to = "fraction") %>% 
  mutate(gene_name = str_split(index, ";")[[1]][2])

ggplot(shortlisted_elements, aes(offset, fraction, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~gene_name, scales = "free_x") +
  ggtitle("High in U2AF1_Q157A, low in WT") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

##### Compare MUT-sgRBM10 to MUT-sgCh3-1 #####
included_reads_wide <- included_reads %>% 
  filter(condition %in% c("MUT-sgRBM10", "MUT-sgCh3-1")) %>%
  pivot_wider(names_from = condition, values_from = fraction, values_fill = NA) %>% 
  mutate(diff = `MUT-sgRBM10` - `MUT-sgCh3-1`) %>% 
  mutate(ratio = `MUT-sgRBM10` / `MUT-sgCh3-1`) %>% 
  mutate(log_ratio = log2(ratio))

high_in_RBM10 <- included_reads_wide %>%
  filter(offset != "0:0:0") %>% 
  # filter(offset_mid_start == 0 & offset_mid_end == 0 & offset_down_start == 0) %>% 
  filter(abs(diff) > 0.1 & abs(log_ratio) > 2) %>% 
  arrange(desc(ratio))

shortlisted_elements <- included_reads_wide %>%
  filter(index %in% high_in_RBM10$index) %>%
  pivot_longer(cols = c(`MUT-sgCh3-1`, `MUT-sgRBM10`), names_to = "condition", values_to = "fraction") %>% 
  mutate(gene_name = str_split(index, ";")[[1]][2])

ggplot(shortlisted_elements, aes(offset, fraction, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~gene_name, scales = "free_x") +
  ggtitle("High in RBM10, low in WT") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

##### Compare MUT-sgRBM5 to MUT-sgCh3-1 #####
included_reads_wide <- included_reads %>% 
  filter(condition %in% c("MUT-sgRBM5", "MUT-sgCh3-1")) %>%
  pivot_wider(names_from = condition, values_from = fraction, values_fill = NA) %>% 
  mutate(diff = `MUT-sgRBM5` - `MUT-sgCh3-1`) %>% 
  mutate(ratio = `MUT-sgRBM5` / `MUT-sgCh3-1`) %>% 
  mutate(log_ratio = log2(ratio))

high_in_RBM5 <- included_reads_wide %>%
  filter(offset != "0:0:0") %>% 
  # filter(offset_mid_start == 0 & offset_mid_end == 0 & offset_down_start == 0) %>% 
  filter(abs(diff) > 0.1 & abs(log_ratio) > 2) %>% 
  arrange(desc(ratio))

shortlisted_elements <- included_reads_wide %>%
  filter(index %in% high_in_RBM5$index) %>%
  pivot_longer(cols = c(`MUT-sgCh3-1`, `MUT-sgRBM5`), names_to = "condition", values_to = "fraction") %>% 
  mutate(gene_name = str_split(index, ";")[[1]][2])

ggplot(shortlisted_elements, aes(offset, fraction, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~gene_name, scales = "free_x") +
  ggtitle("High in RBM5, low in WT") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
        
# get raw counts
raw_counts <- nova240826_df %>% 
  filter(condition %in% c("MUT-sgRBM5", "MUT-sgCh3-1")) %>%
  filter(index %in% high_in_RBM5$index) 

####### Get the different modes in different samples #######
different_modes <- nova240826_df %>% 
  filter(grepl("MUT", condition)) %>% 
  # filter(condition == "MUT-sgRBM5") %>%
  # filter(count > 5) %>% 
  group_by(condition, index, mode, offset) %>% 
  summarise(count = sum(count)) %>% 
  ungroup()

# sgCh3-1 modes.
Ch3_modes <- different_modes %>% 
  filter(condition == "MUT-sgCh3-1") %>% 
  dplyr::rename(count_Ch3 = count) %>% 
  select(-condition)
############################
# RUBP1 modes. 
RUBP1_modes <- different_modes %>% 
  filter(condition == "MUT-sgRUBP1") %>% 
  dplyr::rename(count_RUBP1 = count) %>% 
  select(-condition)

# Get modes not in sgCh3-1 but in sgRUBP1
RUBP1_to_control <- merge(Ch3_modes,RUBP1_modes, by = c("index", "mode", "offset"), all = T) %>% 
  filter(is.na(count_Ch3) | is.na(count_RUBP1)) %>% 
  filter(!grepl("9999", offset))
############################
# Compare RBM5 to control
RBM5_modes <- different_modes %>% 
  filter(condition == "MUT-sgRBM5") %>% 
  dplyr::rename(count_RBM5 = count) %>% 
  select(-condition)

# Get modes not in sgCh3-1 but in sgRBM5
RBM5_to_control <- merge(Ch3_modes,RBM5_modes, by = c("index", "mode", "offset"), all = T) %>% 
  filter(is.na(count_Ch3) | is.na(count_RBM5)) %>% 
  filter(!grepl("9999", offset))
############################
# Compare RBM10 to control
RBM10_modes <- different_modes %>% 
  filter(condition == "MUT-sgRBM10") %>% 
  dplyr::rename(count_RBM10 = count) %>% 
  select(-condition)

# Get modes not in sgCh3-1 but in sgRBM10
RBM10_to_control <- merge(Ch3_modes,RBM10_modes, by = c("index", "mode", "offset"), all = T) %>% 
  filter(is.na(count_Ch3) | is.na(count_RBM10)) %>% 
  filter(!grepl("9999", offset)) %>% 
  arrange(desc(count_RBM10))

# Check out ENSG00000080298.16;RFX3;chr9-3228846-3228889-3218296-3225280-3248031-3248185 0:-8:4
# This is in RBM10 but not in control 
top_10_in_RBM10 <- RBM10_to_control %>% 
  arrange(desc(count_RBM10)) %>%
  head(12) 

one_seq_subset <- nova240826_df %>% 
  filter(index %in% top_10_in_RBM10$index) %>% 
  filter(grepl("MUT-sg", condition)) %>% 
  filter(!grepl("9999", offset)) %>%
  mutate(full_status = paste(mode, offset, sep = "_"))

# Plot the count for every mode in this sequence.
ggplot(one_seq_subset, aes(full_status, count, fill = condition)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  ggtitle("ENSG00000080298.16;RFX3;chr9-3228846-3228889-3218296-3225280-3248031-3248185") +
  facet_wrap(~index, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


####### Quickly check the K700E data #####
K700E_df <- fread("~/Dropbox (Harvard University)/02Splicing/latest/K700E_umi_count_merged_to_ref_normalized.csv")
# Normalize every sample so the total count is 50 million.
K700E_df <- K700E_df %>% 
  group_by(sample) %>% 
  mutate(total_count = sum(count)) %>% 
  ungroup() %>% 
  mutate(count = count * 5e7 / total_count) %>% 
  mutate(count = as.integer(count)) %>%
  select(-total_count) %>% 
  ungroup() %>% 
  # Also get fraction for each mode per index.
  group_by(sample, index) %>% 
  mutate(total_count = sum(count)) %>% 
  ungroup() %>% 
  mutate(fraction = count / total_count) 

different_modes <- K700E_df %>% 
  group_by(condition, index, mode, offset) %>% 
  summarise(count = sum(count),
            fraction = mean(fraction)) %>% 
  ungroup()


K562_K700E_modes <- different_modes %>%
  filter(condition == "K562_K700E") %>% 
  # dplyr::rename(count_K562_K700E = count) %>% 
  dplyr::rename(fraction_K562_K700E = fraction) %>%
  select(-condition)

K562_WT_modes <- different_modes %>% 
  filter(condition == "K562_WT") %>%
  # dplyr::rename(count_K562_WT = count) %>%
  dplyr::rename(fraction_K562_WT = fraction) %>%
  select(-condition)

K562_to_control <- merge(K562_WT_modes, K562_K700E_modes, by = c("index", "mode", "offset"), all = T) %>% 
  # Filter for at least 1 is NA in fraction
  filter(is.na(fraction_K562_WT) | is.na(fraction_K562_K700E)) %>%
  filter(!grepl("9999", offset)) %>% 
  arrange(desc(fraction_K562_K700E))

top_10_in_K700E <- K562_to_control %>%
  arrange(desc(fraction_K562_K700E)) %>%
  head(12)

one_seq_subset <- K700E_df %>%
  filter(index %in% top_10_in_K700E$index) %>%
  filter(!grepl("9999", offset)) %>%
  mutate(full_status = paste(mode, offset, sep = "_"))

# Plot the count for every mode in this sequence.
ggplot(one_seq_subset, aes(full_status, fraction, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  ggtitle("K562 - WT vs K700E") +
  facet_wrap(~index, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
