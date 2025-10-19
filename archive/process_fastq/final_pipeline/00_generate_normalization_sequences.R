library(tidyverse)
library(vroom)
library(data.table)
library(pheatmap)
library(preprocessCore)

all_files_df_path <- "U:/processed_data/latest/250131_merged_v2/WT_all_samples_raw_counts.csv"

all_files_df <- vroom(all_files_df_path)

# First get all perfect sequences.
perfect_PSI <- all_files_df %>% 
  filter((mode == "INCLUDED" & offset == "0:0:0")|(mode == "SKIPPED" & offset == "0")) %>% 
  group_by(index, sample) %>% 
  mutate(total_count = sum(count_scaled)) %>% 
  ungroup() %>% 
  filter(mode == "INCLUDED") %>% 
  mutate(skipped = total_count - count_scaled) %>%
  mutate(included = count_scaled) %>%
  select(-count, -count_scaled) 

perfect_PSI_pass_filter <- perfect_PSI %>% 
  filter(total_count > 30) %>% 
  mutate(PSI = included/total_count) 

# Get number of pass filter per sample:
num_pass_filter_by_sample <- perfect_PSI_pass_filter %>% 
  group_by(sample) %>% 
  summarise(num_pass = n())

# Plot heatmap of values:
PSI_mat <- perfect_PSI_pass_filter %>% 
  group_by(index, condition) %>% 
  summarise(PSI = mean(PSI)) %>% 
  pivot_wider(names_from = condition, values_from = PSI) %>% 
  column_to_rownames(var = "index") %>% 
  as.matrix()

# Get average of every row. 
seq_means <- rowMeans(PSI_mat)
# Put the top 1% as dataframe.
good_ref_sequences <- perfect_PSI_pass_filter %>% 
  group_by(index) %>% 
  summarise(total_count = sum(total_count), average_PSI = mean(PSI),
            num_cellline_represented = length(unique(condition))) %>% 
  arrange(desc(average_PSI)) %>% 
  filter(!is.na(average_PSI)) %>% 
  filter(num_cellline_represented == 46) %>% 
  filter(average_PSI > 0.95 & total_count > 100000)

# Save this list. 
write_tsv(good_ref_sequences, "U:/processed_data/latest/250131_merged_v2/sequences_used_for_normalization_SUPERSET.tsv")

set.seed(100)
# Subsample like, 100 sequences from this set.
ref_seq_index_subsample <- sample(good_ref_sequences$index, 100)
# Write to file. 
write_tsv(tibble(index = ref_seq_index_subsample), "U:/processed_data/latest/250131_merged_v2/sequences_used_for_normalization.tsv")
