library(tidyverse)
library(vroom)
library(data.table)
library(pheatmap)
library(PAIRADISE)

reverse_complement <- function(dna_seq) {
  complement <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
  nucleotides <- unlist(strsplit(dna_seq, ""))
  complement_nucleotides <- complement[nucleotides]
  reverse_complement_seq <- paste(rev(complement_nucleotides), collapse = "")
  return(reverse_complement_seq)
}

args <- commandArgs(trailingOnly = TRUE)
celltype1 <- args[1]
celltype2 <- args[2]
out_dir <- args[3]

dir.create(out_dir, showWarnings = FALSE)
all_files_df <- fread("~/Dropbox (Harvard University)/02Splicing/library_47k_missplicing/V5_results/umi_count_merged_to_ref_normalized.csv")

###### Look at sample replicates ######
all_sample_reps <- all_files_df %>% 
  filter(mode == "INCLUDED") %>% 
  filter(!condition %in% c("OSRC2", "Kelly_old", "SKNAS_Nuc", "A172")) %>% 
  separate(offset, c("offset_mid_start", "offset_mid_end", "offset_down_start"), sep = ":", convert = T, remove = F) %>% 
  filter(offset_down_start == 0) %>%
  group_by(sample, index) %>% 
  mutate(total_sum = sum(count)) %>% 
  ungroup() %>%
  mutate(other_splice = total_sum - count) 

unique_conditions <- unique(all_sample_reps$condition)

get_paradise_output <- function(df, condition1, condition2) {
  # Filter data frames based on conditions
  condition1_df <- df %>% filter(condition == condition1)
  condition2_df <- df %>% filter(condition == condition2)
  
  # Create a list of the first 10 unique indices
  unique_indices <- unique(condition1_df$index)[1:10]
  
  # Function to pad vectors with zeros
  pad_with_zeros <- function(vec, length = 3) {
    if (length(vec) < length) {
      vec <- c(vec, rep(0, length - length(vec)))
    }
    return(vec[1:length])
  }
  
  # Create the paradise data frame using dplyr operations
  paradise_df <- unique_indices %>%
    lapply(function(idx) {
      condition1_index_df <- condition1_df %>% filter(index == idx) %>% arrange(sample)
      condition2_index_df <- condition2_df %>% filter(index == idx) %>% arrange(sample)
      
      I1 <- pad_with_zeros(condition1_index_df$count)
      S1 <- pad_with_zeros(condition1_index_df$other_splice)
      I2 <- pad_with_zeros(condition2_index_df$count)
      S2 <- pad_with_zeros(condition2_index_df$other_splice)
      
      data.frame(
        ExonID = idx, 
        I1 = paste(I1, collapse = ","), 
        S1 = paste(S1, collapse = ","), 
        I2 = paste(I2, collapse = ","), 
        S2 = paste(S2, collapse = ","), 
        I_len = 1, 
        S_len = 1
      )
    }) %>%
    bind_rows()
  
  return(paradise_df)
}

formatted_df <- get_paradise_output(all_sample_reps, celltype1, celltype2)
pdat <- PDseDataSetFromMat(formatted_df)
pairadise_output <- pairadise(pdat, numCluster = 8)
res <- results(pairadise_output, p.adj = "BH", sig.level = 1)
r1 <- data.frame(res)
# set rowname as index
r1$index <- rownames(r1)
# Merge with test df.
merged_df <- merge(formatted_df, r1, by.x = "ExonID", by.y = "index")

write_tsv(merged_df, file = paste0(out_dir, "/", celltype1, "_", celltype2, "_paradise_output.tsv"))
