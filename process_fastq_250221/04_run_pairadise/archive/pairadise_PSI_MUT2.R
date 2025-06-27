library(tidyverse)
library(vroom)
library(data.table)
library(future)
library(future.apply)

# if pairadise is not installed, install it
if(!requireNamespace("PAIRADISE", quietly = TRUE)) {
  BiocManager::install("PAIRADISE")
}
# Also for BiocParallel
if(!requireNamespace("BiocParallel", quietly = TRUE)) {
  BiocManager::install("BiocParallel")
}
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
output_filename <- args[3]

print(paste("Celltype 1: ", celltype1))
print(paste("Celltype 2: ", celltype2))
print(paste("Output filename: ", output_filename))
all_sample_reps <- fread("/broad/dawnccle/processed_data/reprocess_250221/count_normalized_chimeric_rate_considering_included/MUT2_all_samples_PSI_count_table.csv")
all_sample_reps <- all_sample_reps %>% 
  mutate(index_offset = paste0(index, "__", offset)) 
print("Finished reading data")


get_paradise_output <- function(df, condition1, condition2) {
  # Filter data frames based on conditions
  condition1_df <- df %>% filter(condition == condition1) # %>% filter(offset == "0:0:0")
  condition2_df <- df %>% filter(condition == condition2) # %>% filter(offset == "0:0:0")
  
  # Create a list of all unique indices.
  unique_indices <- unique(condition1_df$index_offset)
  
  # Function to pad vectors with zeros
  pad_with_zeros <- function(vec, target_len = 2) {
    if (length(vec) < target_len) {
      vec <- c(vec, rep(0, target_len - length(vec)))
    }
    return(vec)
  }
  
  
  # Create the paradise data frame using dplyr operations and future.apply for parallel processing
  paradise_df <- lapply(unique_indices, function(idx) {
    condition1_index_df <- condition1_df %>% filter(index_offset == idx) %>% arrange(sample)
    condition2_index_df <- condition2_df %>% filter(index_offset == idx) %>% arrange(sample)
    
    I1 <- pad_with_zeros(condition1_index_df$included_count)
    S1 <- pad_with_zeros(condition1_index_df$skipped_count)
    I2 <- pad_with_zeros(condition2_index_df$included_count)
    S2 <- pad_with_zeros(condition2_index_df$skipped_count)
    
    return (data.frame(
      ExonID = idx, 
      I1 = paste(I1, collapse = ","), 
      S1 = paste(S1, collapse = ","), 
      I2 = paste(I2, collapse = ","), 
      S2 = paste(S2, collapse = ","), 
      I_len = 1, 
      S_len = 1
    ))
  }) %>% bind_rows()
  
  return(paradise_df)
}

formatted_df <- get_paradise_output(all_sample_reps, celltype1, celltype2)
print("Finished formatting data")
print(output_filename)
write_tsv(formatted_df, file = output_filename)
