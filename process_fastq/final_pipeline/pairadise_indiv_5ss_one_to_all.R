library(tidyverse)
library(vroom)
library(data.table)

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

all_sample_reps <- fread("/broad/dawnccle/processed_data/latest/all_sample_reps_5ss.csv")
print("Finished reading data")

get_paradise_output <- function(df, condition1, condition2) {
  # Filter data frames based on conditions
  condition1_df <- df %>% filter(grepl(condition1, condition)) %>% filter(offset == "0:0:0")
  condition2_df <- df %>% filter(!grepl(condition1, condition)) %>% filter(offset == "0:0:0")
  
  # Create a list of the first 10 unique indices
  unique_indices <- unique(condition1_df$index)
  # Function to pad vectors with zeros
  pad_with_zeros <- function(vec, target_len = 3) {
    if (length(vec) < target_len) {
      vec <- c(vec, rep(0, target_len - length(vec)))
    }
    return(vec)
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
print("Finished formatting data")
write_tsv(formatted_df, file = output_filename)

# pdat <- PDseDataSetFromMat(formatted_df)
# pairadise_output <- pairadise(pdat, numCluster = 8)
# res <- results(pairadise_output, p.adj = "BH", sig.level = 1)
# r1 <- data.frame(res)
# # set rowname as index
# r1$index <- rownames(r1)
# # Merge with test df.
# merged_df <- merge(formatted_df, r1, by.x = "ExonID", by.y = "index")
# print("Writing to file")
# print(paste0(out_dir, "/", celltype1, "_", celltype2, "_pairadise_output.tsv"))
# write_tsv(merged_df, file = paste0(out_dir, "/", celltype1, "_", celltype2, "_pairadise_output.tsv"))
# 
# # Also write the subset of p.adj < 0.05.
# write_tsv(merged_df %>% filter(p.adj < 0.05), file = paste0(out_dir, "/", celltype1, "_", celltype2, "_pairadise_output_padj_0.05.tsv"))
