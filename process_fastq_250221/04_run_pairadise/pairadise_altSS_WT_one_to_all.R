library(tidyverse)
library(vroom)
library(data.table)
library(future)
library(future.apply)

# if pairadise is not installed, install it
if (!requireNamespace("PAIRADISE", quietly = TRUE)) {
  BiocManager::install("PAIRADISE")
}
# Also for BiocParallel
if (!requireNamespace("BiocParallel", quietly = TRUE)) {
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
all_sample_reps <- fread("/broad/dawnccle/processed_data/reprocess_250221/count_normalized_v4_merged/WT_all_samples_raw_counts.csv")
all_sample_reps <- all_sample_reps %>% 
  filter(mode == "INCLUDED") %>% 
  mutate(index_offset = paste0(index, "__", offset)) 


all_event_pairs <- fread("/broad/dawnccle/processed_data/reprocess_250221/count_normalized_v4_merged/WT_included_all_unique_event_pairs.csv") %>% 
  select(-index) %>% 
  mutate(index_offset_major_minor = paste0(index_offset_major, "___", index_offset_minor))
print("Finished reading data")

# Create a function to process samples to avoid code duplication
process_samples <- function(df, all_event_pairs, samples) {
  # Pre-filter the data frames
  major_indices <- all_event_pairs$index_offset_major
  minor_indices <- all_event_pairs$index_offset_minor
  valid_major_minor_pairs <- all_event_pairs$index_offset_major_minor
  
  # Filter initial dataframe to reduce size before operations
  sample_df <- df %>% 
    filter(sample %in% samples) %>%
    filter(index_offset %in% c(major_indices, minor_indices)) %>% 
    select(-count) 
  
  # Create smaller combinations tables for major and minor separately
  # This is more memory efficient than one large expand.grid
  major_combinations <- expand.grid(
    index_offset = unique(intersect(sample_df$index_offset, major_indices)),
    sample = samples,
    stringsAsFactors = FALSE
  )

  # Process major indices
  sample_df_major <- major_combinations %>%
    left_join(sample_df, by = c("index_offset", "sample")) %>%
    mutate(count_scaled = coalesce(count_scaled, 0)) %>% 
    separate(index_offset, into = c("index", "offset"), sep = "__", remove = FALSE)
  
  # Clean up to free memory
  rm(major_combinations)
  gc()
  
    
  minor_combinations <- expand.grid(
    index_offset = unique(intersect(sample_df$index_offset, minor_indices)),
    sample = samples,
    stringsAsFactors = FALSE
  )

  # Process minor indices
  sample_df_minor <- minor_combinations %>%
    left_join(sample_df, by = c("index_offset", "sample")) %>%
    select(-condition) %>%
    mutate(count_scaled = coalesce(count_scaled, 0)) %>% 
    separate(index_offset, into = c("index", "offset"), sep = "__", remove = FALSE)
  
  # Clean up to free memory
  rm(minor_combinations, sample_df)
  gc()
  
  # Join and filter
  sample_df_with_counts <- sample_df_major %>% 
    inner_join(sample_df_minor, 
              by = c("index" = "index", "sample" = "sample"),
              relationship = "many-to-many") %>% 
    mutate(index_offset_major_minor = paste0(index_offset.x, "___", index_offset.y)) %>% 
    filter(index_offset_major_minor %in% valid_major_minor_pairs)
  
  return(sample_df_with_counts)
}

get_paradise_output <- function(df, all_event_pairs, condition1, condition2) {
  # Get unique samples for each condition
  unique_samples_condition1 <- df %>% 
    filter(condition == condition1) %>% 
    pull(sample) %>% 
    unique()
  
  if (condition2 == "all") {
    unique_samples_condition2 <- df %>% 
      filter(condition != condition1) %>% 
      pull(sample) %>% 
      unique()
  } else {
    unique_samples_condition2 <- df %>% 
      filter(condition == condition2) %>% 
      pull(sample) %>% 
      unique()
  }
  
  # Process each condition using the function
  condition1_df <- process_samples(df, all_event_pairs, unique_samples_condition1)
  condition2_df <- process_samples(df, all_event_pairs, unique_samples_condition2)
  
  condition1_df_summary <- condition1_df %>% 
    mutate(altSS_1 = count_scaled.x/(count_scaled.x + count_scaled.y)) %>% 
    group_by(index_offset_major_minor) %>% 
    summarise(altSS_1 = mean(altSS_1, na.rm = TRUE)) %>% 
    select(index_offset_major_minor, altSS_1)
  
  condition2_df_summary <- condition2_df %>% 
    mutate(altSS_2 = (count_scaled.x)/(count_scaled.x + count_scaled.y)) %>% 
    group_by(index_offset_major_minor) %>% 
    summarise(altSS_2 = mean(altSS_2, na.rm = TRUE)) %>% 
    select(index_offset_major_minor, altSS_2)
  
  condition1_2_merged <- merge(condition1_df_summary, condition2_df_summary, by = c("index_offset_major_minor")) %>% 
    mutate(altSS_diff = abs(altSS_1 - altSS_2)) %>%
    filter(altSS_diff >= 0.1)
  
  
  # Create a list of all unique indices.
  unique_indices <- unique(condition1_2_merged$index_offset_major_minor)
  
  # Function to pad vectors with zeros
  pad_with_zeros <- function(vec, target_len = 2) {
    if (length(vec) < target_len) {
      vec <- c(vec, rep(0, target_len - length(vec)))
    }
    return(vec)
  }
  
  
  # Create the paradise data frame using dplyr operations and future.apply for parallel processing
  paradise_df <- lapply(unique_indices, function(idx) {
    condition1_index_df <- condition1_df %>% filter(index_offset_major_minor == idx) %>% arrange(sample)
    condition2_index_df <- condition2_df %>% filter(index_offset_major_minor == idx) %>% arrange(sample)
    
    I1 <- pad_with_zeros(condition1_index_df$count_scaled.x)
    S1 <- pad_with_zeros(condition1_index_df$count_scaled.y)
    I2 <- pad_with_zeros(condition2_index_df$count_scaled.x)
    S2 <- pad_with_zeros(condition2_index_df$count_scaled.y)
    
    return(data.frame(
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

formatted_df <- get_paradise_output(all_sample_reps, all_event_pairs, celltype1, celltype2)
print("Finished formatting data")
print(output_filename)
write_tsv(formatted_df, file = output_filename)
