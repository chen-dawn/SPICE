library(tidyverse)
library(vroom)
library(data.table)

out_dir <- "U:/processed_data/reprocess_250221/count_normalized_chimeric_rate/"
dir.create(out_dir, showWarnings = FALSE)

process_samples <- function(input_dir, sample_type, out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  input_filenames <- list.files(path = input_dir, pattern = "umi_dedup_normalized.tsv$", full.names = TRUE)
  
  # Precompute sample and condition for each file
  file_metadata <- tibble(
    filename = input_filenames,
    sample = basename(input_filenames) %>% str_extract(".+(?=_umi_dedup_normalized.tsv)"),
    condition = str_extract(basename(input_filenames), "^.+(?=-rep\\d)")
  )
  
  # Read data and attach metadata
  all_files_df <- map_dfr(seq_along(file_metadata$filename), function(i) {
    df <- vroom(file_metadata$filename[i], delim = ",")
    df$sample <- file_metadata$sample[i]
    df$condition <- file_metadata$condition[i]
    df
  })
  
  fwrite(all_files_df, file.path(out_dir, paste0(sample_type, "_all_samples_raw_counts.csv")))
}

# Define input and output directories
input_output_mapping <- list(
  list(input_dir = "U:/processed_data/reprocess_250221/count_normalized_chimeric_rate/WT", 
       sample_type = "WT", 
       out_dir = "U:/processed_data/reprocess_250221/count_normalized_chimeric_rate/"),
  
  list(input_dir = "U:/processed_data/reprocess_250221/count_normalized_chimeric_rate/MUT", 
       sample_type = "MUT", 
       out_dir = "U:/processed_data/reprocess_250221/count_normalized_chimeric_rate/"),
  
  list(input_dir = "U:/processed_data/reprocess_250221/count_normalized_chimeric_rate/OEx", 
       sample_type = "OEx", 
       out_dir = "U:/processed_data/reprocess_250221/count_normalized_chimeric_rate/")
)

# Process each sample type
walk(input_output_mapping, ~process_samples(.x$input_dir, .x$sample_type, .x$out_dir))

####################################
# Also process for the other folder.
####################################

out_dir <- "U:/processed_data/reprocess_250221/count_normalized_chimeric_rate_considering_included/"
dir.create(out_dir, showWarnings = FALSE)

process_samples <- function(input_dir, sample_type, out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  input_filenames <- list.files(path = input_dir, pattern = "umi_dedup_normalized.tsv$", full.names = TRUE)
  
  # Precompute sample and condition for each file
  file_metadata <- tibble(
    filename = input_filenames,
    sample = basename(input_filenames) %>% str_extract(".+(?=_umi_dedup_normalized.tsv)"),
    condition = str_extract(basename(input_filenames), "^.+(?=-rep\\d)")
  )
  
  # Read data and attach metadata
  all_files_df <- map_dfr(seq_along(file_metadata$filename), function(i) {
    df <- vroom(file_metadata$filename[i], delim = ",")
    df$sample <- file_metadata$sample[i]
    df$condition <- file_metadata$condition[i]
    df
  })
  
  fwrite(all_files_df, file.path(out_dir, paste0(sample_type, "_all_samples_raw_counts.csv")))
}

# Define input and output directories
input_output_mapping <- list(
  list(input_dir = "U:/processed_data/reprocess_250221/count_normalized_chimeric_rate_considering_included/WT", 
       sample_type = "WT", 
       out_dir = "U:/processed_data/reprocess_250221/count_normalized_chimeric_rate_considering_included/"),
  
  list(input_dir = "U:/processed_data/reprocess_250221/count_normalized_chimeric_rate_considering_included/MUT", 
       sample_type = "MUT", 
       out_dir = "U:/processed_data/reprocess_250221/count_normalized_chimeric_rate_considering_included/"),
  
  list(input_dir = "U:/processed_data/reprocess_250221/count_normalized_chimeric_rate_considering_included/OEx", 
       sample_type = "OEx", 
       out_dir = "U:/processed_data/reprocess_250221/count_normalized_chimeric_rate_considering_included/")
)

# Process each sample type
walk(input_output_mapping, ~process_samples(.x$input_dir, .x$sample_type, .x$out_dir))