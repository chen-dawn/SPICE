library(tidyverse)
library(vroom)
library(data.table)

input_dir <- "U:/processed_data/latest/250131_merged_v3/WT"
out_dir <- "U:/processed_data/latest/250131_merged_v3"
# Create outdir.
dir.create(out_dir, showWarnings = FALSE)


# List all files in input dir. 
input_filenames <- list.files(path = input_dir, pattern = "umi_dedup_normalized.tsv$", full.names = TRUE)
# Read in all files.
all_files_df <- vroom(input_filenames, id = "filename", delim = ",") %>% 
  # Strip _umi_dedup_normalized.tsv
  mutate(sample = str_extract(basename(filename), ".+(?=_umi_dedup_normalized.tsv)")) %>%
  # Extract the condition from sample.
  mutate(condition = str_extract(sample, "^.+(?=-rep\\d)")) 

# Write to outdir.
all_files_df <- all_files_df %>% select(-filename)
fwrite(all_files_df, file.path(out_dir, "WT_all_samples_raw_counts.csv"))
