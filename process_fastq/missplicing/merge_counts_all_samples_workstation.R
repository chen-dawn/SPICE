library(tidyverse)
library(vroom)
library(data.table)

# OKAY I SHOULD BE MORE SYSTEMATIC. AIYO.
out_dir <- "U:/processed_data/missplicing_processed_df/V5"
# Create outdir.
dir.create(out_dir, showWarnings = FALSE)
input_dir <- "U:/processed_data/missplicing_test_v5/"
alt_ref_file <- read_tsv("U:/melange/data/guide_library_cleaned/ref_test_alt_ref_dict.tsv")

##### Uncomment to process files again. This takes a while. #####
cellline_filenames <- list.files(path = input_dir, pattern = "fine_grained_idx_formatted.tsv$", full.names = TRUE)

all_files <- vroom(cellline_filenames, id = "filename", delim = "\t")
# fwrite(all_files, file.path(out_dir, "umi_count_all_celltypes.csv"))
all_files_df <- all_files %>%  mutate(sample = str_extract(basename(filename), ".+(?=_umi_dedup)")) %>%
  mutate(condition = str_extract(sample, "^.+(?=-rep\\d)")) %>%
  select(-filename)
fwrite(all_files_df, file.path(out_dir, "umi_count_all_celltypes_formatted.csv"))


#### This is processing for K700E samples only. ####
input_dir <- "U:/processed_data/missplicing_test_v5/K700E/"
cellline_filenames <- list.files(path = input_dir, pattern = "fine_grained_idx_formatted.tsv$", full.names = TRUE)

all_files <- vroom(cellline_filenames, id = "filename", delim = "\t")
# fwrite(all_files, file.path(out_dir, "K700E_umi_count_all_celltypes.csv"))
all_files_df <- all_files %>%  mutate(sample = str_extract(basename(filename), ".+(?=_umi_dedup)")) %>%
  mutate(condition = str_extract(sample, "^.+(?=-H\\d)")) %>%
  select(-filename)
fwrite(all_files_df, file.path(out_dir, "K700E_umi_count_all_celltypes_formatted.csv"))
##### End processing for K700E samples only. ####



replace_exact_match <- function(index_str, alt_ref_df) {
  alt_to_ref <- setNames(alt_ref_df$ref, alt_ref_df$alt)
  return(ifelse(index_str %in% names(alt_to_ref), alt_to_ref[index_str], index_str))
}

all_files_df <- fread(file.path(out_dir, "umi_count_all_celltypes_formatted.csv"))
K700E_df <- fread(file.path(out_dir, "K700E_umi_count_all_celltypes_formatted.csv"))


# Apply the function to the "index" column of K700E_df
K700E_df_to_ref <- K700E_df %>%
  mutate(index = sapply(index, replace_exact_match, alt_ref_file)) %>%
  group_by(sample, condition, index, mode, offset) %>% 
  summarise(count = sum(count)) %>% 
  ungroup()
fwrite(K700E_df_to_ref, file.path(out_dir, "K700E_umi_count_merged_to_ref_normalized.csv"))

# Also apply it to the all_files_df
all_files_df_to_ref <- all_files_df %>%
  mutate(index = sapply(index, replace_exact_match, alt_ref_file)) %>%
  group_by(sample, condition, index, mode, offset) %>% 
  summarise(count = sum(count)) %>% 
  ungroup()
fwrite(all_files_df_to_ref, file.path(out_dir, "umi_count_merged_to_ref_normalized.csv"))
replace_exact_match <- function(index_str, alt_ref_df) {
  alt_to_ref <- setNames(alt_ref_df$ref, alt_ref_df$alt)
  return(ifelse(index_str %in% names(alt_to_ref), alt_to_ref[index_str], index_str))
}

# Apply the function to the "index" column of K700E_df
K700E_df_to_ref <- K700E_df %>%
  mutate(index = sapply(index, replace_exact_match, alt_ref_file)) %>%
  group_by(sample, condition, index, mode, offset) %>% 
  summarise(count = sum(count)) %>% 
  ungroup()
fwrite(K700E_df_to_ref, file.path(out_dir, "K700E_umi_count_merged_to_ref_normalized.csv"))

# Also apply it to the all_files_df
all_files_df_to_ref <- all_files_df %>%
  mutate(index = sapply(index, replace_exact_match, alt_ref_file)) %>%
  group_by(sample, condition, index, mode, offset) %>% 
  summarise(count = sum(count)) %>% 
  ungroup()
fwrite(all_files_df_to_ref, file.path(out_dir, "umi_count_merged_to_ref_normalized.csv"))