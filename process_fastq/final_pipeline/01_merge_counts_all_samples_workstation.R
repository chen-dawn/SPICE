library(tidyverse)
library(vroom)
library(data.table)

# OKAY I SHOULD BE MORE SYSTEMATIC. AIYO.
out_dir <- "U:/processed_data/latest/updated/"
# Create outdir.
dir.create(out_dir, showWarnings = FALSE)
input_dir <- "U:/processed_data/latest/raw_47celltype/"
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
input_dir <- "U:/processed_data/latest/raw_K700E/"
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

####### Make the all_samples file that's cleaned. ########
all_files_df <- fread(file.path(out_dir, "umi_count_merged_to_ref_normalized.csv"))
# Look at sample replicates for 5ss
all_sample_reps_3ss <- all_files_df %>%
  filter(mode == "INCLUDED") %>%
  filter(!condition %in% c("OSRC2", "Kelly_old", "SKNAS_Nuc", "A172")) %>%
  separate(offset, c("offset_mid_start", "offset_mid_end", "offset_down_start"), sep = ":", convert = T, remove = F) %>%
  filter(offset_down_start == 0) %>%
  group_by(sample, index) %>%
  mutate(total_sum = sum(count)) %>%
  ungroup() %>%
  mutate(other_splice = total_sum - count)
fwrite(all_sample_reps_5ss, file = file.path(out_dir, "all_sample_reps_3ss.csv"))

# Look at sample replicates for exon skipping PSI
all_sample_reps_PSI <- all_files_df %>%
  filter((mode == "INCLUDED" & offset == "0:0:0") | (mode == "SKIPPED" & offset == "0")) %>%
  filter(!condition %in% c("OSRC2", "Kelly_old", "SKNAS_Nuc", "A172")) %>%
  group_by(sample, index) %>%
  mutate(total_sum = sum(count)) %>%
  ungroup() %>%
  mutate(skipped = total_sum - count) %>%
  filter(mode == "INCLUDED") %>%
  select(-total_sum)
fwrite(all_sample_reps_PSI, file = file.path(out_dir, "all_sample_reps_PSI.csv"))
