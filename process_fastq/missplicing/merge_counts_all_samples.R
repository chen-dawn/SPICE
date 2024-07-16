library(tidyverse)
library(vroom)
library(data.table)
library(Biostrings)

reverse_complement <- function(dna_seq) {
  complement <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
  nucleotides <- unlist(strsplit(dna_seq, ""))
  complement_nucleotides <- complement[nucleotides]
  reverse_complement_seq <- paste(rev(complement_nucleotides), collapse = "")
  return(reverse_complement_seq)
}

# OKAY I SHOULD BE MORE SYSTEMATIC. AIYO.

# ##### Uncomment to process files again. This takes a while. #####
# cellline_filenames <- list.files(path = input_dir, pattern = "fine_grained_idx_formatted.tsv$", full.names = TRUE)
# 
# all_files <- vroom(cellline_filenames, id = "filename")
# # fwrite(all_files, file.path(out_dir, "umi_count_all_celltypes.csv"))
# all_files_df <- all_files %>%  mutate(sample = str_extract(basename(filename), ".+(?=_umi_dedup)")) %>%
#   mutate(condition = str_extract(sample, "^.+(?=-rep\\d)")) %>%
#   select(-filename)
# fwrite(all_files_df, file.path(out_dir, "umi_count_all_celltypes_formatted.csv"))
# 
# #### This is processing for K700E samples only. ####
# input_dir <- "/Volumes/broad_dawnccle/processed_data/missplicing_test/K700E/"
# cellline_filenames <- list.files(path = input_dir, pattern = "fine_grained_idx_formatted.tsv$", full.names = TRUE)
# 
# all_files <- vroom(cellline_filenames, id = "filename")
# # fwrite(all_files, file.path(out_dir, "K700E_umi_count_all_celltypes.csv"))
# all_files_df <- all_files %>%  mutate(sample = str_extract(basename(filename), ".+(?=_umi_dedup)")) %>%
#   mutate(condition = str_extract(sample, "^.+(?=-H\\d)")) %>%
#   dplyr::select(-filename)
# fwrite(all_files_df, file.path(out_dir, "K700E_umi_count_all_celltypes_formatted.csv"))
# ##### End processing for K700E samples only. ####

# Read in the data.
out_dir <- "~/Dropbox (Harvard University)/02Splicing/library_47k_missplicing/V2_results/"
# Create outdir.
dir.create(out_dir, showWarnings = FALSE)
input_dir <- "/Volumes/broad_dawnccle/processed_data/missplicing_test/"
all_files_df <- fread("~/Dropbox (Harvard University)/02Splicing/library_47k_missplicing/V2_results/umi_count_all_celltypes_formatted.csv")
K700E_df <- fread("~/Dropbox (Harvard University)/02Splicing/library_47k_missplicing/V2_results/K700E_umi_count_all_celltypes_formatted.csv")
alt_ref_file <- read_tsv("~/melange/data/guide_library_cleaned/ref_test_alt_ref_dict.tsv")
barcodes <- read_csv("/Users/dawnxi/melange/data/guide_library/20230130_twist_library_v3.csv") %>%
  mutate(barcodeRevcomp = sapply(barcode, reverse_complement))

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

############# Also process the V3 files ##################
out_dir <- "~/Dropbox (Harvard University)/02Splicing/library_47k_missplicing/V3_results/"
all_files_df <- fread("~/Dropbox (Harvard University)/02Splicing/library_47k_missplicing/V3_results/umi_count_all_celltypes_formatted.csv")
K700E_df <- fread("~/Dropbox (Harvard University)/02Splicing/library_47k_missplicing/V3_results/K700E_umi_count_all_celltypes_formatted.csv")
alt_ref_file <- read_tsv("~/melange/data/guide_library_cleaned/ref_test_alt_ref_dict.tsv")
barcodes <- read_csv("/Users/dawnxi/melange/data/guide_library/20230130_twist_library_v3.csv") %>%
  mutate(barcodeRevcomp = sapply(barcode, reverse_complement))

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
fwrite(all_files_df_to_ref, file.path(out_dir, "umi_count_merged_to_ref_normalized.csv"), nThread = 12)

############# Also process the V4 files ##################
out_dir <- "~/Dropbox (Harvard University)/02Splicing/library_47k_missplicing/V4_results/"
all_files_df <- fread("~/Dropbox (Harvard University)/02Splicing/library_47k_missplicing/V4_results/umi_count_all_celltypes_formatted.csv")
K700E_df <- fread("~/Dropbox (Harvard University)/02Splicing/library_47k_missplicing/V4_results/K700E_umi_count_all_celltypes_formatted.csv")
alt_ref_file <- read_tsv("~/melange/data/guide_library_cleaned/ref_test_alt_ref_dict.tsv")
barcodes <- read_csv("/Users/dawnxi/melange/data/guide_library/20230130_twist_library_v3.csv") %>%
  mutate(barcodeRevcomp = sapply(barcode, reverse_complement))

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
fwrite(all_files_df_to_ref, file.path(out_dir, "umi_count_merged_to_ref_normalized.csv"), nThread = 12)

############# Also process the V5 files ##################
out_dir <- "~/Dropbox (Harvard University)/02Splicing/library_47k_missplicing/V5_results/"
all_files_df <- fread("~/Dropbox (Harvard University)/02Splicing/library_47k_missplicing/V5_results/umi_count_all_celltypes_formatted.csv")
K700E_df <- fread("~/Dropbox (Harvard University)/02Splicing/library_47k_missplicing/V5_results/K700E_umi_count_all_celltypes_formatted.csv")
alt_ref_file <- read_tsv("~/melange/data/guide_library_cleaned/ref_test_alt_ref_dict.tsv")
barcodes <- read_csv("/Users/dawnxi/melange/data/guide_library/20230130_twist_library_v3.csv") %>%
  mutate(barcodeRevcomp = sapply(barcode, reverse_complement))

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

fwrite(all_files_df_to_ref, file.path(out_dir, "umi_count_merged_to_ref_normalized.csv"), nThread = 12)

  
