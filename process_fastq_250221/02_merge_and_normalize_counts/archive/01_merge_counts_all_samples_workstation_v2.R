library(tidyverse)
library(vroom)
library(data.table)

####### Systematic processing. ########
out_dir <- "U:/processed_data/reprocess_250221/count_normalized_chimeric_rate/"
# Create outdir.
dir.create(out_dir, showWarnings = FALSE)


# First we look at all the samples that are available. 
input_dir1 <- "U:/processed_data/reprocess_250221/nova230516/"
input_dir2 <- "U:/processed_data/reprocess_250221/nova240826/"
input_dir3 <- "U:/processed_data/reprocess_250221/nova241106/"
# input_dir4 <- "U:/processed_data/latest/raw_K700E/"

# Get all full filenames in all the input dirs. 
cellline_filenames1 <- list.files(path = input_dir1, pattern = "umi_dedup_fine_grained_idx.csv$", full.names = TRUE)
cellline_filenames2 <- list.files(path = input_dir2, pattern = "umi_dedup_fine_grained_idx.csv$", full.names = TRUE)
cellline_filenames3 <- list.files(path = input_dir3, pattern = "umi_dedup_fine_grained_idx.csv$", full.names = TRUE)
# cellline_filenames4 <- list.files(path = input_dir4, pattern = "fine_grained_idx_formatted.tsv$", full.names = TRUE)

# Make a df of these filenames.
cellline_paths <- data.frame(filename = c(cellline_filenames1, cellline_filenames2, cellline_filenames3))
cellline_paths <- cellline_paths %>% 
  mutate(basename = basename(filename)) %>%
  mutate(sample = str_extract(basename(filename), ".+(?=_umi_dedup)")) %>%
  mutate(condition = str_extract(sample, "^.+(?=-rep\\d)")) %>% 
  # mutate condition to strip the following: _100tfx, _150tfx, _1ugNuc, _2ugNuc.
  mutate(condition = str_replace(condition, "_\\d+tfx", "")) %>%
  mutate(condition = str_replace(condition, "_\\d+ugNuc", "")) %>% 
  # Filter out these conditions: "OSRC2", "Kelly_old", "SKNAS_Nuc"
  filter(!condition %in% c("OSRC2", "Kelly_old", "SKNAS_Nuc")) %>% 
  # Filter out filename that contain raw_47celltype/A172 or raw_47celltype/KMRC1 or 47celltype/K562_1ugNuc
  filter(!str_detect(filename, "raw_47celltype/A172|raw_47celltype/KMRC1|47celltype/K562_1ugNuc")) %>%
  # Strip SKNAS_tfx of the _tfx.
  mutate(condition = str_replace(condition, "_tfx", "")) %>% 
  # ALso extract _Nuc.
  mutate(condition = str_replace(condition, "_Nuc", "")) %>%
  # Change the DBTR05MG to DBTRG05MG
  mutate(condition = str_replace(condition, "DBTR05MG", "DBTRG05MG")) %>%
  # Change MEWO to MeWo
  mutate(condition = str_replace(condition, "MEWO", "MeWo")) %>%
  # Change JHOM to JHOM1
  mutate(condition = str_replace(condition, "JHOM", "JHOM1")) %>%
  mutate(rep_old = str_extract(basename(filename), "rep\\d")) %>%
  # Create rep new, which is rep1 = rep1, rep2 = rep2, rep3 = rep3, rep4 = rep1, rep5 = rep2, rep6 = rep3.
  mutate(rep_new = case_when(
    rep_old == "rep1" ~ "rep1",
    rep_old == "rep2" ~ "rep2",
    rep_old == "rep3" ~ "rep3",
    str_detect(basename(filename), "rep4") ~ "rep1",
    str_detect(basename(filename), "rep5") ~ "rep2",
    str_detect(basename(filename), "rep6") ~ "rep3"
  )) %>% 
  mutate(sample_new = paste0(condition, "-", rep_new))


# Extract all metrics and append to cellline_paths
metrics_list <- lapply(1:nrow(cellline_paths), function(i) {
  filename_full <- cellline_paths$filename[i]
  sample_name <- cellline_paths$sample[i]
  filename_dir <- dirname(filename_full)
  stats_log_path <- file.path(filename_dir, paste0(sample_name, "_stats_log_fine_grained_idx.txt"))
  
  if (file.exists(stats_log_path)) {
    stats_log <- read_csv(stats_log_path, col_names = FALSE)
    colnames(stats_log) <- c("metric", "count")
    
    # Convert metrics to a named list
    metric_values <- as.list(setNames(stats_log$count, stats_log$metric))
    return(metric_values)
  } else {
    return(NULL)
  }
})

# Combine metrics with cellline_paths
metrics_df <- bind_rows(metrics_list) %>%
  mutate_all(as.character)  # Ensure consistent data type

cellline_paths <- bind_cols(cellline_paths, metrics_df)

# Filter out samples that have < 1M aligned reads.
cellline_paths_filtered <- cellline_paths %>% 
  mutate(total_aligned_reads = as.integer(total_aligned_reads)) %>% 
  filter(total_aligned_reads >= 1e6) 

# Write this metadata to file.
write_csv(cellline_paths_filtered, "U:/melange/process_fastq_250221/02_merge_and_normalize_counts/cellline_sample_metadata.csv")

# Get all unique sample_new
unique_samples <- cellline_paths_filtered %>% group_by(sample_new) %>% summarise(n=n()) %>% ungroup() 
unique_sample_names <- unique(cellline_paths_filtered$sample_new)

for (sample_tmp in unique_sample_names){
  # Get all the filepaths with that sample name.
  sample_filepaths <- cellline_paths_filtered %>% filter(sample_new == sample_tmp) %>% pull(filename)
  tmp_out <- data.frame()
  for (filepath in sample_filepaths) {
    base::print(paste("Processing", filepath))
    # Read in the tsv file.
    tmp <- vroom(filepath, id = "filename", delim = ",")
    
    # Get chimeric rate from the metadata table. 
    chimeric_rate <- cellline_paths_filtered %>% 
      filter(filename == filepath) %>% 
      pull(perc_chimera_reads)
    base::print(paste("chimeric_rate:", chimeric_rate))
    
    # Separate the "index" column into id, mode, offset, insert_size based on __ separator.
    tmp_to_ref <- tmp %>% 
      separate(index, into = c("index", "mode", "offset", "insert_size"), sep = "__", remove = FALSE)
    
    # Split into included (no need adjustment), and everything else (need adjustment).
    # We just directly adjust the counts for non-included as count_scaled = count * (1 - chimeric_rate).
    tmp_included_sequences_only <- tmp_to_ref %>% filter(mode == "INCLUDED") %>% 
      mutate(count_scaled = count)
    tmp_everything_else <- tmp_to_ref %>% filter(mode != "INCLUDED") %>% 
      mutate(count_scaled = count * (1 - as.numeric(chimeric_rate))) %>% 
      mutate(count_scaled = as.integer(count_scaled)) 
    
    tmp_final <- bind_rows(tmp_included_sequences_only, tmp_everything_else)
    tmp_out <- bind_rows(tmp_out, tmp_final)
  }
  
  # Group the tmp_out.
  tmp_grouped <- tmp_out %>% group_by(index, mode, offset) %>%
    summarise(count = sum(count), 
              count_scaled = sum(count_scaled)) %>% 
    ungroup()
  
  # Write to outdir. 
  base::print(paste("Writing to", file.path(out_dir, paste0(sample_tmp, "_umi_dedup_normalized.tsv"))))
  fwrite(tmp_grouped, file.path(out_dir, paste0(sample_tmp, "_umi_dedup_normalized.tsv")))
}
