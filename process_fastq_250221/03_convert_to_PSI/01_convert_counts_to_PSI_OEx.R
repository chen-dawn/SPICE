library(tidyverse)
library(vroom)
library(data.table)
library(pheatmap)
library(preprocessCore)
library(purrr)
library(RColorBrewer)

all_files_df_path <- "U:/processed_data/reprocess_250221/count_normalized_chimeric_rate_considering_included/OEx_all_samples_raw_counts.csv"
all_files_df <- fread(all_files_df_path)

# > head(all_files_df)
# index
# <char>
#   1: ENSG00000000003.15;TSPAN6;chrX-100632484-100632568-100630758-100630866-100633404-100633539
# 2: ENSG00000000003.15;TSPAN6;chrX-100632484-100632568-100630758-100630866-100633404-100633539
# 3: ENSG00000000003.15;TSPAN6;chrX-100632484-100632568-100630758-100630866-100633404-100633539
# 4: ENSG00000000003.15;TSPAN6;chrX-100632484-100632568-100630758-100630866-100633404-100633539
# 5: ENSG00000000003.15;TSPAN6;chrX-100633930-100634029-100632484-100632568-100635177-100635252
# 6: ENSG00000000003.15;TSPAN6;chrX-100633930-100634029-100632484-100632568-100635177-100635252
# mode offset count count_scaled    sample condition
# <char> <char> <int>        <int>    <char>    <char>
#   1:  INCLUDED 0:-1:0     1            1 769P-rep1      769P
# 2:  INCLUDED  0:0:0    29           29 769P-rep1      769P
# 3:   SKIPPED      0     3            0 769P-rep1      769P
# 4: UNSPLICED    113     3            0 769P-rep1      769P
# 5:  INCLUDED  0:0:0    88           88 769P-rep1      769P
# 6:   SKIPPED      0    15            1 769P-rep1      769P

# First generate a set of all included events. 
all_files_df <- all_files_df %>% 
  group_by(sample,condition, index) %>%
  mutate(total_counts = sum(count_scaled),
         frac_of_total = count_scaled/total_counts) %>% 
  ungroup() %>% 
  # filter(total_counts >= 30) %>% 
  group_by(sample,condition, index, mode) %>%
  mutate(
    total_counts_per_mode = sum(count_scaled),
    frac_of_mode = count_scaled/sum(count_scaled)) %>% 
  ungroup() 


##### Now generate a set of all possible spliced events. #####
# First get a list of all possible perfect splicing events given all the barcodes. 
barcode_table <- read_csv("U:/melange/data/guide_library_cleaned/20240605_twist_library_v3_ID_barcode_ROUT_filtered.csv")
all_perfect_events <- data.frame(index = barcode_table$ID, mode = "INCLUDED", offset = "0:0:0") %>% as_tibble()

# Imperfect events.
all_possible_imperfect_events <- all_files_df %>%
  filter(mode == "INCLUDED" & offset != "0:0:0" & total_counts_per_mode >= 30) %>% 
  group_by(index, mode, offset) %>% 
  summarise(n = n(), 
            max_perc = max(frac_of_mode), 
            max_reads = max(count_scaled),
            .groups = "drop" ) %>% 
  separate(offset, into = c("skipped_exon_start", "skipped_exon_end", "downstream_exon_start"), sep = ":", remove = FALSE) %>%
  filter(skipped_exon_start < 150 & skipped_exon_end < 150)

all_possible_imperfect_events_filtered <- all_possible_imperfect_events %>% 
  filter(max_reads >= 5 & max_perc > 0.01) %>% 
  select(index, mode, offset)

# All possible included events. 
all_possible_included_events <- rbind(all_perfect_events, all_possible_imperfect_events_filtered) %>% 
  arrange(index, mode, offset) 

# Write the set of all possible included events to a file.
write_csv(all_possible_included_events, "U:/melange/process_fastq_250221/03_convert_to_PSI/OEx_all_included_events.csv")

# Now also find the set of all possible exon skipping events. This is just the barcode list. 
all_possible_skipped_events <- data.frame(index = barcode_table$ID, mode = "SKIPPED", offset = "0") %>% as_tibble()

# Now also find a set of all possible unspliced events. Probably this doesn't need to use the adjusted reads. 
all_possible_unspliced_events <- all_files_df %>% 
  filter(mode == "UNSPLICED" ) %>% 
  group_by(mode, offset) %>%
  summarise(n = n(), 
            max_perc = max(frac_of_mode, na.rm = T), 
            max_reads = max(count, na.rm = T),
            .groups = "drop" ) %>% 
  mutate(offset = as.integer(offset)) 

####### Now generate a PSI table for each event. #######
# Convert to data.table for faster lookups
all_files_df <- fread("U:/processed_data/reprocess_250221/count_normalized_chimeric_rate_considering_included/OEx_all_samples_raw_counts.csv")
setDT(all_files_df)
setDT(all_possible_included_events)

# Create fast lookup tables for INCLUDED and SKIPPED modes
included_lookup <- all_files_df[mode == "INCLUDED", .(count_scaled), keyby = .(sample, index, offset)]
skipped_lookup <- all_files_df[mode == "SKIPPED", .(count_scaled), keyby = .(sample, index)]

# Extract sample-condition pairs
sample_conditions <- unique(all_files_df[, .(sample, condition)])
setDT(sample_conditions)  # Ensure it's a data.table

####### Duplicate Events for Each Sample #######
# Initialize empty list to store duplicated tables
expanded_events_list <- vector("list", length = nrow(sample_conditions))  
# Loop over each sample and duplicate `all_possible_included_events`
for (i in seq_len(nrow(sample_conditions))) {
  sample_name <- sample_conditions$sample[i]
  condition <- sample_conditions$condition[i]
  
  # Duplicate events for this sample
  tmp_events <- copy(all_possible_included_events)
  tmp_events[, sample := sample_name]
  tmp_events[, condition := condition]
  
  # Store in list
  expanded_events_list[[i]] <- tmp_events
}

# Combine all duplicated tables into a single data.table
expanded_events <- rbindlist(expanded_events_list)

####### PSI Computation #######
# Join precomputed tables for efficient count retrieval
final_psi_table <- expanded_events[
  included_lookup, on = .(sample, index, offset), included_count := count_scaled][
    skipped_lookup, on = .(sample, index), skipped_count := count_scaled]

# Replace NAs with 0 (indicating no reads found for an event)
final_psi_table[, c("included_count", "skipped_count") := 
                  lapply(.SD, function(x) fifelse(is.na(x), 0, x)), 
                .SDcols = c("included_count", "skipped_count")]
# Write to file
fwrite(final_psi_table, "U:/processed_data/reprocess_250221/count_normalized_chimeric_rate_considering_included/OEx_all_samples_PSI_count_table.csv")

# Also filter to sequences that have total reads >= 30. 
final_psi_table_filtered <- final_psi_table %>% 
  mutate(total_reads = included_count + skipped_count) %>%
  filter(total_reads >= 30) 
fwrite(final_psi_table_filtered, "U:/processed_data/reprocess_250221/count_normalized_chimeric_rate_considering_included/OEx_all_samples_PSI_count_table_filtered.csv")
