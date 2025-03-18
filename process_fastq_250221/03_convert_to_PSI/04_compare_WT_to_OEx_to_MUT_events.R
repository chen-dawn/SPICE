library(tidyverse)
library(pheatmap)
library(data.table)

WT_included_events <- read_csv("U:/melange/process_fastq_250221/03_convert_to_PSI/WT_all_included_events.csv")
MUT_included_events <- read_csv("U:/melange/process_fastq_250221/03_convert_to_PSI/MUT_all_included_events.csv")
OEx_included_events <- read_csv("U:/melange/process_fastq_250221/03_convert_to_PSI/OEx_all_included_events.csv")

# > WT_included_events
# # A tibble: 96,278 × 3
# index                                                    mode  offset
# <chr>                                                    <chr> <chr> 
#   1 ENSG00000000003.15;TSPAN6;chrX-100632484-100632568-1006… INCL… 0:-1:0
# 2 ENSG00000000003.15;TSPAN6;chrX-100632484-100632568-1006… INCL… 0:-69…
# 3 ENSG00000000003.15;TSPAN6;chrX-100632484-100632568-1006… INCL… 0:0:0 
# 4 ENSG00000000003.15;TSPAN6;chrX-100632484-100632568-1006… INCL… 15:0:0
# 5 ENSG00000000003.15;TSPAN6;chrX-100633930-100634029-1006… INCL… -1:0:0
# 6 ENSG00000000003.15;TSPAN6;chrX-100633930-100634029-1006… INCL… 0:-1:0
# 7 ENSG00000000003.15;TSPAN6;chrX-100633930-100634029-1006… INCL… 0:0:0 
# 8 ENSG00000000003.15;TSPAN6;chrX-100633930-100634029-1006… INCL… 11:0:0
# 9 ENSG00000000003.15;TSPAN6;chrX-100635177-100635252-1006… INCL… 0:0:0 
# 10 ENSG00000000419.14;DPM1;chr20-50945736-50945762-5094203… INCL… 0:0:0 

# Get events that are in mut but not in WT.
MUT_only <- anti_join(MUT_included_events, WT_included_events, by = c("index", "mode", "offset")) %>% 
  mutate(index_offset = paste(index, offset, sep = "__"))

# Get events in OEx only.
OEx_only <- anti_join(OEx_included_events, WT_included_events, by = c("index", "mode", "offset")) %>% 
  mutate(index_offset = paste(index, offset, sep = "__"))

final_psi_table_filtered <- fread("U:/processed_data/reprocess_250221/count_normalized_chimeric_rate_considering_included/OEx_all_samples_PSI_count_table.csv")
final_psi_table_filtered <- final_psi_table_filtered %>% 
  mutate(index_offset = paste(index, offset, sep = "__")) %>% 
  # filter(offset == "0:0:0") %>% 
  select(-mode) %>%
  mutate(total_count = included_count + skipped_count) %>%
  filter(total_count >= 20) %>%
  mutate(PSI = included_count/(included_count + skipped_count)) %>% 
  # Split the offset.
  separate(offset, into = c("skipped_exon_start", "skipped_exon_end", "downstream_exon_start"), sep = ":", remove = FALSE) %>% 
  mutate(skipped_exon_start = as.integer(skipped_exon_start),
         skipped_exon_end = as.integer(skipped_exon_end)) %>%
  filter(skipped_exon_start < 150 & skipped_exon_end < 150) %>%
  filter(abs(skipped_exon_start) != 1 & abs(skipped_exon_end) != 1) 


final_psi_table_pivot <- final_psi_table_filtered %>%
  select(sample, index_offset, PSI) %>%
  pivot_wider(names_from = c(sample), values_from = PSI) 

shortlist_MUT_only <- final_psi_table_pivot %>% 
  filter(index_offset %in% OEx_only$index_offset) 

shortlist_MUT_only_mat <- as.matrix(shortlist_MUT_only %>% select(-index_offset))
rownames(shortlist_MUT_only_mat) <- shortlist_MUT_only$index_offset
# Filter to rows with < 20% missing data. and columns with < 20% missing data.
shortlist_MUT_only_mat <- shortlist_MUT_only_mat[rowMeans(is.na(shortlist_MUT_only_mat)) < 0.2, colMeans(is.na(shortlist_MUT_only_mat)) < 0.2]


# Plot pheatmap
heatmap_colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)  # Better Blue-White-Red

pheatmap(shortlist_MUT_only_mat, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         color = heatmap_colors, 
         breaks = seq(0, 1, length.out = 100), 
         main = "OEx Only Events")

final_psi_table_filtered <- final_psi_table_filtered %>% 
  group_by(sample,condition, index) %>% 
  mutate(total_included = sum(included_count), num_included = n()) %>% 
  mutate(included_perc = included_count/total_included) %>% 
  ungroup() %>% 
  mutate(skipped_perc = skipped_count/(total_included + skipped_count))

# Filter for num_included > 4
high_num_included_index <- final_psi_table_filtered %>% 
  filter(num_included > 3) %>% pull(index) %>% unique()

high_num_included <- final_psi_table_filtered %>% 
  # filter(index %in% high_num_included_index) %>%
  mutate(has_dox = ifelse(grepl("no_dox", condition), "no_dox", "dox")) %>% 
  mutate(rbp = str_extract(condition, "rbp\\d+")) %>% 
  mutate(rep = str_extract(sample, "rep\\d+")) %>%
  mutate(sample = paste0(rbp, "_", has_dox, "_", rep)) 

delta_psi <- high_num_included %>% 
  group_by(rbp, has_dox, index, offset) %>% 
  summarise(mean_included_perc = mean(included_perc)) %>% 
  pivot_wider(names_from = c(has_dox), values_from = mean_included_perc) %>%
  mutate(delta_included_perc = dox - no_dox) %>%
  filter(abs(delta_included_perc) > 0.2) %>% 
  filter(rbp == "rbp3")
  # filter(rbp != "rbp7" & rbp != "rbp8")

# Pivot data for heatmap
high_num_included_pivot <- high_num_included %>% 
  filter(index %in% delta_psi$index) %>%
  select(sample, index_offset, included_perc) %>% 
  pivot_wider(names_from = c(sample), values_from = included_perc)

# Convert to matrix
high_num_included_mat <- as.matrix(high_num_included_pivot %>% select(-index_offset))
rownames(high_num_included_mat) <- high_num_included_pivot$index_offset

# Filter for missing data threshold
high_num_included_mat <- high_num_included_mat[
  rowMeans(is.na(high_num_included_mat)) < 0.2, 
  colMeans(is.na(high_num_included_mat)) < 0.2
]

# Arrange columns by RBP groups
col_order <- order(colnames(high_num_included_mat))
high_num_included_mat <- high_num_included_mat[, col_order]

# Create annotation dataframe for column groups (RBP groups)
rbp_groups <- data.frame(RBP = str_extract(colnames(high_num_included_mat), "rbp\\d+"))
rownames(rbp_groups) <- colnames(high_num_included_mat)

# Define colors for RBP groups (optional)
rbp_group_colors <- setNames(
  rainbow(length(unique(rbp_groups$RBP))), 
  unique(rbp_groups$RBP)
)

# Plot heatmap with vertical gaps between RBP groups
pheatmap(
  high_num_included_mat, 
  cluster_rows = T, 
  cluster_cols = FALSE, 
  color = heatmap_colors, 
  breaks = seq(0, 1, length.out = 100), 
  main = "High Num Included Events",
  show_rownames = T,
  annotation_col = rbp_groups,  # Add vertical gaps between RBP groups
  gaps_col = which(diff(as.numeric(as.factor(rbp_groups$RBP))) != 0)  # Add gaps between RBP groups
)

