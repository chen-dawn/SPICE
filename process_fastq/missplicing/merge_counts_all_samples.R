library(tidyverse)
library(vroom)
library(data.table)

# OKAY I SHOULD BE MORE SYSTEMATIC. AIYO.
out_dir <- "~/Dropbox (Harvard University)/02Splicing/library_47k_missplicing/"
# Create outdir.
dir.create(out_dir, showWarnings = FALSE)
input_dir <- "/Volumes/broad_dawnccle/processed_data/missplicing_test/"

##### Uncomment to process files again. This takes a while. #####
cellline_filenames <- list.files(path = input_dir, pattern = "fine_grained_idx_formatted.tsv$", full.names = TRUE)

all_files <- vroom(cellline_filenames, id = "filename")
fwrite(all_files, file.path(out_dir, "umi_count_all_celltypes.csv"))
all_files_df <- all_files %>%  mutate(sample = str_extract(basename(filename), ".+(?=_umi_dedup)")) %>%
  mutate(condition = str_extract(sample, "^.+(?=-rep\\d)")) %>%
  select(-filename)
fwrite(all_files_df, file.path(out_dir, "umi_count_all_celltypes_formatted.csv"))

#### This is processing for K700E samples only. ####
input_dir <- "/Volumes/broad_dawnccle/processed_data/missplicing_test_v3/K700E/"
cellline_filenames <- list.files(path = input_dir, pattern = "fine_grained_idx_formatted.tsv$", full.names = TRUE)

all_files <- vroom(cellline_filenames, id = "filename")
fwrite(all_files, file.path(out_dir, "K700E_umi_count_all_celltypes.csv"))
all_files_df <- all_files %>%  mutate(sample = str_extract(basename(filename), ".+(?=_umi_dedup)")) %>%
  mutate(condition = str_extract(sample, "^.+(?=-H\\d)")) %>%
  dplyr::select(-filename)
fwrite(all_files_df, file.path(out_dir, "K7000E_umi_count_all_celltypes_formatted.csv"))
##### End processing for K700E samples only. ####

# Read in the data.
all_files_df <- fread("~/Dropbox (Harvard University)/02Splicing/library_47k_missplicing/V3umi_count_all_celltypes_formatted.csv")
K700E_df <- fread("~/Dropbox (Harvard University)/02Splicing/library_47k_missplicing/K7000E_umi_count_all_celltypes_formatted.csv")

##### Let's just check out the K700E samples first #####
unspliced <- K700E_df %>% filter(mode == "UNSPLICED") %>% mutate(offset = as.integer(offset))
ggplot(unspliced, aes(offset, count)) + geom_point() + facet_wrap(~sample, scales = "free_y")

# For every index take the average of the unsplice counts for condition.
unspliced_avg <- unspliced %>% 
  group_by(index, condition, offset) %>% 
  summarise(avg_count = mean(count)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = condition, values_from = avg_count, values_fill = 0) 

ggplot(unspliced_avg, aes(K562_WT, K562_K700E)) + geom_point()



# Get included reads only.
included_reads <- K700E_df %>% filter(mode == "INCLUDED")

# For every offset calculate the fraction of included reads.
included_reads <- K700E_df %>% filter(mode == "INCLUDED") %>% 
  group_by(sample, index) %>% 
  mutate(total_sum = sum(count)) %>% 
  filter(total_sum > 30) %>% 
  ungroup() %>% 
  mutate(fraction = count / total_sum) %>% 
  # Average by condition.
  group_by(condition, index, offset) %>%
  summarise(fraction = mean(fraction)) 

# Pivot wider by count. 
included_reads_wide <- included_reads %>% 
  pivot_wider(names_from = condition, values_from = fraction, values_fill = NA) %>% 
  mutate(diff = K562_K700E - K562_WT)

high_in_K700E <- included_reads_wide %>% 
  filter(offset == "0:0:0") %>% 
  filter(diff < -0.3 & K562_WT > 0.9) %>% 
  arrange(diff)

# Get the ones that are high in K700E and check the entire counts.
shortlisted_elements <- included_reads_wide %>% 
  filter(index %in% high_in_K700E$index) %>%
  pivot_longer(cols = c(K562_WT, K562_K700E), names_to = "condition", values_to = "fraction") 

ggplot(shortlisted_elements, aes(offset, fraction, fill = condition)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  facet_wrap(~index, scales = "free_x") + 
  ggtitle("High in K700E, low in WT") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

test <-shortlisted_elements %>% group_by(index, condition) %>% summarise(total = sum(fraction, na.rm = T))



# K700E old and new
