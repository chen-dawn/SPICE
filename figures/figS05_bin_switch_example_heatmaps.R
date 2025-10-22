library(tidyverse)
library(vroom)
library(data.table)
library(future.apply)
library(pheatmap)
library(dplyr)
library(stringr)
library(purrr)
library(broom)
library(patchwork)   
library(tidyr)
library(ggplot2)
library(ggseqlogo)

calculate_upsilon <- function(row) {
  non_na_row <- as.numeric(row[!is.na(row)])
  if (length(non_na_row) < 3) return(NA)
  non_na_row <- non_na_row + 1
  norm_row <- non_na_row / max(non_na_row)
  tau <- sum(1 - norm_row) / (length(non_na_row) - 1)
  return(tau * 2)
}
calculate_reverse_upsilon <- function(row) {
  non_na_row <- as.numeric(row[!is.na(row)])
  if (length(non_na_row) < 3) return(NA)
  non_na_row <- 2 - non_na_row
  norm_row <- non_na_row / max(non_na_row)
  tau <- sum(1 - norm_row) / (length(non_na_row) - 1)
  return(tau * 2)
}
# Define colors
color_palette2 <- c(
  "#4575B4",  # deep blue
  "#85B6D6",  # slightly lighter/more even blue
  "#E2EFF2",  # less stark pastel blue
  "#FFE3B0",  # warmer, slightly less saturated yellow
  "#EF9651",  # softer orange
  "#D83629"   # red
)
color_palette2_custom <- colorRampPalette(color_palette2)(100)
color_palette2_custom_rev <- colorRampPalette(rev(color_palette2))(100)

output_filepath_supp <- "/Volumes/broad_dawnccle/melange/figures_outputs/figS05"


######## Import Data: Bin switching designed sequences  ###########
PSI_file <- read.csv("/Volumes/broad_dawnccle/melange/data/designed_bin_switch_count_table_normalized.csv")
reference_file_bin_switching <- read.csv("/Volumes/broad_dawnccle/melange/data/designed_bin_switch_reference.csv")

raw_df_clean <- PSI_file %>%
  mutate(index_offset = paste0(index, "__", offset_initial)) %>% 
  select(-filename, - index, - offset_initial) %>%
  filter(mode %in% c("INCLUDED", "SKIPPED")) %>% 
  group_by(sample, condition, index_offset, mode, offset, design) %>%
  summarise(count = sum(count)) %>% 
  arrange(sample, condition, index_offset, design)

raw_df_to_psi <- raw_df_clean %>% 
  group_by(sample, condition, index_offset, mode, design) %>%
  summarise(count = sum(count)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = mode, values_from = count, values_fill = 0) %>%
  mutate(total_counts = INCLUDED + SKIPPED) %>%
  filter(total_counts >= 30) %>%
  mutate(PSI = INCLUDED / total_counts)



######## filter to sequences that correctly bin switch ########
merged_df <- raw_df_to_psi %>%
  group_by(index_offset, design, condition) %>%
  summarise(PSI = mean(PSI)) %>%
  mutate(index_offset_design = paste(index_offset, design, sep = "__")) %>%
  left_join(reference_file_bin_switching, by = "index_offset_design") %>%
  select(-barcode, -twist_seq) 


merged_df <- merged_df %>%
  mutate(
    category = if_else(str_detect(design, "_"),
                       str_replace(design, "_[^_]+$", ""),   # everything before last "_"
                       "parent"),
    R1design = if_else(str_detect(design, "_"),
                       str_extract(design, "[^_]+$"),        # everything after last "_"
                       design)
  ) 



# 1) Get parent_PSI per (index_offset, condition)
parent_vals <- merged_df %>%
  filter(category == "parent") %>%
  group_by(index_offset, condition) %>%
  summarise(parent_PSI = mean(PSI, na.rm = TRUE), .groups = "drop")

fallback_vals <- parent_vals %>%
  group_by(index_offset) %>%
  summarise(fallback_parent = mean(parent_PSI, na.rm = TRUE), .groups = "drop")

# Step 3. Join both and coalesce
merged_df <- merged_df %>%
  left_join(parent_vals, by = c("index_offset","condition")) %>%
  left_join(fallback_vals, by = "index_offset") %>%
  mutate(parent_PSI = if_else(is.na(parent_PSI), fallback_parent, parent_PSI)) %>%
  select(-fallback_parent) %>%
  # 2) Split `design` into 4 tokens; keep token #2 (original) and #4 (final)
  mutate(
    toks     = str_split_fixed(category, "_", 4),
    original = if_else(category == "parent", NA_character_, toks[,2]),
    final    = if_else(category == "parent", NA_character_, toks[,4])
  ) %>%
  # 3) Original check against parent_PSI window
  mutate(
    original_pass = case_when(
      category == "parent" ~ TRUE,  # keep parent rows
      original == "0"       & !is.na(parent_PSI) & parent_PSI >= 0.2 & parent_PSI <= 0.8 ~ TRUE,
      original == "10"      & !is.na(parent_PSI) & parent_PSI >  0.8 & parent_PSI <= 1.0 ~ TRUE,
      original == "minus10" & !is.na(parent_PSI) & parent_PSI >= 0.0 & parent_PSI <  0.2 ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  # 4) Final check on this row's PSI window (same bounds as before)
  mutate(
    final_pass = case_when(
      category == "parent" ~ TRUE,
      final == "0"        & PSI >= 0.2 & PSI <= 0.8 ~ TRUE,
      final == "10"       & PSI >  0.8 & PSI <= 1.0 ~ TRUE,
      final == "minus10"  & PSI >= 0.0 & PSI <  0.2 ~ TRUE,
      TRUE ~ FALSE
    ),
    correct = original_pass & final_pass
  ) %>%
  select(-toks) 

merged_df_filtered <- merged_df %>%
  filter(correct == TRUE)

######## Figure S5C: find an example sequences that bin switches properly #########
make_heatmap <- function(index_offset, cats, outfile) {
  df <- merged_df %>%
    filter(index_offset == !!index_offset) %>%
    filter(category %in% !!cats) %>%
    group_by(design, condition) %>%
    summarise(PSI = mean(PSI, na.rm = TRUE), .groups = "drop") %>%
    mutate(trailing_num = as.numeric(str_extract(design, "\\d+$"))) %>%
    mutate(design = factor(design,
                           levels = rev(unique(design[order(trailing_num)])),
                           ordered = TRUE))
  
  p <- ggplot(df, aes(x = condition, y = design, fill = PSI)) +
    geom_tile(color = "white") +
    scale_fill_gradientn(colors = color_palette2_custom, limits = c(0, 1)) +
    theme_minimal() +
    labs(title = "PSI heatmap by design and condition",
         x = "Condition", y = "Design", fill = "PSI")
  
  ggsave(file.path(output_filepath_supp, outfile), p, width = 5, height = 4, dpi = 300)
}

jobs <- tibble::tibble(
  index_offset = c(
    "ENSG00000073711.11;PPP2R3A;chr3-136070477-136070552-136049258-136049361-136078366-136078453__0:0:0",
    "ENSG00000131828.14;PDHA1;chrX-19357651-19357719-19355581-19355757-19358915-19359024__0:0:0",
    "ENSG00000228983.9;SLC47A1P1;chr17-19591764-19591860-19590743-19590788-19594209-19594334__0:0:0",
    "ENSG00000083857.14;FAT1;chr4-186590367-186590403-186588391-186589220-186592692-186592752__0:0:0"
  ),
  cats = list(
    c("vae_10_to_minus10", "parent"),
    c("vae_10_to_minus10", "parent"),
    c("vae_minus10_to_10", "parent"),
    c("vae_minus10_to_10", "parent")
  ),
  outfile = c(
    "figS05_binswitch_example_PPP2R3A.pdf",
    "figS05_binswitch_example_PDHA1.pdf",
    "figS05_binswitch_example_SLC47A1P1.pdf",
    "figS05_binswitch_example_FAT1.pdf"
  )
)

for (i in seq_len(nrow(jobs))) {
  make_heatmap(
    index_offset = jobs$index_offset[i],
    cats        = jobs$cats[[i]],
    outfile     = jobs$outfile[i]
  )
}



