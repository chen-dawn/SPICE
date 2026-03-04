library(tidyverse)
library(vroom)
library(data.table)
library(pheatmap)
library(dplyr)
library(stringr)
library(purrr)
library(broom)
library(patchwork)   
library(tidyr)
library(ggplot2)
library(ggseqlogo)
library(stringdist)


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

output_filepath_main <- "/Volumes/broad_dawnccle/melange/figures_outputs/fig05/"
output_filepath_supp <- "/Volumes/broad_dawnccle/melange/figures_outputs/figS05/"



######## Import Data: KELLY-specific designed sequences ########
KELLY_designed_data_file <- read.csv("/Volumes/broad_dawnccle/melange/data/designed_KELLY_specific_count_table_normalized_PSI.csv")
KELLY_designed_reference <- read.csv("/Volumes/broad_dawnccle/melange/data/designed_KELLY_specific_reference.csv")


KELLY_designed_data <- KELLY_designed_data_file %>% 
  filter(grepl("KELLYspecific", condition)) %>% 
  group_by(condition, index_offset, design) %>% 
  summarise(PSI = mean(PSI)) %>% 
  ungroup() %>% 
  select(condition, index_offset, design, PSI) %>%
  pivot_wider(names_from = condition, values_from = PSI) %>%
  mutate(index_offset_design = paste(index_offset, design, sep = "__")) %>%
  filter(index_offset_design %in% KELLY_designed_reference$ID) 
  

  

######## Import Pairadise data for KELLY designed sequences and add in significance labels to PSI chart ##########
kelly_pairadise <- read_tsv("/Volumes/broad_dawnccle/melange/data/designed_KELLY_specific_rMATS_results.txt", col_types = cols(S1 = col_character(), S2 = col_character()))

kelly_pairadise_PSI <- kelly_pairadise %>%
  mutate(
    I1_sum = sapply(str_split(I1, ","), function(x) sum(as.numeric(x))),
    I2_sum = sapply(str_split(I2, ","), function(x) sum(as.numeric(x))),
    S1_sum = sapply(str_split(S1, ","), function(x) sum(as.numeric(x))),
    S2_sum = sapply(str_split(S2, ","), function(x) sum(as.numeric(x))),
    PSI1 = I1_sum / (I1_sum + S1_sum),
    PSI2 = I2_sum / (I2_sum + S2_sum),
    deltaPSI = PSI1 - PSI2
  )

kelly_pairadise_PSI <- kelly_pairadise_PSI %>%
  mutate(FDR = as.numeric(FDR)) %>%
  mutate(significance = if_else(FDR < 0.05 & abs(deltaPSI) > 0.1, TRUE, FALSE)) 


KELLY_designed_data <- KELLY_designed_data %>%
  left_join(kelly_pairadise_PSI, by = c("index_offset_design" = "ExonID")) %>%
  mutate(
    category = if_else(str_detect(design, "_"),
                       str_replace(design, "_[^_]+$", ""),   # everything before last "_"
                       "parent"),
    R1design = if_else(str_detect(design, "_"),
                       str_extract(design, "[^_]+$"),        # everything after last "_"
                       design)
  ) 

######## Figure S5E: Make volcano plot from KELLY-specific designed sequences #####

# Filter data to remove parent sequences
plot_data <- KELLY_designed_data %>%
  filter(design != "R2design0") 

# Volcano plot
plot <- ggplot(plot_data, aes(x = deltaPSI, y = -log10(FDR))) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed") +
  labs(
    title = "Volcano Plot of KELLY Pairadise PSI",
    x = expression(Delta*"PSI"),
    y = expression(-log[10]*"(FDR)")
  ) +
  ylim(0, 15) +
  xlim(-0.75, 0.75) + 
  theme_minimal()
outfile <- file.path(output_filepath_supp, "figS05_kelly_designed_sequences_volcano.pdf")
ggsave(outfile, plot, width = 3, height = 3, dpi = 300)



######## Figure 5g: Make a heatmap of the top sequences ######

cols_of_interest <- c("KELLYspecific-KELLY", 
                      "KELLYspecific-HEK", 
                      "KELLYspecific-A375", 
                      "KELLYspecific-T47D")

df_sig <- KELLY_designed_data %>%
  filter(significance == TRUE) %>%
  filter(design != "R2design0")

heatmap_matrix <- df_sig[, cols_of_interest] %>%
  as.matrix()

plot <- pheatmap(
  heatmap_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = color_palette2_custom,
  breaks = seq(0, 1, length.out = 101),  # force scale from 0 → 1
  na_col = "grey90",
  main = "KELLY PSI Heatmap (Significant Only)",
  legend = TRUE
)

outfile <- file.path(output_filepath_main, "fig05_KELLY_specific_significant_hits_heatmap.pdf")
ggsave(outfile, plot, width = 5.5, height = 5.0, dpi = 300)


 
######## Figure 5h: Make examples of sequences that are KELLY-specific ##########


make_heatmap <- function(index_offset, cats, outfile) {
  cond_cols <- c("KELLYspecific-A375", "KELLYspecific-HEK",
                 "KELLYspecific-KELLY", "KELLYspecific-T47D")
  cond_cols <- intersect(cond_cols, names(KELLY_designed_data))  
  
  df <- KELLY_designed_data %>%
    filter(index_offset == !!index_offset) %>%
    filter(category %in% !!cats) %>%
    mutate(significance = tidyr::replace_na(significance, FALSE)) %>%
    filter(significance | category == "parent") %>%
    pivot_longer(cols = all_of(cond_cols),
                 names_to = "condition", values_to = "PSI") %>%
    mutate(trailing_num = as.numeric(str_extract(design, "\\d+$"))) %>%
    mutate(design = factor(design,
                           levels = unique(design[order(-trailing_num)]),
                           ordered = TRUE))
  
  p <- ggplot(df, aes(x = condition, y = design, fill = PSI)) +
    geom_tile(color = "white") +
    scale_fill_gradientn(colors = color_palette2_custom, limits = c(0, 1)) +
    theme_minimal() +
    labs(title = "PSI heatmap by design and condition",
         x = "Condition", y = "Design", fill = "PSI")
  
  ggsave(file.path(output_filepath_main, outfile), p, width = 5, height = 4, dpi = 300)
}

jobs <- tibble::tibble(
  index_offset = c(
    "ENSG00000164024.12;METAP1;chr4-99028866-99028918-98995658-98995867-99034229-99034342__0:0:0",
    "ENSG00000106683.15;LIMK1;chr7-74120582-74120638-74115367-74115958-74120891-74121029__0:0:0",
    "ENSG00000121064.13;SCPEP1;chr17-56985377-56985427-56981081-56981230-56991098-56991171__0:0:0",
    "ENSG00000102385.13;DRP2;chrX-101242903-101242982-101242324-101242471-101245016-101245077__0:0:0"
  ),
  cats = list(
    c("vae_Kelly_1_to_pos", "parent"),
    c("vae_Kelly_1_to_pos", "parent"),
    c("vae_Kelly_0_to_pos", "parent"),
    c("vae_Kelly_0_to_pos", "parent")
  ),
  outfile = c(
    "fig05_KELLY_specific_example_METAP1.pdf",
    "fig05_KELLY_specific_example_LIMK1.pdf",
    "fig05_KELLY_specific_example_SCPEP1.pdf",
    "fig05_KELLY_specific_example_DRP2.pdf"
  )
)

for (i in seq_len(nrow(jobs))) {
  make_heatmap(
    index_offset = jobs$index_offset[i],
    cats        = jobs$cats[[i]],
    outfile     = jobs$outfile[i]
  )
}

######## Text: Check if the designed sequences are very similar to any of the KELLY-specific positive hits that are in the training data #######

designed_df <- KELLY_designed_data %>%
  filter(significance == TRUE) %>%
  filter(design != "R2design0") %>%
  left_join(KELLY_designed_reference, by = c("index_offset_design" = "ID")) %>%
  mutate(seq = toupper(seq)) 

df47k  <- read.csv("/Volumes/broad_dawnccle/melange/figures_outputs/fig02/fig02_num_cell_type_celltype_specific.csv")
twist_lib <- read_csv("/Volumes/broad_dawnccle/melange/data/20230130_twist_library_v3.csv")

df47k <- df47k %>%
  filter(target_cell_type == "Kelly") %>%
  mutate(index_offset = sub("__.*$", "", index_offset)) %>%
  left_join(twist_lib, by = c("index_offset" = "ID")) %>%
  mutate(librarySequence = toupper(librarySequence))

lib_seqs <- df47k$librarySequence

dist_df <- designed_df %>%
  rowwise() %>%
  mutate(min_ed = min(stringdist(seq, lib_seqs, method = "lv"))) %>%
  ungroup()

write_csv(dist_df, file.path(output_filepath_supp, "designed_vs_47k_edit_distance_summary.csv"))

p <- ggplot(dist_df, aes(x = "", y = min_ed)) +
  geom_boxplot() +
  theme_classic(base_size = 10) +
  xlab(NULL) + ylab("Min edit distance to KELLY-specific sequence from unbiased library")

ggsave(file.path(output_filepath_supp, "min_edit_distance_boxplot.png"), p, width = 3, height = 5, dpi = 300)
ggsave(file.path(output_filepath_supp, "min_edit_distance_boxplot.pdf"), p, width = 3, height = 5)







