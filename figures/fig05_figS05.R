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

output_filepath_bin_switch <- "/Users/mjim/Dropbox/02Splicing/designed_sequences_validation/bin_switch/"
output_filepath_kelly_specific <- "/Users/mjim/Dropbox/02Splicing/designed_sequences_validation/KELLY_specific/"


######## Import Data: Bin switching designed sequences  ###########
PSI_file <- read.csv("/Volumes/broad_dawnccle/processed_data/reprocess_250221/count_normalized_satmutv2_dawn/satmutv2_dawn_all_samples_raw_counts.csv")
reference_file_bin_switching <- read.csv("/Users/mjim/Dropbox/02Splicing/library_satmutv2/mutagenesis_twist_library_250604_combined.csv")

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
  mutate(ID = paste(index_offset, design, sep = "__")) %>%
  left_join(reference_file_bin_switching, by = "ID") %>%
  select(-upstreamIntronSeq, -skippedExonSeq, -downstreamIntronSeq, -mut_pos, -mut_base, -mutated_region, -barcode, -twist_seq) 


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

######## Figure S5A and B: find an example sequences that bin switches properly #########
# ENSG00000073711.11;PPP2R3A;chr3-136070477-136070552-136049258-136049361-136078366-136078453__0:0:0__vae_10_to_minus10_R1design1
merged_df_example2 <- merged_df %>%
  filter(index_offset == "ENSG00000073711.11;PPP2R3A;chr3-136070477-136070552-136049258-136049361-136078366-136078453__0:0:0") %>%
  filter(category %in% c("vae_10_to_minus10", "parent")) %>%
  group_by(design, condition) %>%
  summarise(PSI = mean(PSI, na.rm = TRUE), .groups = "drop") %>%
  mutate(trailing_num = as.numeric(str_extract(design, "\\d+$"))) %>%
  # reorder factor levels of design by the trailing number
  mutate(design = factor(design, 
                         levels = rev(unique(design[order(trailing_num)])), 
                         ordered = TRUE))

plot <- ggplot(merged_df_example2, aes(x = condition, y = design, fill = PSI)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = color_palette2_custom, limits = c(0,1)) +
  theme_minimal() +
  labs(title = "PSI heatmap by design and condition",
       x = "Condition",
       y = "Design",
       fill = "PSI")
outfile <- file.path(output_filepath_bin_switch, "binswitch_example_heatmap1.pdf")
ggsave(outfile, plot, width = 5, height = 4, dpi = 300)

# ENSG00000131828.14;PDHA1;chrX-19357651-19357719-19355581-19355757-19358915-19359024__0:0:0__vae_10_to_minus10_R1design1
merged_df_example2 <- merged_df %>%
  filter(index_offset == "ENSG00000131828.14;PDHA1;chrX-19357651-19357719-19355581-19355757-19358915-19359024__0:0:0") %>%
  filter(category %in% c("vae_10_to_minus10", "parent")) %>%
  group_by(design, condition) %>%
  summarise(PSI = mean(PSI, na.rm = TRUE), .groups = "drop") %>%
  mutate(trailing_num = as.numeric(str_extract(design, "\\d+$"))) %>%
  # reorder factor levels of design by the trailing number
  mutate(design = factor(design, 
                         levels = rev(unique(design[order(trailing_num)])), 
                         ordered = TRUE))

plot <- ggplot(merged_df_example2, aes(x = condition, y = design, fill = PSI)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = color_palette2_custom, limits = c(0,1)) +
  theme_minimal() +
  labs(title = "PSI heatmap by design and condition",
       x = "Condition",
       y = "Design",
       fill = "PSI")
outfile <- file.path(output_filepath_bin_switch, "binswitch_example_heatmap2.pdf")
ggsave(outfile, plot, width = 5, height = 4, dpi = 300)

# ENSG00000228983.9;SLC47A1P1;chr17-19591764-19591860-19590743-19590788-19594209-19594334__0:0:0__vae_minus10_to_10_R1design1


merged_df_example2 <- merged_df %>%
  filter(index_offset == "ENSG00000228983.9;SLC47A1P1;chr17-19591764-19591860-19590743-19590788-19594209-19594334__0:0:0") %>%
  filter(category %in% c("vae_minus10_to_10", "parent")) %>%
  group_by(design, condition) %>%
  summarise(PSI = mean(PSI, na.rm = TRUE), .groups = "drop") %>%
  mutate(trailing_num = as.numeric(str_extract(design, "\\d+$"))) %>%
  # reorder factor levels of design by the trailing number
  mutate(design = factor(design, 
                         levels = rev(unique(design[order(trailing_num)])), 
                         ordered = TRUE))

plot <- ggplot(merged_df_example2, aes(x = condition, y = design, fill = PSI)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = color_palette2_custom, limits = c(0,1)) +
  theme_minimal() +
  labs(title = "PSI heatmap by design and condition",
       x = "Condition",
       y = "Design",
       fill = "PSI")
outfile <- file.path(output_filepath_bin_switch, "binswitch_example_heatmap3.pdf")
ggsave(outfile, plot, width = 5, height = 4, dpi = 300)

# ENSG00000083857.14;FAT1;chr4-186590367-186590403-186588391-186589220-186592692-186592752__0:0:0__vae_minus10_to_10_R1design2

merged_df_example2 <- merged_df %>%
  filter(index_offset == "ENSG00000083857.14;FAT1;chr4-186590367-186590403-186588391-186589220-186592692-186592752__0:0:0") %>%
  filter(category %in% c("vae_minus10_to_10", "parent")) %>%
  group_by(design, condition) %>%
  summarise(PSI = mean(PSI, na.rm = TRUE), .groups = "drop") %>%
  mutate(trailing_num = as.numeric(str_extract(design, "\\d+$"))) %>%
  # reorder factor levels of design by the trailing number
  mutate(design = factor(design, 
                         levels = rev(unique(design[order(trailing_num)])), 
                         ordered = TRUE))

plot <- ggplot(merged_df_example2, aes(x = condition, y = design, fill = PSI)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = color_palette2_custom, limits = c(0,1)) +
  theme_minimal() +
  labs(title = "PSI heatmap by design and condition",
       x = "Condition",
       y = "Design",
       fill = "PSI")
outfile <- file.path(output_filepath_bin_switch, "binswitch_example_heatmap4.pdf")
ggsave(outfile, plot, width = 5, height = 4, dpi = 300)



######## make graphs of position and frequency of mutations ######




seq_to_chars <- function(s) strsplit(s, "", fixed = TRUE)[[1]]

# Build tidy mutation fractions vs parent for a set of sequences in one category
# Returns tibble: position, base, frac (fraction of sequences mutated to 'base' at pos)
build_mut_tidy <- function(parent_seq, seq_vec) {
  L <- nchar(parent_seq)
  p_chars <- seq_to_chars(parent_seq)
  n <- length(seq_vec)
  if (n == 0) {
    return(tibble(position = integer(), base = character(), frac = double()))
  }
  # gather differences
  diffs <- map(seq_vec, function(s) {
    if (nchar(s) != L) return(NULL)
    v <- seq_to_chars(s)
    idx <- which(v != p_chars & v %in% c("A","C","G","T"))
    if (!length(idx)) return(NULL)
    tibble(position = idx, base = v[idx])
  }) |> list_rbind()
  
  if (is.null(diffs) || nrow(diffs) == 0) {
    return(tibble(position = integer(), base = character(), frac = double()))
  }
  
  diffs |>
    count(position, base, name = "count") |>
    mutate(frac = count / n, .keep = "unused")
}

# One plot per (index_offset, condition, category != parent)
plot_group_bars <- function(df_group) {
  stopifnot(nrow(df_group) >= 1)
  
  parent_seq <- df_group |>
    filter(category == "parent") |>
    pull(full_seq) |>
    unique()
  
  if (length(parent_seq) != 1) return(NULL)
  parent_seq <- parent_seq[[1]]
  L <- nchar(parent_seq)
  parent_chars <- seq_to_chars(parent_seq)
  
  cats <- df_group |>
    filter(category != "parent") |>
    distinct(category) |>
    pull(category)
  
  if (!length(cats)) return(NULL)
  
  map(cats, function(cat) {
    seqs <- df_group |>
      filter(category == cat) |>
      pull(full_seq)
    
    mut_tidy <- build_mut_tidy(parent_seq, seqs)
    
    # ensure all positions exist so x-axis spans full length
    base_levels <- c("A","C","G","T")
    mut_tidy_full <-
      tibble(position = 1:L) |>
      left_join(mut_tidy, by = "position") |>
      replace_na(list(frac = 0)) |>
      mutate(base = factor(base, levels = base_levels)) |>
      arrange(position, base)
    
    # --- TOP: stacked bars showing fractions (0..1)
    p_top <- ggplot(mut_tidy_full, aes(x = position, y = frac, fill = base)) +
      geom_col(width = 0.9) +
      scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
      scale_x_continuous(expand = expansion(mult = c(0,0))) +
      labs(
        title = paste0("index_offset: ", df_group$index_offset[1],
                       "   |   condition: ", df_group$condition[1],
                       "   |   category: ", cat),
        y = "Fraction mutated", x = NULL
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 9),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "top",
        legend.title = element_blank()
      )
    
    # --- BOTTOM: parent sequence letters in grey, equal height, aligned
    parent_df <- tibble(
      position = 1:L,
      base = parent_chars
    )
    
    p_bottom <- ggplot(parent_df, aes(x = position, y = 0.1, label = base)) +
      geom_text(color = "grey50", size = 3) +
      scale_y_continuous(limits = c(0, 0.2), expand = c(0,0)) +
      scale_x_continuous(limits = c(1, L), expand = expansion(mult = c(0,0))) +
      labs(x = "Position", y = NULL) +
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(t = 0, r = 5, b = 5, l = 5)
      )
    
    # stack vertically; give top more space
    p_top / p_bottom + plot_layout(heights = c(3, 1))
  })
}

# --- driver over all (index_offset, condition) groups ---
make_all_mutation_plots <- function(dat) {
  dat |>
    group_by(index_offset, condition) |>
    group_split() |>
    map(plot_group_bars) |>
    flatten() |>
    keep(~ !is.null(.x))
}




# `merged_df` must include: index_offset, condition, category, full_seq
merged_df_filtered_example1 <- merged_df_filtered %>%
  filter(index_offset == "ENSG00000057468.7;MSH4;chr1-75810696-75810807-75806980-75807141-75815020-75815136__0:0:0")
plot <- make_all_mutation_plots(merged_df_filtered_example1)
outfile <- file.path(output_filepath_kelly_specific, "0_to_1_example.pdf")
ggsave(outfile, plot[[1]], width = 10, height = 2, dpi = 300)

merged_df_filtered_example2 <- merged_df_filtered %>%
  filter(index_offset == "ENSG00000109814.12;UGDH;chr4-39503874-39503985-39498754-39500253-39504416-39504508__0:0:0")
plot <- make_all_mutation_plots(merged_df_filtered_example2)
outfile <- file.path(output_filepath_kelly_specific, "1_to_0_example.pdf")
ggsave(outfile, plot[[1]], width = 20, height = 2, dpi = 300)

all_plots <- make_all_mutation_plots(merged_df_filtered)

# Example: display the first few
print(all_plots[[7]]); 
print(all_plots[[8]])
print(all_plots[[9]]); 
print(all_plots[[10]])
print(all_plots[[11]]); 
print(all_plots[[12]])

# If you want to save them:
# for (i in seq_along(all_plots)) {
#   ggsave(sprintf("mutation_logo_%03d.png", i), all_plots[[i]], width = 12, height = 4, dpi = 200)
# }







######## Import Data: KELLY-specific designed sequences ########
raw_df <- read.csv("/Volumes/broad_dawnccle/processed_data/reprocess_250221/count_normalized_denovoR2/denovoR2_all_samples_raw_counts.csv")

reference_file_kelly_designed <- read.csv("/Users/mjim/Dropbox/02Splicing/for_kai/R2_validation/output/250718_R2_Kelly_all_sequences_to_order.csv")

# Get the unique conditions.
unique_conditions <- unique(raw_df$condition)

# Get the unique samples.
unique_samples <- unique(raw_df$sample)

raw_df_clean <- raw_df %>%
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

raw_df_to_psi_kelly <- raw_df_to_psi %>% 
  filter(grepl("KELLYspecific", condition))
raw_df_to_psi_kelly <- raw_df_to_psi_kelly %>% 
  group_by(condition, index_offset, design) %>% 
  summarise(PSI = mean(PSI)) %>% 
  ungroup() %>% 
  select(condition, index_offset, design, PSI) %>% 
  pivot_wider(names_from = condition, values_from = PSI) %>% 
  mutate(index_offset_design = paste0(index_offset, "___", design)) %>%
  mutate(KELLY_PSI_minus_HEK_PSI = `KELLYspecific-KELLY` - `KELLYspecific-HEK`) %>% 
  mutate(KELLY_PSI_minus_A375_PSI = `KELLYspecific-KELLY` - `KELLYspecific-A375`) %>% 
  mutate(KELLY_PSI_minus_T47D_PSI = `KELLYspecific-KELLY` - `KELLYspecific-T47D`)

######## Import Pairadise data for kelly designed sequences ##########
kelly_pairadise <- read_tsv("/Volumes/broad_dawnccle/processed_data/reprocess_250221/pairadise_denovo_KELLY_PSI/WT_KELLYspecific-KELLY_all/KELLYspecific-KELLY_all_rMATS_Result_P.FDR.txt", col_types = cols(S1 = col_character(), S2 = col_character()))

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


# Volcano plot
plot <- ggplot(kelly_pairadise_PSI, aes(x = deltaPSI, y = -log10(FDR))) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # FDR cutoff
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "blue") + # effect size cutoffs (example)
  labs(
    title = "Volcano Plot of KELLY Pairadise PSI",
    x = expression(Delta*"PSI"),
    y = expression(-log[10]*"(FDR)")
  ) +
  ylim(0, 15) +
  xlim(-0.75, 0.75) + 
  theme_minimal()
outfile <- file.path(output_filepath_kelly_specific, "kelly_designed_volcano.pdf")
ggsave(outfile, plot, width = 5.5, height = 5.0, dpi = 300)
kelly_pairadise_PSI$ExonID <- sub("__(?!.*__)", "___", kelly_pairadise_PSI$ExonID, perl = TRUE)

 
########add in pairadise stats ##########
KELLY_PSI_filtered <- raw_df_to_psi_kelly %>%
  mutate(across(c(`KELLYspecific-KELLY`, `KELLYspecific-HEK`, `KELLYspecific-A375`, `KELLYspecific-T47D`), ~ as.numeric(.))) %>%
  rowwise() %>%
  mutate(upsilon = calculate_upsilon(c(`KELLYspecific-KELLY`, `KELLYspecific-HEK`, `KELLYspecific-A375`, `KELLYspecific-T47D`))) %>%
  mutate(reverse_upsilon = calculate_reverse_upsilon(c(`KELLYspecific-KELLY`, `KELLYspecific-HEK`, `KELLYspecific-A375`, `KELLYspecific-T47D`))) %>%
  ungroup() %>%
  left_join(kelly_pairadise_PSI, by = c("index_offset_design" = "ExonID"))

KELLY_PSI_filtered <- KELLY_PSI_filtered %>%
  mutate(
    category = if_else(str_detect(design, "_"),
                       str_replace(design, "_[^_]+$", ""),   # everything before last "_"
                       "parent"),
    R1design = if_else(str_detect(design, "_"),
                       str_extract(design, "[^_]+$"),        # everything after last "_"
                       design)
  )
######## make examples of sequnces that are kelly specific ##########

# Example 1: ENSG00000164024.12;METAP1;chr4-99028866-99028918-98995658-98995867-99034229-99034342__0:0:0___vae_Kelly_1_to_pos_R2design6
merged_df_example1 <- KELLY_PSI_filtered %>%
  filter(index_offset == "ENSG00000164024.12;METAP1;chr4-99028866-99028918-98995658-98995867-99034229-99034342__0:0:0") %>%
  filter(category %in% c("vae_Kelly_1_to_pos", "parent")) %>%
  filter(significance | category == "parent") %>%
  # pivot the four condition columns into long format
  pivot_longer(cols = c(`KELLYspecific-A375`, `KELLYspecific-HEK`,
                        `KELLYspecific-KELLY`, `KELLYspecific-T47D`),
               names_to = "condition", values_to = "PSI") %>%
  # extract trailing number from design
  mutate(trailing_num = as.numeric(str_extract(design, "\\d+$"))) %>%
  # reorder design in descending order
  mutate(design = factor(design,
                         levels = unique(design[order(-trailing_num)]),
                         ordered = TRUE))

plot <- ggplot(merged_df_example1, aes(x = condition, y = design, fill = PSI)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = color_palette2_custom, limits = c(0,1)) +
  theme_minimal() +
  labs(title = "PSI heatmap by design and condition",
       x = "Condition",
       y = "Design",
       fill = "PSI")
outfile <- file.path(output_filepath_kelly_specific, "example1_heatmap.pdf")
ggsave(outfile, plot, width = 5, height = 4, dpi = 300)

merged_df_example1 <- merged_df_example1 %>%
  left_join(reference, by = c("ExonID" = "ID"))



# Example 2: ENSG00000106683.15;LIMK1;chr7-74120582-74120638-74115367-74115958-74120891-74121029__0:0:0___vae_Kelly_1_to_pos_R2design3
merged_df_example1 <- KELLY_PSI_filtered %>%
  filter(index_offset == "ENSG00000106683.15;LIMK1;chr7-74120582-74120638-74115367-74115958-74120891-74121029__0:0:0") %>%
  filter(category %in% c("vae_Kelly_1_to_pos", "parent")) %>%
  filter(significance | category == "parent") %>%
  # pivot the four condition columns into long format
  pivot_longer(cols = c(`KELLYspecific-A375`, `KELLYspecific-HEK`,
                        `KELLYspecific-KELLY`, `KELLYspecific-T47D`),
               names_to = "condition", values_to = "PSI") %>%
  # extract trailing number from design
  mutate(trailing_num = as.numeric(str_extract(design, "\\d+$"))) %>%
  # reorder design in descending order
  mutate(design = factor(design,
                         levels = unique(design[order(-trailing_num)]),
                         ordered = TRUE))

plot <- ggplot(merged_df_example1, aes(x = condition, y = design, fill = PSI)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = color_palette2_custom, limits = c(0,1)) +
  theme_minimal() +
  labs(title = "PSI heatmap by design and condition",
       x = "Condition",
       y = "Design",
       fill = "PSI")
outfile <- file.path(output_filepath_kelly_specific, "example2_heatmap.pdf")
ggsave(outfile, plot, width = 5, height = 4, dpi = 300)

merged_df_example1 <- merged_df_example1 %>%
  left_join(reference, by = c("ExonID" = "ID"))

# Example 3: #ENSG00000121064.13;SCPEP1;chr17-56985377-56985427-56981081-56981230-56991098-56991171__0:0:0___vae_Kelly_0_to_pos_R2design6

merged_df_example1 <- KELLY_PSI_filtered %>%
  filter(index_offset == "ENSG00000121064.13;SCPEP1;chr17-56985377-56985427-56981081-56981230-56991098-56991171__0:0:0") %>%
  filter(category %in% c("vae_Kelly_0_to_pos", "parent")) %>%
  filter(significance | category == "parent") %>%
  # pivot the four condition columns into long format
  pivot_longer(cols = c(`KELLYspecific-A375`, `KELLYspecific-HEK`,
                        `KELLYspecific-KELLY`, `KELLYspecific-T47D`),
               names_to = "condition", values_to = "PSI") %>%
  # extract trailing number from design
  mutate(trailing_num = as.numeric(str_extract(design, "\\d+$"))) %>%
  # reorder design in descending order
  mutate(design = factor(design,
                         levels = unique(design[order(-trailing_num)]),
                         ordered = TRUE))

plot <- ggplot(merged_df_example1, aes(x = condition, y = design, fill = PSI)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = color_palette2_custom, limits = c(0,1)) +
  theme_minimal() +
  labs(title = "PSI heatmap by design and condition",
       x = "Condition",
       y = "Design",
       fill = "PSI")
outfile <- file.path(output_filepath_kelly_specific, "example3_heatmap.pdf")
ggsave(outfile, plot, width = 5, height = 4, dpi = 300)

merged_df_example1 <- merged_df_example1 %>%
  left_join(reference, by = c("ExonID" = "ID"))

# Example 4: ENSG00000102385.13;DRP2;chrX-101242903-101242982-101242324-101242471-101245016-101245077__0:0:0___vae_Kelly_0_to_pos_R2design4
merged_df_example1 <- KELLY_PSI_filtered %>%
  filter(index_offset == "ENSG00000102385.13;DRP2;chrX-101242903-101242982-101242324-101242471-101245016-101245077__0:0:0") %>%
  filter(category %in% c("vae_Kelly_0_to_pos", "parent")) %>%
  filter(significance | category == "parent") %>%
  # pivot the four condition columns into long format
  pivot_longer(cols = c(`KELLYspecific-A375`, `KELLYspecific-HEK`,
                        `KELLYspecific-KELLY`, `KELLYspecific-T47D`),
               names_to = "condition", values_to = "PSI") %>%
  # extract trailing number from design
  mutate(trailing_num = as.numeric(str_extract(design, "\\d+$"))) %>%
  # reorder design in descending order
  mutate(design = factor(design,
                         levels = unique(design[order(-trailing_num)]),
                         ordered = TRUE))

plot <- ggplot(merged_df_example1, aes(x = condition, y = design, fill = PSI)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = color_palette2_custom, limits = c(0,1)) +
  theme_minimal() +
  labs(title = "PSI heatmap by design and condition",
       x = "Condition",
       y = "Design",
       fill = "PSI")
outfile <- file.path(output_filepath_kelly_specific, "example4_heatmap.pdf")
ggsave(outfile, plot, width = 5, height = 4, dpi = 300)

merged_df_example1 <- merged_df_example1 %>%
  left_join(reference, by = c("ExonID" = "ID"))


######## calculate upsilon and reverse upsilon scores ##########
# plot upsilon scores here: 
KELLY_PSI_filtered <- KELLY_PSI_filtered %>%
  rowwise() %>%
  mutate(
    max_cell_line = colnames(across(`KELLYspecific-A375`:`KELLYspecific-T47D`))[ 
      which.max(c_across(`KELLYspecific-A375`:`KELLYspecific-T47D`)) 
    ],
    min_cell_line = colnames(across(`KELLYspecific-A375`:`KELLYspecific-T47D`))[ 
      which.min(c_across(`KELLYspecific-A375`:`KELLYspecific-T47D`)) 
    ]
  ) %>%
  ungroup() 
plot <- ggplot(KELLY_PSI_filtered, aes(x = upsilon)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", color = "white") +
  labs(
    x = "Upsilon Score",
    y = "Count",
    title = "Distribution of Upsilon Scores"
  ) +
  scale_y_log10() +
  theme_minimal()
outfile <- file.path(output_filepath_kelly_specific, "upsilon_score_histogram.pdf")
ggsave(outfile, plot, width = 5.5, height = 5.0, dpi = 300)

######## make a heatmap of the top sequences ######

# Select the relevant columns
cols_of_interest <- c("KELLYspecific-KELLY", 
                      "KELLYspecific-HEK", 
                      "KELLYspecific-A375", 
                      "KELLYspecific-T47D")

# Filter for significant rows only
df_sig <- KELLY_PSI_filtered %>%
  filter(significance == TRUE)

# Extract and convert to matrix
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

outfile <- file.path(output_filepath_kelly_specific, "significant_hits_heatmap.pdf")
ggsave(outfile, plot, width = 5.5, height = 5.0, dpi = 300)
######## histogram of upsilon scores ###########



### overlap the parent sequences on top
plot <- ggplot(KELLY_PSI_filtered, aes(x = upsilon)) +
  # background distribution
  geom_histogram(
    data = KELLY_PSI_filtered,
    binwidth = 0.05,
    fill = "steelblue",
    alpha = 0.6
  ) +
  # highlight R2design0
  geom_histogram(
    data = subset(KELLY_PSI_filtered, design == "R2design0"),
    binwidth = 0.05,
    fill = "red",
    alpha = 0.6
  ) +
  labs(
    x = "Upsilon Score",
    y = "Count",
    title = "Distribution of Upsilon Scores"
  ) +
  scale_y_log10() +
  theme_minimal()
outfile <- file.path(output_filepath_kelly_specific, "upsilon_score_histogram_parent_overlap.pdf")
ggsave(outfile, plot, width = 5.5, height = 5.0, dpi = 300)

plot <- ggplot(KELLY_PSI_filtered, aes(x = reverse_upsilon)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", color = "white") +
  labs(
    x = "Reverse Upsilon Score",
    y = "Count",
    title = "Distribution of Reverse Upsilon Scores"
  ) +
  scale_y_log10() +
  theme_minimal()
outfile <- file.path(output_filepath_kelly_specific, "reverse_upsilon_score_histogram.pdf")
ggsave(outfile, plot, width = 5.5, height = 5.0, dpi = 300)

### overlap the parent sequences on top
plot <- ggplot(KELLY_PSI_filtered, aes(x = reverse_upsilon)) +
  # background distribution
  geom_histogram(
    data = KELLY_PSI_filtered,
    binwidth = 0.05,
    fill = "steelblue",
    alpha = 0.6
  ) +
  # highlight R2design0
  geom_histogram(
    data = subset(KELLY_PSI_filtered, design == "R2design0"),
    binwidth = 0.05,
    fill = "red",
    alpha = 0.6
  ) +
  labs(
    x = "reverse Upsilon Score",
    y = "Count",
    title = "Distribution of reverse Upsilon Scores"
  ) +
  scale_y_log10() +
  theme_minimal()
outfile <- file.path(output_filepath_kelly_specific, "reverse_upsilon_score_histogram_parent_overlap.pdf")
ggsave(outfile, plot, width = 5.5, height = 5.0, dpi = 300)

######## add in sequences and look for RBFOX1 #######

## add in sequences from reference
KELLY_PSI_filtered_sequences <- KELLY_PSI_filtered %>%
  mutate(index_offset_design = str_replace_all(index_offset_design, "___", "__")) %>%
  left_join(reference_file_kelly_designed, by = c("index_offset_design" = "ID"))


KELLY_PSI_filtered_sequences <- KELLY_PSI_filtered_sequences %>%
  mutate(
    motif_count = str_count(seq, "TGCATG"),
    group = case_when(
      design == "R2design0" ~ "R2design0",
      significance == TRUE ~ "Sig",
      TRUE ~ "Other"
    )
  )

# Summarize fraction with the motif
plot_df <- KELLY_PSI_filtered_sequences %>%
  group_by(group) %>%
  summarise(
    frac_with_motif = mean(motif_count > 0),  # fraction of sequences that have ≥1 motif
    .groups = "drop"
  )

# Plot
ggplot(plot_df, aes(x = group, y = frac_with_motif, fill = group)) +
  geom_col(width = 0.6) +
  labs(
    x = "Design",
    y = "Fraction with RBFOX1 motif",
    title = "RBFOX1 motif frequency"
  ) +
  theme_minimal() +
  theme(legend.position = "none")


######## Do k-mer analysis ########

# make a binary grouping: Sig vs Background
df_for_enrich <- KELLY_PSI_filtered_sequences %>%
  mutate(group_bin = if_else(group == "Sig", "Sig", "Background"))

# Helper to extract all kmers from a sequence
get_kmers <- function(seq, k) {
  n <- nchar(seq)
  if (n < k) return(character(0))
  sapply(1:(n - k + 1), function(i) substr(seq, i, i + k - 1))
}

# Function to build enrichment table for given k
kmer_enrichment <- function(df, k) {
  
  kmer_df <- df %>%
    rowwise() %>%
    mutate(kmer = list(unique(get_kmers(seq, k)))) %>%
    unnest(kmer) %>%
    ungroup()
  
  counts <- kmer_df %>%
    group_by(group_bin, kmer) %>%
    summarise(seq_with_kmer = n(), .groups = "drop") %>%
    pivot_wider(names_from = group_bin, values_from = seq_with_kmer, values_fill = 0)
  
  total_sig <- df %>% filter(group_bin == "Sig") %>% nrow()
  total_bg  <- df %>% filter(group_bin == "Background") %>% nrow()
  
  results <- counts %>%
    rowwise() %>%
    mutate(
      ft = list(
        fisher.test(
          matrix(c(Sig,
                   total_sig - Sig,
                   Background,
                   total_bg - Background),
                 nrow = 2)
        )
      ),
      pval = ft$p.value,
      odds_ratio = unname(ft$estimate)
    ) %>%
    ungroup() %>%
    mutate(
      fdr = p.adjust(pval, method = "BH"),
      log2_or = log2(odds_ratio),
      neglog10_fdr = -log10(fdr)
    ) %>%
    select(kmer, Sig, Background, pval, fdr, odds_ratio, log2_or, neglog10_fdr)
  
  return(results)
}


# Volcano plotting function
plot_kmer_volcano <- function(enrich_df, k, fdr_cutoff = 0.05) {
  
  # pick top kmers based on significance + direction
  sig_df <- enrich_df %>% filter(fdr < fdr_cutoff)
  
  # top 5 enriched, top 5 depleted
  top_enriched <- sig_df %>% arrange(desc(log2_or)) %>% head(15)
  top_depleted <- sig_df %>% arrange(log2_or) %>% head(15)
  label_df <- bind_rows(top_enriched, top_depleted)
  
  ggplot(enrich_df, aes(x = log2_or, y = neglog10_fdr)) +
    geom_point(aes(color = fdr < fdr_cutoff), alpha = 0.7) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey70")) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_text_repel(data = label_df, aes(label = kmer),
                    size = 3, max.overlaps = 50) +
    labs(
      x = "log2(Odds Ratio)",
      y = "-log10(FDR)",
      title = paste0(k, "-mer Enrichment Volcano Plot"),
      color = paste0("FDR <", fdr_cutoff)
    ) +
    theme_minimal()
}

# --- Run for 3–6mers ---
df_for_enrich <- KELLY_PSI_filtered_sequences %>%
  mutate(group_bin = if_else(group == "Sig", "Sig", "Background"))

enrich3 <- kmer_enrichment(df_for_enrich, 3)
enrich4 <- kmer_enrichment(df_for_enrich, 4)
enrich5 <- kmer_enrichment(df_for_enrich, 5)
enrich6 <- kmer_enrichment(df_for_enrich, 6)

plot3 <- plot_kmer_volcano(enrich3, 3)
plot4 <- plot_kmer_volcano(enrich4, 4)
plot5 <- plot_kmer_volcano(enrich5, 5)
plot6 <- plot_kmer_volcano(enrich6, 6)


outfile <- file.path(output_filepath_kelly_specific, "3mer_volcano.pdf")
ggsave(outfile, plot3, width = 5.5, height = 5.0, dpi = 300)
outfile <- file.path(output_filepath_kelly_specific, "4mer_volcano.pdf")
ggsave(outfile, plot4, width = 5.5, height = 5.0, dpi = 300)
outfile <- file.path(output_filepath_kelly_specific, "5mer_volcano.pdf")
ggsave(outfile, plot5, width = 5.5, height = 5.0, dpi = 300)
outfile <- file.path(output_filepath_kelly_specific, "6mer_volcano.pdf")
ggsave(outfile, plot6, width = 5.5, height = 5.0, dpi = 300)
######## export sequences for streame ########

## add in individual sequences split by region
downstream_path <- "/Users/mjim/Dropbox/02Splicing/for_kai/R2_validation/output/RNA_denovo_R2_250718_downstreamIntronSeq.fasta"
skipped_path    <- "/Users/mjim/Dropbox/02Splicing/for_kai/R2_validation/output/RNA_denovo_R2_250718_skippedExonSeq.fasta"
upstream_path   <- "/Users/mjim/Dropbox/02Splicing/for_kai/R2_validation/output/RNA_denovo_R2_250718_upstreamIntronSeq.fasta"
read_fasta_to_df <- function(path, colname) {
  seqs <- readDNAStringSet(path)
  tibble(
    index_offset_design = names(seqs),
    !!colname := as.character(seqs)
  )
}
downstream_df <- read_fasta_to_df(downstream_path, "downstreamIntronSeq")
skipped_df    <- read_fasta_to_df(skipped_path, "skippedExonSeq")
upstream_df   <- read_fasta_to_df(upstream_path, "upstreamIntronSeq")
KELLY_PSI_filtered_sequences <- KELLY_PSI_filtered_sequences %>%
  left_join(downstream_df, by = "index_offset_design") %>%
  left_join(skipped_df,    by = "index_offset_design") %>%
  left_join(upstream_df,   by = "index_offset_design")

### make fasta files to put into STREME

# Filter sequences with upsilon > 0.2
high_upsilon <- KELLY_PSI_filtered_sequences %>%
  filter(significance == TRUE)%>%
  filter(design != "R2design0") %>%
  group_by(index_offset) %>%
  slice(1) %>%
  ungroup()

# Create FASTA for high upsilon sequences
high_fasta <- DNAStringSet(high_upsilon$seq)
names(high_fasta) <- paste0("seq_", seq_along(high_fasta))
writeXStringSet(high_fasta, "/Users/mjim/Dropbox/02Splicing/designed_sequences_validation/high_upsilon_sequences.fasta")

# Find one control (R2design0) sequence per high-upsilon row
controls <- KELLY_PSI_filtered_sequences %>%
  filter(design == "R2design0" & index_offset %in% high_upsilon$index_offset) 

# Create FASTA for control sequences
control_fasta <- DNAStringSet(controls$seq)
names(control_fasta) <- paste0("control_", seq_along(control_fasta))
writeXStringSet(control_fasta, "/Users/mjim/Dropbox/02Splicing/designed_sequences_validation/control_sequences.fasta")


#### check for frequency of STREME-identified motifs in the dataset

#motifs_RNA <- c("UUGCACGC", "AGUUAGUU", "CCUAGGGA", "RUAGUCAG", "MAGGUAAA", "GCAGGRAG", "UUUCU", "GCAGUCA", "GACUAAU", "AUUUUUCU")
motifs_DNA <- c("TTGCACGC", "AGTTAGTT", "CCTAGGGA", "RTAGTCAG", "MAGGTAAA", "GCAGGRAG", "TTTCT", "GCAGTCA", "GACTAAT", "ATTTTTCT")

# IUPAC code mapping to regex
iupac_map <- c(
  A="A", C="C", G="G", T="T", U="T",
  R="[AG]", Y="[CT]", K="[GT]", M="[AC]", S="[GC]", W="[AT]",
  B="[GTC]", D="[GAT]", H="[ACT]", V="[GCA]", N="[ACGT]"
)

# Convert motif to regex pattern
motif_to_regex <- function(motif) {
  paste0(sapply(strsplit(motif, "")[[1]], function(base) iupac_map[[base]]), collapse = "")
}

# Triangular smoothing weights
triangular_weights <- function(n) {
  mid <- ceiling(n / 2)
  w <- c(1:mid, rev(1:(n - mid)))
  w / sum(w)
}

# Frequency by position with IUPAC support + triangular smoothing
motif_freq_by_pos <- function(df, region_col, motifs) {
  seqs <- df[[region_col]]
  max_len <- max(nchar(seqs), na.rm = TRUE)
  window <- max(3, ceiling(max_len * 0.1))  # 5% width, min size 3
  weights <- triangular_weights(window)
  
  map_dfr(motifs, function(motif) {
    motif_regex <- motif_to_regex(motif)
    raw_counts <- map_dbl(1:max_len, ~{
      sum(str_detect(str_sub(seqs, ., . + nchar(motif) - 1), paste0("^", motif_regex)), na.rm = TRUE)
    })
    
    smoothed_counts <- zoo::rollapply(
      raw_counts,
      width = window,
      FUN = function(x) sum(x * weights),
      align = "center",
      fill = NA
    )
    
    tibble(
      position = 1:max_len,
      freq = smoothed_counts,
      motif = motif,
      region = region_col
    )
  })
}

# Combine all three regions
freq_df <- bind_rows(
  motif_freq_by_pos(KELLY_PSI_filtered_sequences, "upstreamIntronSeq", motifs_DNA),
  motif_freq_by_pos(KELLY_PSI_filtered_sequences, "skippedExonSeq", motifs_DNA),
  motif_freq_by_pos(KELLY_PSI_filtered_sequences, "downstreamIntronSeq", motifs_DNA)
)

# Plot smoothed motif frequency by position
ggplot(freq_df, aes(x = position, y = freq, color = motif)) +
  geom_line(size = 0.8, alpha = 0.9) +
  facet_wrap(~region, scales = "free_x") +
  labs(
    x = "Position",
    y = "Smoothed Motif Frequency",
    title = "Smoothed Motif Frequency by Position (IUPAC-aware)",
    color = "Motif"
  ) +
  theme_minimal()

## t test ##
# IUPAC regex map
iupac_map <- c(
  A="A", C="C", G="G", T="T", U="T",
  R="[AG]", Y="[CT]", K="[GT]", M="[AC]", S="[GC]", W="[AT]",
  B="[GTC]", D="[GAT]", H="[ACT]", V="[GCA]", N="[ACGT]"
)

motif_to_regex <- function(motif) {
  paste0(sapply(strsplit(motif, "")[[1]], function(base) iupac_map[[base]]), collapse = "")
}

# For each motif, run t-test of upsilon values in motif+ vs motif-
t_results <- map_dfr(motifs_DNA, function(motif) {
  motif_regex <- motif_to_regex(motif)
  
  df_counts <- KELLY_PSI_filtered_sequences %>%
    mutate(has_motif = str_detect(seq, motif_regex))
  
  high <- df_counts$upsilon[df_counts$has_motif]
  low  <- df_counts$upsilon[!df_counts$has_motif]
  
  # Skip if not enough variation
  if (length(unique(df_counts$has_motif)) < 2) {
    return(tibble(motif = motif, t_stat = NA, p_value = NA))
  }
  
  test <- t.test(high, low, alternative = "greater")
  
  tibble(
    motif = motif,
    t_stat = unname(test$statistic),
    p_value = test$p.value
  )
}) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

print(t_results)

# Prepare data for plotting: mean upsilon for motif+ and motif-
plot_df <- map_dfr(motifs_DNA, function(motif) {
  motif_regex <- motif_to_regex(motif)
  
  KELLY_PSI_filtered_sequences %>%
    mutate(
      has_motif = if_else(str_detect(seq, motif_regex), "Motif+", "Motif-")
    ) %>%
    group_by(has_motif) %>%
    summarise(mean_upsilon = mean(upsilon, na.rm = TRUE), .groups = "drop") %>%
    mutate(motif = motif)
})

# Bar plot: upsilon values by motif presence
ggplot(plot_df, aes(x = motif, y = mean_upsilon, fill = has_motif)) +
  geom_col(position = "dodge") +
  labs(
    x = "Motif",
    y = "Mean Upsilon Score",
    title = "Upsilon Scores by Motif Presence"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



######## make barchart of the frequency of each mutation on the kelly-specific sequences ########
KELLY_PSI_filtered_sequences <- KELLY_PSI_filtered_sequences %>%
  mutate(
    category = if_else(str_detect(design, "_"),
                       str_replace(design, "_[^_]+$", ""),   # everything before last "_"
                       "parent"),
    R1design = if_else(str_detect(design, "_"),
                       str_extract(design, "[^_]+$"),        # everything after last "_"
                       design)
  ) 


seq_to_chars <- function(s) strsplit(s, "", fixed = TRUE)[[1]]

# Build tidy mutation fractions vs parent for a set of sequences in one category
# Returns tibble: position, base, frac (fraction of sequences mutated to 'base' at pos)
build_mut_tidy <- function(parent_seq, seq_vec) {
  L <- nchar(parent_seq)
  p_chars <- seq_to_chars(parent_seq)
  n <- length(seq_vec)
  if (n == 0) {
    return(tibble(position = integer(), base = character(), frac = double()))
  }
  # gather differences
  diffs <- map(seq_vec, function(s) {
    if (nchar(s) != L) return(NULL)
    v <- seq_to_chars(s)
    idx <- which(v != p_chars & v %in% c("A","C","G","T"))
    if (!length(idx)) return(NULL)
    tibble(position = idx, base = v[idx])
  }) |> list_rbind()
  
  if (is.null(diffs) || nrow(diffs) == 0) {
    return(tibble(position = integer(), base = character(), frac = double()))
  }
  
  diffs |>
    count(position, base, name = "count") |>
    mutate(frac = count / n, .keep = "unused")
}

# One plot per (index_offset, condition, category != parent)
plot_group_bars <- function(df_group) {
  stopifnot(nrow(df_group) >= 1)
  
  parent_seq <- df_group |>
    filter(category == "parent") |>
    pull(full_seq) |>
    unique()
  
  if (length(parent_seq) != 1) return(NULL)
  parent_seq <- parent_seq[[1]]
  L <- nchar(parent_seq)
  parent_chars <- seq_to_chars(parent_seq)
  
  cats <- df_group |>
    filter(category != "parent") |>
    distinct(category) |>
    pull(category)
  
  if (!length(cats)) return(NULL)
  
  map(cats, function(cat) {
    seqs <- df_group |>
      filter(category == cat) |>
      pull(full_seq)
    
    mut_tidy <- build_mut_tidy(parent_seq, seqs)
    
    # ensure all positions exist so x-axis spans full length
    base_levels <- c("A","C","G","T")
    mut_tidy_full <-
      tibble(position = 1:L) |>
      left_join(mut_tidy, by = "position") |>
      replace_na(list(frac = 0)) |>
      mutate(base = factor(base, levels = base_levels)) |>
      arrange(position, base)
    
    # --- TOP: stacked bars showing fractions (0..1)
    p_top <- ggplot(mut_tidy_full, aes(x = position, y = frac, fill = base)) +
      geom_col(width = 0.9) +
      scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
      scale_x_continuous(expand = expansion(mult = c(0,0))) +
      labs(
        title = paste0("index_offset: ", df_group$index_offset[1],
                       "   |   condition: ", df_group$condition[1],
                       "   |   category: ", cat),
        y = "Fraction mutated", x = NULL
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 9),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "top",
        legend.title = element_blank()
      )
    
    # --- BOTTOM: parent sequence letters in grey, equal height, aligned
    parent_df <- tibble(
      position = 1:L,
      base = parent_chars
    )
    
    p_bottom <- ggplot(parent_df, aes(x = position, y = 0.1, label = base)) +
      geom_text(color = "grey50", size = 3) +
      scale_y_continuous(limits = c(0, 0.2), expand = c(0,0)) +
      scale_x_continuous(limits = c(1, L), expand = expansion(mult = c(0,0))) +
      labs(x = "Position", y = NULL) +
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(t = 0, r = 5, b = 5, l = 5)
      )
    
    # stack vertically; give top more space
    p_top / p_bottom + plot_layout(heights = c(3, 1))
  })
}

# --- driver over all (index_offset, condition) groups ---
make_all_mutation_plots <- function(dat) {
  dat |>
    group_by(index_offset) |>
    group_split() |>
    map(plot_group_bars) |>
    flatten() |>
    keep(~ !is.null(.x))
}




# `merged_df` must include: index_offset, condition, category, full_seq
# (Your earlier pipeline already had these.)
all_plots <- make_all_mutation_plots(KELLY_PSI_filtered_sequences)

# Example: display the first few
print(all_plots[[20]]); print(all_plots[[21]])








######## ARCHIEVE ##########
######## Figure S5C, D example cell type specific designed sequences  - OLD ############

merged_df_example2 <- KELLY_PSI_filtered %>%
  filter(index_offset == "ENSG00000258216.7;RP11-654D12.2;chr12-90094435-90094479-89989223-89989289-90100702-90100763__0:0:0") %>%
  filter(category %in% c("vae_Kelly_1_to_pos", "parent")) %>%
  filter(significance | category == "parent") %>%
  # pivot the four condition columns into long format
  pivot_longer(cols = c(`KELLYspecific-A375`, `KELLYspecific-HEK`,
                        `KELLYspecific-KELLY`, `KELLYspecific-T47D`),
               names_to = "condition", values_to = "PSI") %>%
  # extract trailing number from design
  mutate(trailing_num = as.numeric(str_extract(design, "\\d+$"))) %>%
  # reorder design in descending order
  mutate(design = factor(design,
                         levels = unique(design[order(-trailing_num)]),
                         ordered = TRUE))

plot <- ggplot(merged_df_example2, aes(x = condition, y = design, fill = PSI)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = color_palette2_custom, limits = c(0,1)) +
  theme_minimal() +
  labs(title = "PSI heatmap by design and condition",
       x = "Condition",
       y = "Design",
       fill = "PSI")
outfile <- file.path(output_filepath_kelly_specific, "RP11-654D12_example_barchart.pdf")
ggsave(outfile, plot, width = 5, height = 4, dpi = 300)

## add in the sequences themselves
merged_df_example2 <- merged_df_example2 %>%
  mutate(ID = sub("___", "__", index_offset_design)) %>%
  left_join(reference_file_kelly_designed, by = c("ID"))

