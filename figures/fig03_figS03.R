library(dplyr)
library(ggplot2)
library(tidyr)
library(vroom)
library(stringr)
library(tidyverse)
library(vroom)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(pheatmap)
library(patchwork)

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


######## Data Processing: Import original MPRA data ########
original_MPRA_data <- read_csv("/Volumes/broad_dawnccle/melange/data/250527_FINAL_PSI_TABLE_by_condition.csv")
original_MPRA_data <- original_MPRA_data %>%
  dplyr::rename(PSI_original = PSI)




######## Data Processing: Import all the KELLY specific sequences and make merged_PSI df  ########
PSI_file <- read.csv("/Volumes/broad_dawnccle/melange/data/satmut_KELLY_specific_count_table_normalized.csv")
reference_file <- read.csv("/Volumes/broad_dawnccle/melange/data/satmut_KELLY_specific_reference.csv")

PSI_file_clean <- PSI_file %>% 
  mutate(index_offset = paste0(index, "__", offset_initial)) %>% 
  filter(offset_initial == offset | offset == 0) %>%  ####new added in 250721
  select(-filename, - index, - offset_initial) %>%
  filter(mode %in% c("INCLUDED", "SKIPPED")) %>% 
  group_by(sample, condition, index_offset, mode, offset, location, base) %>%
  summarise(count = sum(count))

PSI_file_clean_to_PSI <- PSI_file_clean %>% 
  group_by(sample, condition, index_offset, mode, location, base) %>%
  summarise(count = sum(count)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = mode, values_from = count, values_fill = 0) %>%
  mutate(total_counts = INCLUDED + SKIPPED) %>%
  filter(total_counts >= 30) %>%
  mutate(PSI = INCLUDED / (INCLUDED + SKIPPED)) %>%
  mutate(position = str_extract(base, "\\d+")) %>% 
  mutate(position = as.integer(position)) %>%
  mutate(nucleotide = str_extract(base, "[A-Z]$")) 

PSI_file_clean_to_PSI_grouped <- PSI_file_clean_to_PSI %>%
  group_by(index_offset, condition, position, nucleotide, location) %>% 
  summarise(PSI = mean(PSI), .groups = "drop")

# Calculate the average PSI for each gene, intron, and position group
parent_psi <- PSI_file_clean_to_PSI_grouped %>%
  filter(location == "extra3" | is.na(position)) %>%
  mutate(par_PSI = PSI) %>%
  select(condition, index_offset, par_PSI)


non_parent_psi <- PSI_file_clean_to_PSI_grouped %>% 
  filter(location != "extra3") 

merged_PSI <- merge(non_parent_psi, parent_psi, by = c("index_offset", "condition"), all.x = T) %>% 
  mutate(location = factor(location, levels = c("upstream", "exon", "downstream", "extra3"))) %>%
  mutate(plot_position = ifelse(location == "upstream", -position, position)) %>%
  mutate(index_offset_label = sapply(strsplit(as.character(index_offset), ";"), `[`, 2))




######## Figure S3A - KELLY-specific heatmap from original MPRA dataset ##########
original_MPRA_data_filtered <- original_MPRA_data %>%
  pivot_wider(names_from = condition, values_from = PSI_original, values_fill = NA) %>%
  filter(index_offset %in% c(
    "ENSG00000138821.14;SLC39A8;chr4-102285948-102285973-102267871-102268079-102304316-102304481__0:0:0",
    "ENSG00000135365.16;PHF21A;chr11-45946075-45946098-45938156-45938312-45948885-45948946__0:0:0",
    "ENSG00000172292.15;CERS6;chr2-168766321-168766345-168765591-168765748-168769509-168775134__0:0:0",
    "ENSG00000111850.11;SMIM8;chr6-87330691-87330712-87322587-87322632-87337008-87337166__0:0:0",
    "ENSG00000075711.21;DLG1;chr3-197075836-197075870-197066703-197066754-197076585-197076682__0:0:0",
    "ENSG00000258826.6;LINC01147;chr14-88023674-88023776-88018180-88018703-88026087-88026334__11:0:0",
    "ENSG00000157107.14;FCHO2;chr5-73074741-73074853-73068649-73068779-73077337-73077440__0:0:0",
    "ENSG00000136436.15;CALCOCO2;chr17-48856578-48856616-48856100-48856187-48860313-48860449__0:0:0",
    "ENSG00000088538.13;DOCK3;chr3-51333000-51333027-51330137-51330223-51333157-51333253__0:0:0",
    "ENSG00000196739.15;COL27A1;chr9-114231821-114231858-114231078-114231132-114235598-114235652__0:0:0",
    "ENSG00000136717.15;BIN1;chr2-127053421-127053445-127051153-127051243-127053904-127054012__0:0:0",
    "ENSG00000100592.16;DAAM1;chr14-59338384-59338414-59331812-59331920-59340073-59340180__0:0:0",
    "ENSG00000176986.16;SEC24C;chr10-73762097-73762166-73760712-73760849-73763489-73763555__0:0:0",
    "ENSG00000169224.13;GCSAML;chr1-247538675-247538791-247526939-247527054-247549044-247549220__0:0:0",
    "ENSG00000156931.16;VPS8;chr3-184838713-184838746-184834648-184834742-184849070-184849195__0:0:0"
  ))
original_MPRA_data_filtered <- original_MPRA_data_filtered %>%
  mutate(
    parts = str_split(index_offset, ";"),
    gene_name = map_chr(parts, 2),
    last_group = map_chr(parts, ~ str_extract(.x[3], "chr[^_]+")),  # e.g., chr4-102285948-...
    chr_parts = str_split(last_group, "-"),
    chr = map_chr(chr_parts, 1),
    start = map_chr(chr_parts, 2),
    end = map_chr(chr_parts, 3),
    row_label = paste0(gene_name, " (", chr, ": ", start, "-", end, ")")
  )

heatmap_matrix <- original_MPRA_data_filtered %>%
  column_to_rownames("index_offset") %>%  
  select(where(is.numeric)) %>%
  as.matrix()

plot <- pheatmap(
  heatmap_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = color_palette2_custom,
  na_col = "grey90",
  main = "PSI Heatmap",
  angle_col = 90
)
ggsave("/Volumes/broad_dawnccle/melange/figures_outputs/figS03/figS03_heatmap_kelly_specific_unbiasedMPRA.pdf", plot = plot, width = 15, height = 4, dpi = 300)




######## Figure 3B - Correlation between PSI in the KELLY-specific satmut data and the unbiased screen using the unmutated sequences ########


parent_psi_comparison <- original_MPRA_data %>%
  filter(condition %in% c("A375", "T47D", "HEK", "KELLY")) %>%
  left_join(parent_psi, by = c("index_offset", "condition")) %>%
  filter(!is.na(par_PSI))

# Pearson correlation
pearson_cor <- cor(parent_psi_comparison$PSI_original,
                   parent_psi_comparison$par_PSI,
                   method = "pearson", use = "complete.obs")

# Spearman correlation
spearman_cor <- cor(parent_psi_comparison$PSI_original,
                    parent_psi_comparison$par_PSI,
                    method = "spearman", use = "complete.obs")

# Print results
cat("Pearson correlation:", round(pearson_cor, 3), "\n")
cat("Spearman correlation:", round(spearman_cor, 3), "\n")

plot <- ggplot(parent_psi_comparison, aes(x = PSI_original, y = par_PSI, color = condition)) +
  geom_point(alpha = 1) +
  geom_smooth(method = "lm", color = "black") +
  xlim(0, 1) +
  ylim(0, 1) +
  labs(
    title = "Correlation between PSI and PSI_original",
    x = "Original PSI",
    y = "SatMut PSI"
  ) +
  theme_bw()
ggsave("/Volumes/broad_dawnccle/melange/figures_outputs/fig03/fig03_parent_psi_correlation.pdf", plot = plot, width = 4, height = 3, dpi = 300)

######## Figure S3B - Replicate correlation ########
psi_wide <- PSI_file_clean_to_PSI %>%
  mutate(full_name = paste0(index_offset, base, location)) %>%
  select(sample, full_name, PSI) %>%
  pivot_wider(names_from = sample, values_from = PSI)

psi_matrix <- psi_wide %>%
  select(-full_name) %>%
  as.data.frame() %>%
  as.matrix()

cor_matrix <- cor(psi_matrix, use = "pairwise.complete.obs")

plot <- pheatmap(
  cor_matrix,
  main = "Sample-to-Sample Correlation (PSI)",
  color = color_palette2_custom,
  border_color = NA
)
ggsave("/Volumes/broad_dawnccle/melange/figures_outputs/figS03/figS03_satmut_KELLY_specific_replicate_correlation.pdf", plot = plot, width = 4, height = 3, dpi = 300)



######## Figure 3C - Histogram of deltaPSI  #######

delta_PSI_df <- merged_PSI %>%
  mutate(deltaPSI = abs(par_PSI - PSI))
plot <- ggplot(delta_PSI_df %>% filter(deltaPSI >= 0 & deltaPSI <= 1),
               aes(x = deltaPSI)) +
  geom_histogram(binwidth = 0.05, fill = "blue", color = "black") +
  xlim(0, 1) +
  ylim(0,10000) + 
  scale_y_log10() +
  labs(x = "deltaPSI", y = "Count") +
  theme_minimal()

ggsave("/Volumes/broad_dawnccle/melange/figures_outputs/fig03/fig03_satmut_KELLY_specific_delta_PSI_histogram.pdf", plot = plot, width = 4, height = 3, dpi = 300)
cat("Proportion of deltaPSI < 0.1:", mean(delta_PSI_df$deltaPSI < 0.1), "\n")
cat("Proportion of deltaPSI > 0.5:", mean(delta_PSI_df$deltaPSI > 0.5), "\n")



######## Figure 3D - deltaPSI by position ########
delta_PSI_df <- merged_PSI %>%
  mutate(deltaPSI = PSI - par_PSI)

delta_PSI_df_kellyonly <- delta_PSI_df %>%
  filter(condition == "KELLY") %>%
  filter(par_PSI > 0.3) %>%
  filter(plot_position < 51 | location %in% c("upstream", "exon")) %>%
  filter(plot_position > -51 | location %in% c("downstream", "exon"))


# set global positions
delta_PSI_df_kellyonly_split <- delta_PSI_df_kellyonly %>%
  group_by(index_offset) %>%
  mutate(
    exon_cutoff = max(plot_position[location == "exon"], na.rm = TRUE) -10,
    plot_group = case_when(
      location == "upstream" ~ "Upstream + Exon(1–10)",
      location == "exon" & plot_position <= 10 ~ "Upstream + Exon(1–10)",
      location == "exon" & plot_position > exon_cutoff ~ "Exon(-10–end) + Downstream",
      location == "downstream" ~ "Exon(-10–end) + Downstream",
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup() %>%
  filter(!is.na(plot_group)) %>%
  mutate(
    global_position = case_when(
      location == "upstream"   ~ plot_position,
      location == "exon" & plot_position <= 10       ~ plot_position,
      location == "exon" & plot_position > exon_cutoff       ~ 10 + plot_position - exon_cutoff,
      location == "downstream" ~ plot_position + 20
    )
  ) %>%
  filter(global_position <= 70)


# make box and whisker plot
plot <- ggplot(delta_PSI_df_kellyonly_split, aes(x = global_position, y = deltaPSI)) +
  geom_boxplot(aes(group = global_position), outlier.size = 0.5, width = 0.8, na.rm = TRUE) +
  facet_wrap(~ plot_group, nrow = 2, scales = "fixed") +
  labs(
    title = expression("ΔPSI per Position (KELLY only)"),
    x = "Position",
    y = expression(Delta*"PSI")
  ) +
  theme_minimal()

# make bar graph data 
bar_data <- delta_PSI_df_kellyonly_split %>%
  group_by(global_position, plot_group) %>%
  summarise(n_variants = n(), .groups = "drop")

# bar plot
bar_plot <- ggplot(bar_data, aes(x = global_position, y = n_variants)) +
  geom_col(width = 0.8, fill = "gray40") +
  facet_wrap(~ plot_group, nrow = 2, scales = "fixed") +
  labs(y = "Variant Count", x = NULL) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(b = 0),
    strip.text = element_blank()  
  )
combined_plot <- bar_plot / plot + plot_layout(heights = c(1, 3))
ggsave("/Volumes/broad_dawnccle/melange/figures_outputs/fig03/fig03_satmut_KELLY_specific_delta_PSI_by_position.pdf", plot = combined_plot, width = 10, height = 6, dpi = 300)


# find z score statistics for overall distribution 
baseline_stats <- delta_PSI_df_kellyonly_split %>%
  summarise(
    baseline_mean = mean(deltaPSI, na.rm = TRUE),
    baseline_sd   = sd(deltaPSI, na.rm = TRUE)
  )

mu <- baseline_stats$baseline_mean
sigma <- baseline_stats$baseline_sd

# find adjusted p-value from z score test in each position
pos_summary <- delta_PSI_df_kellyonly_split %>%
  group_by(global_position) %>%
  summarise(mean_delta = mean(deltaPSI, na.rm = TRUE), n = n()) %>%
  ungroup() %>%
  mutate(
    z = (mean_delta - mu) / (sigma / sqrt(n)),  # one-sample z against baseline
    pval = 2 * pnorm(-abs(z))                   # two-tailed p-value
  ) %>%
  mutate(pval_adj = p.adjust(pval, method = "BH"))

print(pos_summary)





######## Figure 3G and H - SLC39A8 upstream + COL27A1 downstream dot plot ######
combos <- tibble::tribble(
  ~location,      ~index_offset_label,
  "upstream", "SLC39A8",
  "downstream", "COL27A1",
)

# Loop and plot each combination
for(i in seq_len(nrow(combos))) {
  loc <- combos$location[i]
  label <- combos$index_offset_label[i]
  
  df <- merged_PSI %>%
    filter(location == loc, index_offset_label == label) %>%
    group_by(plot_position, condition) %>%  # average across nucleotide
    summarise(PSI = mean(PSI, na.rm = TRUE), .groups = "drop")
  
  p <- ggplot(df, aes(plot_position, PSI, color = condition)) +
    geom_point(size = 1) +
    labs(
      title = paste(label, "in", loc),
      subtitle = paste0("location = ", loc, ", index_offset_label = ", label),
      x = "Position",
      y = "PSI"
    ) +
    theme_minimal()
  
  filename <- paste0("/Volumes/broad_dawnccle/melange/figures_outputs/fig03/fig03_", label, "_", loc, ".pdf") %>% str_replace_all(" ", "_")
  ggsave(filename, plot = p, width = 5, height = 3, dpi = 300)
}





######## Figure 3I - RBFOX1 expression barchart ########

cell_lines_names <- c(
  "786O", "769P", "A172", "A375", "ACHN", "CAL120", "OC314", "OC316", "COGN278",
  "COLO783", "DAOY", "DBTRG05MG", "EFO27", "GB1", "GAMG", "GI1", "HCC1428",
  "HCC38", "IGR37", "IPC298", "JHH6", "K562", "KELLY", "KMRC1", "KMRC20",
  "8MGBA", "MCF7", "MDAMB231", "MELHO", "MeWo", "PLCPRF5", "SF126", 
  "SKNAS", "SNU398", "SNU423", "SNU449", "SNUC4", "T47D", "TOV21G", "U251MG", "VMRCRCZ") # without HEK or MEL202


# GEX from Broad DepMap
all_gex <- read.csv("/Volumes/broad_dawnccle/melange/data/CCLE_expression.csv", sep = ",")
all_gex <- all_gex %>% 
  rename_with(~str_extract(., "^[^\\s]+")) 
colnames(all_gex)[1] <- "DepMap_ID"


gene_name <- "RBFOX1"
gene_filtered <- all_gex %>% select(DepMap_ID, starts_with(gene_name))

# Pivot longer to be DepMap_ID and "gene_name"
gene_filtered_long <- gene_filtered %>% pivot_longer(-DepMap_ID, names_to = "gene_name", values_to = "gex")



cellline_metadata <- read.csv("/Volumes/broad_dawnccle/melange/data/cellline_data_full_metadata.csv") %>% 
  select(DepMap_ID, StrippedName) %>% distinct()
# Rename V1 in gex based on cellline_metadata
gex_formatted <- gene_filtered_long %>% 
  left_join(cellline_metadata, by = "DepMap_ID") %>% 
  select(-DepMap_ID) %>% 
  filter(StrippedName %in% cell_lines_names) %>%
  mutate(StrippedName = factor(StrippedName, levels = sort(unique(StrippedName))))

# Create bar graph
plot <- ggplot(gex_formatted, aes(x = StrippedName, y = gex)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Cell Line", y = "GEX", title = "Gene Expression by Cell Line") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("/Volumes/broad_dawnccle/melange/figures_outputs/fig03/fig03_RBFOX1_by_cellline.pdf", plot = plot, width = 5, height = 3, dpi = 300)



######## Figure 3E + 3F, S3C-O - Prepare csv files to make lolliplots + z-score filtering test #########

merged_PSI_new <- merged_PSI %>%
  group_by(index_offset) %>%
  mutate(
    upstream_len   = max(plot_position[location == "upstream"], na.rm = TRUE),
    exon_len       = max(plot_position[location == "exon"], na.rm = TRUE),
    position = case_when(
      location == "upstream"   ~ plot_position,
      location == "exon"       ~ upstream_len + plot_position,
      location == "downstream" ~ upstream_len + exon_len + plot_position
    )
  ) %>%
  ungroup() %>%
  select(-location, -plot_position, -upstream_len, -exon_len) %>%
  group_by(index_offset) %>%
  mutate(
    position = position - min(position, na.rm = TRUE) + 1
  ) %>%
  ungroup() %>%
  mutate(delta_PSI = PSI - par_PSI) %>%
  select(index_offset, condition, position, nucleotide, delta_PSI) 

# Non-kelly cells
merged_PSI_nonKELLY <- merged_PSI_new %>%
  filter(condition != "KELLY") %>%
  mutate(condition = "nonKELLY")
write.csv(
  merged_PSI_nonKELLY,
  file = "/Volumes/broad_dawnccle/melange/data/satmut_KELLY_specific_nonKELLY_deltaPSI.csv",
  row.names = FALSE
  )

# Kelly cells 
merged_PSI_nonKELLY <- merged_PSI_new %>%
  filter(condition == "KELLY")
write.csv(
  merged_PSI_nonKELLY,
  file = "/Volumes/broad_dawnccle/melange/data/satmut_KELLY_specific_KELLY_deltaPSI.csv",
  row.names = FALSE
)


