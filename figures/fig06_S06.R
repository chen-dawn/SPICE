library(dplyr)
library(ggplot2)
library(vroom)
library(stringr)
library(tidyverse)
library(data.table)
library(tidyr)
library(patchwork)
library(ggrepel)
library(pheatmap)
library(ggseqlogo)
library(preprocessCore)
library(purrr)
library(RColorBrewer)
library(ggpubr)
library(scales)
library(ggrastr)
library(ragg)
library(Biostrings)



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
output_filepath_main <- "/Volumes/broad_dawnccle/melange/figures_outputs/fig06"
output_filepath_supp <- "/Volumes/broad_dawnccle/melange/figures_outputs/figS06"

library_47k_reference <- read_csv("/Volumes/broad_dawnccle/melange/data/20230130_twist_library_v3.csv")




######### Import data and build Figure 6C, S6D: significant PSI events and alt 3'SS volcano plots  #######

# Define file paths and output names
rmats_files <- list(
  "U2AF1S34F_PSI" = "/Volumes/broad_dawnccle/melange/data/rMATS_cancer_unbiased_library/PSI_U2AF1-WT_U2AF1-S34F_rMATS_Result_P.FDR.txt",
  "SF3B1_K700E_PSI" = "/Volumes/broad_dawnccle/melange/data/rMATS_cancer_unbiased_library/PSI_K562WT_K562K700E_rMATS_Result_P.FDR.txt",  
  "FUBP1_PSI" = "/Volumes/broad_dawnccle/melange/data/rMATS_cancer_unbiased_library/PSI_CH3-1_FUBP1_rMATS_Result_P.FDR.txt",
  "RBM5_PSI" = "/Volumes/broad_dawnccle/melange/data/rMATS_cancer_unbiased_library/PSI_CH3-1_RBM5_rMATS_Result_P.FDR.txt",
  "RBM10_PSI" = "/Volumes/broad_dawnccle/melange/data/rMATS_cancer_unbiased_library/PSI_CH3-1_RBM10_rMATS_Result_P.FDR.txt",
  "ZRSR2_PSI" = "/Volumes/broad_dawnccle/melange/data/rMATS_cancer_unbiased_library/PSI_CH3-1_ZRSR2_rMATS_Result_P.FDR.txt",
  "U2AF1_S34F_altSS" = "/Volumes/broad_dawnccle/melange/data/rMATS_cancer_unbiased_library/altSS_splicelib_U2AF1_WT_splicelib_U2AF1_S34F_rMATS_Result_P.FDR.txt",
  "SF3B1_K700E_altSS" = "/Volumes/broad_dawnccle/melange/data/rMATS_cancer_unbiased_library/altSS_K562WT_K562K700E_rMATS_Result_P.FDR.txt", 
  "FUBP1_altSS" = "/Volumes/broad_dawnccle/melange/data/rMATS_cancer_unbiased_library/altSS_CH3-1__FUBP1__rMATS_Result_P.FDR.txt",
  "RBM5_altSS" = "/Volumes/broad_dawnccle/melange/data/rMATS_cancer_unbiased_library/altSS_CH3-1__RBM5__rMATS_Result_P.FDR.txt",
  "RBM10_altSS" = "/Volumes/broad_dawnccle/melange/data/rMATS_cancer_unbiased_library/altSS_CH3-1__RBM10__rMATS_Result_P.FDR.txt",
  "ZRSR2_altSS" = "/Volumes/broad_dawnccle/melange/data/rMATS_cancer_unbiased_library/altSS_CH3-1__ZRSR2__rMATS_Result_P.FDR.txt"
)

PSI_names <- c("U2AF1S34F_PSI", "SF3B1_K700E_PSI", "FUBP1_PSI", "RBM5_PSI", "RBM10_PSI", "ZRSR2_PSI")


# Process into a single dataframe
all_data <- map2_dfr(names(rmats_files), rmats_files, function(name, filepath) {
  group <- if (name %in% PSI_names) "PSI" else "altSS"
  
  read_tsv(filepath, show_col_types = FALSE) %>%
    mutate(
      I1_sum = sapply(str_split(I1, ","), function(x) sum(as.numeric(x))),
      I2_sum = sapply(str_split(I2, ","), function(x) sum(as.numeric(x))),
      S1_sum = sapply(str_split(S1, ","), function(x) sum(as.numeric(x))),
      S2_sum = sapply(str_split(S2, ","), function(x) sum(as.numeric(x))),
      PSI1 = I1_sum / (I1_sum + S1_sum),
      PSI2 = I2_sum / (I2_sum + S2_sum),
      deltaPSI = PSI2 - PSI1,
      significant = ifelse(abs(deltaPSI) > 0.1 & FDR < 0.05, "yes", "no"),
      name = name,
      group = group
    )
})

# Assign plotting colors
all_data <- all_data %>%
  mutate(sig_color = case_when(
    significant == "no" ~ "not_significant",
    deltaPSI > 0        ~ "up",
    deltaPSI < 0        ~ "down"
  ))

color_values <- c(
  "up" = "#D83629",
  "down" = "#4575B4",
  "not_significant" = "gray"
)

# Plot definition
plot <- ggplot(all_data, aes(x = deltaPSI, y = -log10(FDR), color = sig_color)) +
  geom_point(alpha = 1, size = 0.6) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  scale_color_manual(values = color_values) +
  facet_wrap(~ group + name, ncol = 3, scales = "fixed") +
  scale_x_continuous(limits = c(-0.6, 0.6), oob = squish) +
  scale_y_continuous(limits = c(0, 20), oob = squish) +
  labs(
    x = "\u0394PSI",
    y = "-log10(FDR)",
    title = "Volcano Plots: PSI and altSS mutant cell lines",
    color = "Significant"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),     
    panel.grid.minor = element_blank(),     
    axis.ticks = element_line(color = "black"),  
    axis.ticks.length = unit(0.2, "cm"),    
    axis.text = element_text(color = "black"),
    axis.line = element_line(color = "black"),   # <-- add this
    strip.text = element_text(face = "bold", size = 10)
  ) 

# save to both the main and supp
outfile <- file.path(output_filepath_main, "fig06_cancer_unbiased_sequences_volcano.pdf")
ggsave(outfile, plot, width = 12, height = 8, dpi = 300)

outfile <- file.path(output_filepath_supp, "figS06_cancer_unbiased_sequences_volcano.pdf")
ggsave(outfile, plot, width = 12, height = 8, dpi = 300)



######### Figure 6D: number of significant PSI events ######

combined_counts <- all_data %>%
  filter(significant == "yes") %>%
  filter(group == "PSI") %>%
  count(name)
plot <- ggplot(combined_counts, aes(x = name, y = n)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "Folder",
    y = "Number of Events",
    title = "Alternative Splicing Events per Folder",
    fill = "Data Source"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
outfile <- file.path(output_filepath_main, "fig06_rmats_significant_events_PSI_bargraph.pdf")
ggsave(outfile, plot, width = 6, height = , dpi = 300)

######### Figure S6E: number of significant altSS events ######

combined_counts <- all_data %>%
  filter(significant == "yes") %>%
  filter(group == "altSS") %>%
  count(name)
plot <- ggplot(combined_counts, aes(x = name, y = n)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "Folder",
    y = "Number of Events",
    title = "Alternative Splicing Events per Folder",
    fill = "Data Source"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
outfile <- file.path(output_filepath_supp, "figS06_rmats_significant_events_altSS_bargraph.pdf")
ggsave(outfile, plot, width = 6, height = , dpi = 300)



######### Import designed cancer sequences rMATs results  #########
FDR_threshold <- 0.05
min_abs_delta_psi <- 0.1 


calculate_ratio <- function(I, S) {
  I_values <- as.numeric(unlist(strsplit(I, ",")))
  S_values <- as.numeric(unlist(strsplit(S, ",")))
  ratio <- I_values / (I_values + S_values)
  return(paste(round(ratio,3), collapse = ","))
}

calculate_average <- function(PSI){
  PSI_values <- as.numeric(unlist(strsplit(PSI, ",")))
  average <- mean(PSI_values)
  return(round(average, 3))
}

calculate_average_count_sum <- function(I, S){
  I_values <- as.numeric(unlist(strsplit(I, ",")))
  S_values <- as.numeric(unlist(strsplit(S, ",")))
  total_sum <- I_values + S_values
  average_count_sum <- mean(total_sum)
  return(round(average_count_sum, 0))
}

# loop through each file
sum_field <- function(x) {
  if (is.na(x) || x == "") return(NA_real_)
  nums <- suppressWarnings(as.numeric(str_split(as.character(x), ",")[[1]]))
  if (all(is.na(nums))) return(NA_real_)
  sum(nums, na.rm = TRUE)
}

label_from <- function(path) {
  parent <- basename(dirname(path))
  fname  <- basename(path)
  core   <- str_remove(fname, "\\.txt$")           
  core   <- str_remove(core, "_rMATS_Result_P\\.FDR$") 
  if (!is.na(parent) && nzchar(parent)) parent else core
}

files <- c("RBM5" = "/Volumes/broad_dawnccle/melange/data/rMATS_cancer_designed_library/PSI_designed_CH3-1_RBM5_rMATS_Result_P.FDR.txt",
           "RBM10" = "/Volumes/broad_dawnccle/melange/data/rMATS_cancer_designed_library/PSI_designed_CH3-1_RBM10_rMATS_Result_P.FDR.txt")

read_one <- function(f) {
  read_tsv(
    f, show_col_types = FALSE,
    col_types = cols(
      .default = col_guess(),
      I1 = col_character(), I2 = col_character(),
      S1 = col_character(), S2 = col_character()
    )
  ) %>%
    rowwise() %>%
    mutate(
      I1_sum = sum_field(I1),
      I2_sum = sum_field(I2),
      S1_sum = sum_field(S1),
      S2_sum = sum_field(S2),
      PSI1   = ifelse((I1_sum + S1_sum) > 0, I1_sum / (I1_sum + S1_sum), NA_real_),
      PSI2   = ifelse((I2_sum + S2_sum) > 0, I2_sum / (I2_sum + S2_sum), NA_real_),
      deltaPSI = PSI2 - PSI1,   # match axis label: PSI2 - PSI1
      significant = ifelse(abs(deltaPSI) > min_abs_delta_psi & FDR < FDR_threshold, "yes", "no"),
      facet_label = names(files)[match(f, files)],
      file = f
    ) %>%
    ungroup()
}

all_df <- map_dfr(files, read_one) %>%
  separate(ExonID, into = c("index", "offset", "design"), sep = "__", remove = FALSE) 

######### Figure S6F: Make volcano plots of significant events for designed sequence #######
# Assign plotting colors and filter out non-parent sequences
sig_df <- all_df  %>% 
  filter(design != "R2design0") %>%
  mutate(sig_color = case_when(
    significant == "no" ~ "not_significant",
    deltaPSI > 0        ~ "up",
    deltaPSI < 0        ~ "down"
  ))

color_values <- c(
  "up" = "#D83629",
  "down" = "#4575B4",
  "not_significant" = "gray"
)

p <- ggplot(sig_df, aes(x = deltaPSI, y = -log10(FDR)), color = sig_color) +
  geom_point(alpha = 1, size = 0.6) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  facet_wrap(~ facet_label, scales = "fixed") +
  scale_x_continuous(limits = c(-0.4, 0.4), oob = squish) +
  scale_y_continuous(limits = c(0, 10), oob = squish) +
  scale_color_manual(values = color_values) +
  labs(
    title = "Volcano plots across all pairs/files",
    x = "ΔPSI (PSI2 − PSI1)",
    y = "−log10(FDR)"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),     # remove major gridlines
    panel.grid.minor = element_blank(),     # remove minor gridlines
    axis.ticks = element_line(color = "black"),  
    axis.ticks.length = unit(0.2, "cm"),    
    axis.text = element_text(color = "black"),
    axis.line = element_line(color = "black"),   # <-- add this
    strip.text = element_text(face = "bold", size = 10)
  ) 


outfile <- file.path(output_filepath_supp, "figS06_designed_RBM5_RMB10_volcano.pdf")
ggsave(outfile, plot = p, width = 8, height = 4, dpi = 300)  


######### Figure 6E: Make heatmaps of significant designed sequences ###########
sig_df <- all_df %>%
  mutate(is_sig = significant == "yes" & design != "R2design0") %>%
  filter(is_sig) %>%
  ungroup()

outfile_names <- c(
  "RBM5" = "fig06_RBM5_designed_sig_hits_heatmap.pdf",
  "RBM10" = "fig06_RBM10_designed_sig_hits_heatmap.pdf"
)

sig_df %>%
  group_split(facet_label) %>%
  walk(function(df) {
    this_pair <- unique(df$facet_label)
    
    heat_mat <- df %>%
      select(ExonID, PSI1, PSI2) %>%
      column_to_rownames("ExonID") %>%
      as.matrix()
    
    outfile <- file.path(
      output_filepath_main,
      outfile_names[[this_pair]])
    
    pdf(outfile, width = 15, height = 5)
    pheatmap(
      heat_mat,
      scale = "none",
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      clustering_method = "complete",
      color = colorRampPalette(color_palette2_custom)(100),
      main = paste("PSI values –", this_pair)
    )
    dev.off()
    
    message("Saved: ", outfile)
  })

######### Figure 6I: Make heatmaps for some of the example positive hits hits to compare to the parent sequences ##########
sig_df <- all_df %>%
  mutate(is_sig = significant == "yes") %>%
  ungroup() %>%
  mutate(
    category = if_else(str_detect(design, "_"),
                       str_replace(design, "_[^_]+$", ""),   
                       "parent"),
    R1design = if_else(str_detect(design, "_"),
                       str_extract(design, "[^_]+$"),        
                       design)
  )


make_heatmap <- function(index_value, cats, outfile) {
  df <- sig_df %>%
    filter(index == !!index_value) %>%   
    mutate(significant = tidyr::replace_na(significant, FALSE)) %>%
    filter(is_sig | category == "parent") %>%    
    pivot_longer(
      cols = c("PSI1", "PSI2"),
      names_to = "condition",
      values_to = "PSI"
    ) %>%
    mutate(
      trailing_num = as.numeric(str_extract(design, "\\d+$")),
      design = factor(design,
                      levels = unique(design[order(-trailing_num)]),
                      ordered = TRUE)
    )

    p <- ggplot(df, aes(x = condition, y = design, fill = PSI)) +
    geom_tile(color = "white") +
    scale_fill_gradientn(colors = color_palette2_custom, limits = c(0, 1)) +
    theme_minimal(base_size = 12) +
    labs(
      title = paste0("PSI heatmap – ", paste(cats, collapse = ", ")),
      x = "Condition",
      y = "Design",
      fill = "PSI"
    )
  
  ggsave(
    filename = file.path(output_filepath_main, outfile),
    plot = p,
    width = 5,
    height = 4,
    dpi = 300
  )
  
  message("Saved: ", outfile)
}

jobs <- tibble::tibble(
  index = c(
    "ENSG00000020922.13;MRE11;chr11-94435831-94435899-94429910-94429986-94437176-94437235",
    "ENSG00000169410.10;PTPN9;chr15-75479847-75479914-75473688-75473767-75490207-75490301",
    "ENSG00000172071.15;EIF2AK3;chr2-88558916-88558979-88556740-88557936-88570873-88571041",
    "ENSG00000153982.11;GDPD1;chr17-59248739-59248785-59245413-59245549-59267040-59267174"
  ),
  cats = list(
    c("RBM5"),
    c("RBM5"),
    c("RBM10"),
    c("RBM10")
  ),
  outfile = c(
    "fig06_cancer_specific_example_MRE11.pdf",
    "fig06_cancer_specific_example_PTPN9.pdf",
    "fig06_cancer_specific_example_EIF2AK3.pdf",
    "fig06_cancer_specific_example_GDPD1.pdf"
  )
)

for (i in seq_len(nrow(jobs))) {
  make_heatmap(
    index_value = jobs$index[i],
    cats         = jobs$cats[[i]],
    outfile      = jobs$outfile[i]
  )
}



