library(dplyr)
library(tidyr)
library(vroom)
library(stringr)
library(tidyverse)
library(data.table)
library(ggplot2)
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

safe_substr <- function(seq, start, stop) {
  ifelse(nchar(seq) >= stop, substr(seq, start, stop), NA)
}
output_filepath <- "/Volumes/broad_dawnccle/melange/figures_outputs/fig04"

whitelist_genes <- c("MEP1A", "MFSD1", "BTN2A3P", "UBE2G1", "CDK15", "CTPS1", "KAT5", "ADCY5", "PRPF8", "RNF5", "AL928654.7", 
                     "RP11-654O1.1", "VEZT", "SNX5", "PPP2CB", "TIA1", "KTN1", "INTS7", "TBCK", "NCAM1", "DCHS2", "COL28A1", "CCDC26")

######### Data Processing: Import predicted satmut data #######
library_47k_reference <- read_csv("/Volumes/broad_dawnccle/melange/data/20230130_twist_library_v3.csv")
predicted_PSI_file <- read_csv("/Volumes/broad_dawnccle/melange/data/satmut_constitutive_predicted.csv")

predicted_PSI_df <- predicted_PSI_file %>%
  left_join(library_47k_reference, by = c("ID")) %>%
  filter(sample == "gene_barcode") %>%
  mutate(gene = str_trim(str_split_fixed(ID, ";", 3)[, 2])) %>%
  filter(gene %in% whitelist_genes)


######### Data Processing: Import unbiased library data #######

# we will use these results for parent values if there is no parent data in the observed sat mut data
unbiased_library_PSI <- read_csv("/Volumes/broad_dawnccle/melange/data/250527_FINAL_PSI_TABLE_by_condition.csv")
unbiased_library_PSI_filtered <- unbiased_library_PSI %>%
  separate(index_offset, into = c("index", "offset"), sep = "__") %>%
  filter(condition %in% c("8MGBA","HCC1428","KELLY","KMRC20","HCC38","JHH6","MCF7","T47D")) %>%
  filter(index %in% parent_reference$ExonID) %>%
  filter(offset == "0:0:0") %>%
  mutate(index_offset = paste(index, offset, sep = "__"))


######### Data Processing: Import observed satmut data ########

PSI_file <- read.csv("/Volumes/broad_dawnccle/melange/data/satmut_constitutive_count_table_normalized.csv")
parent_reference <- read.csv("/Volumes/broad_dawnccle/melange/data/satmut_constitutive_index_switch_reference.csv")

PSI_file_clean <- PSI_file %>% 
  filter(offset %in% c("0:0:0", "0")) %>%  
  select(-filename) %>%
  filter(index %in% whitelist_genes) %>%
  filter(condition != "HCC1428") %>%
  mutate(condition = str_replace(condition, "Kelly", "KELLY")) %>%
  filter(mode %in% c("INCLUDED", "SKIPPED")) %>% 
  group_by(sample, condition, index, mode, offset, location, base) %>%
  summarise(count = sum(count))

PSI_file_clean_to_PSI <- PSI_file_clean %>% 
  group_by(sample, condition, index, mode, location, base) %>%
  summarise(count = sum(count)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = mode, values_from = count, values_fill = 0) %>%
  mutate(total_counts = INCLUDED + SKIPPED) %>%
  filter(total_counts >= 30) %>%
  mutate(PSI = INCLUDED / (INCLUDED + SKIPPED)) %>%
  mutate(position = str_extract(base, "\\d+")) %>% 
  mutate(position = as.integer(position)) %>%
  mutate(nucleotide = str_extract(base, "[A-Z]$")) 

PSI_file_clean_to_PSI <- PSI_file_clean_to_PSI %>%
  group_by(index, condition, position, nucleotide, location) %>% 
  summarise(PSI = mean(PSI), .groups = "drop")

# Find parent sequences' PSI values in the sat mut data
parent_psi <- PSI_file_clean_to_PSI %>%
  mutate(location = trimws(as.character(location))) %>%
  filter(is.na(location) | location == "") %>%
  group_by(condition, index) %>%
  summarise(par_PSI = mean(PSI, na.rm = TRUE), .groups = "drop")

# Find parent sequences' PSI values in the unbiased library
unbiased_parents <- unbiased_library_PSI %>%
  separate(index_offset, into = c("index", "offset"), sep = "__", remove = FALSE) %>%
  filter(condition %in% c("8MGBA","HCC1428","KELLY","KMRC20","HCC38","JHH6","MCF7","T47D"),
         index %in% parent_reference$ExonID,
         offset == "0:0:0") %>%
  mutate(index = as.character(index)) %>%
  group_by(condition, index) %>%
  summarise(unbiased_par_PSI = mean(PSI, na.rm = TRUE), .groups = "drop") %>%
  mutate(index = str_trim(str_split_fixed(index, ";", 3)[, 2]))

# for sequences that don't have a parent PSI in the sat mut data, use the PSI values from the unbiased library
parent_psi_final <- parent_psi %>%
  full_join(unbiased_parents, by = c("condition", "index")) %>%
  mutate(par_PSI = coalesce(par_PSI, unbiased_par_PSI),
         par_PSI = ifelse(!is.na(par_PSI), par_PSI + 0.0001, NA_real_)) %>%   # add epsilon only if we have a value
  select(condition, index, par_PSI)


non_parent_psi <- PSI_file_clean_to_PSI %>% 
  filter(!is.na(location))

observed_PSI_df <- merge(non_parent_psi, parent_psi_final, by = c("index", "condition"), all.x = T) %>% 
  mutate(location = factor(location, levels = c("upstream", "exon", "downstream", "extra3"))) %>%
  mutate(plot_position = ifelse(location == "upstream", -position, position)) %>%
  left_join(parent_reference, by = c("index" = "gene")) %>%
  select(-index) 

# per-location ranks
observed_PSI_df <- observed_PSI_df %>%
  mutate(position = as.numeric(position)) %>%
  group_by(location) %>%
  mutate(
    rank_up = if_else(location == "upstream",   dplyr::dense_rank(-position), NA_integer_),
    rank_ex = if_else(location == "exon",       dplyr::dense_rank( position), NA_integer_),
    rank_dn = if_else(location == "downstream", dplyr::dense_rank( position), NA_integer_)
  ) %>%
  ungroup()



observed_PSI_df <- observed_PSI_df %>%
  left_join(library_47k_reference, by = c("ExonID" = "ID")) %>%
  mutate(
    rank_up = as.integer(rank_up),
    rank_ex = as.integer(rank_ex),
    rank_dn = as.integer(rank_dn),
    len_up  = nchar(upstreamIntronSeq),   
    len_ex  = nchar(skippedExonSeq),      
    overall_index = case_when(
      location == "upstream"   ~ rank_up + len_up - 60,
      location == "exon"       ~ rank_ex + len_up,
      location == "downstream" ~ rank_dn + len_up + len_ex,
      TRUE ~ NA_integer_
    )
  ) %>%
  select(-rank_up, -rank_ex, -rank_dn) 


######### Figure 4I: satmut constitutive observed vs predicted plot correlation ########

predicted_df <- predicted_PSI_df %>%
  group_by(ID, mutation_pos, mutation_base, Celltype) %>%
  mutate(average_PSI = mean(Predicted_PSI_val))
observed_df <- observed_PSI_df %>%
  group_by(ExonID, overall_index, nucleotide, condition) %>%
  mutate(average_PSI = mean(PSI))

# Rename observed keys to match predicted keys for an easy join
observed_df2 <- observed_df %>%
  dplyr::rename(
    ID            = ExonID,
    mutation_pos  = overall_index,
    mutation_base = nucleotide,
    Celltype = condition
  )

key_cols <- c("ID", "mutation_pos", "mutation_base", "Celltype")
pred_val_col <- "average_PSI"
obvs_val_col <- "average_PSI"

# Join and prep plotting df
plot_df <- predicted_df %>%
  inner_join(observed_df2, by = key_cols, suffix = c(".pred", ".obvs")) %>%
  transmute(
    ID, mutation_pos, mutation_base,
    predicted = .data[[paste0(pred_val_col, ".pred")]],
    observed  = .data[[paste0(obvs_val_col,  ".obvs")]]
  ) %>%
  drop_na(predicted, observed) %>%
  distinct()
# Correlation
pearson_r <- cor(plot_df$predicted, plot_df$observed, use = "complete.obs", method = "pearson")
cat(sprintf("Pearson r = %.4f (n = %d)\n", pearson_r, nrow(plot_df)))
pearson_r <- cor(plot_df$predicted, plot_df$observed, use = "complete.obs", method = "spearman")
cat(sprintf("Spearman r = %.4f (n = %d)\n", pearson_r, nrow(plot_df)))




p <- ggplot(plot_df, aes(x = observed, y = predicted)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_bin2d(binwidth = c(0.02, 0.02)) +  # Adjust bin size as needed
  scale_fill_gradientn(
    colors = color_palette2_custom,
    trans = "log10",
    name = "Log10 Count",
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", math_format(10^.x))
  ) +
  labs(title = "Predicted vs observed PSI (Density Boxes)") +
  theme_classic()
outfile <- file.path(output_filepath, "fig04_predicted_vs_observed_density.pdf")
ggsave(outfile, p, width = 5.5, height = 5.0, dpi = 300)


######### Figure 4J,K: Make csv files for lolliplots ##########


# predicted effect sizes #
effect_sizes_predicted <- predicted_PSI_df %>%
  mutate(index_offset = paste(ID, offset, sep = "__")) %>%
  left_join(unbiased_library_PSI_filtered, by = c("Celltype" = "condition", "index_offset" = "index_offset")) %>%
  mutate(deltaPSI = Predicted_PSI_val - PSI) %>%
  select(index_offset, Celltype, mutation_pos, mutation_base, deltaPSI) %>%
  group_by(index_offset, mutation_pos, mutation_base) %>%
  summarise(delta_PSI = mean(deltaPSI, na.rm = TRUE), .groups = "drop") %>%
  arrange(index_offset, mutation_pos, mutation_base) %>%
  mutate(Celltype = "predicted")
write_csv(effect_sizes_predicted, file.path("/Volumes/broad_dawnccle/melange/data/satmut_constitutive_effect_sizes_predicted.csv"))


# observed effect sizes #
effect_sizes_observed <- observed_PSI_df %>%
  mutate(index_offset = paste(ExonID, "0:0:0", sep = "__")) %>%
  mutate(deltaPSI = PSI - par_PSI) %>%
  select(index_offset, condition, overall_index, nucleotide, deltaPSI) %>%
  group_by(index_offset, overall_index, nucleotide) %>%
  summarise(delta_PSI = mean(deltaPSI, na.rm = TRUE), .groups = "drop") %>%
  arrange(index_offset, overall_index, nucleotide) %>%
  mutate(Celltype = "observed")


names(effect_sizes_observed)[names(effect_sizes_observed) == "nucleotide"] <- "mutation_base"
names(effect_sizes_observed)[names(effect_sizes_observed) == "overall_index"] <- "mutation_pos"

write_csv(effect_sizes_observed, file.path("/Volumes/broad_dawnccle/melange/data/satmut_constitutive_effect_sizes_observed.csv"))


# parent reference #
unbiased_library_sequences <- read_csv("/Volumes/broad_dawnccle/melange/data/20230130_twist_library_v3.csv")
parent_reference_updated <- parent_reference %>%
  left_join(unbiased_library_sequences, by = c("ExonID" = "ID")) %>%
  mutate(index_offset = paste(ExonID, "0:0:0", sep = "__")) %>%
  select(index_offset, librarySequence) 
write_csv(parent_reference_updated, file.path("/Volumes/broad_dawnccle/melange/data/satmut_constitutive_lolliplot_reference.csv"))

######### make bar graph of number of splice sites accurately predicted #########
data <- data.frame(
  category = c("3' measured", "5' measured", "5' predicted", "3' predicted"),
  fraction = c(1, 1, 1, 11/14)
)

# Create the plot
p <- ggplot(data, aes(x = category, y = fraction)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  ylim(0, 1.1) +
  ylab("Fraction sensitive to mutation") +
  xlab("") +
  ggtitle("Splice Site Sensitivity to Mutation") +
  theme(
    axis.ticks = element_line(color = "black"),     # Enable tick lines
    axis.ticks.length = unit(0.25, "cm"),           # Length of ticks
    axis.line = element_line(color = "black")       # Add axis lines
  ) +
  scale_y_continuous(breaks = seq(0, 1.0, by = 0.2)) # Add y-axis tick marks

# Save the plot as a PDF
ggsave("/Volumes/broad_dawnccle/melange/figures_outputs/fig04/fig04_satmut_constitutive_splicesite_accuracy_barplot.pdf", plot = p, width = 6, height = 4, units = "in")



