library(tidyverse)
library(vroom)
library(data.table)
library(Biostrings)
library(ggpointdensity)
library(pheatmap)
library(gridExtra)
library(ggrepel)

combined_psi <- fread("~/Dropbox (Harvard University)/02Splicing/latest/Nova241106_PSI_combined_output_indiv.tsv")
# combined_psi <- fread("~/Dropbox (Harvard University)/02Splicing/latest/Nova240826_3ss_combined_output_indiv.tsv")
calculate_ratio <- function(I, S) {
  I_values <- as.numeric(unlist(strsplit(I, ",")))
  S_values <- as.numeric(unlist(strsplit(S, ",")))
  ratio <- I_values / (I_values + S_values)
  # Check if any of the values are NaN
  if (any(is.nan(ratio))) {
    ratio[is.nan(ratio)] <- NA
  }
  return(paste(round(ratio,3), collapse = ","))
}

calculate_average <- function(PSI){
  PSI_values <- as.numeric(unlist(strsplit(PSI, ",")))
  average <- mean(PSI_values, na.rm = T)
  return(round(average, 3))
}

get_max_PSI <- function(PSI){
  PSI_values <- as.numeric(unlist(strsplit(PSI, ",")))
  max_PSI <- max(PSI_values, na.rm = T)
  return(round(max_PSI, 3))
}

get_min_PSI <- function(PSI){
  PSI_values <- as.numeric(unlist(strsplit(PSI, ",")))
  min_PSI <- min(PSI_values, na.rm = T)
  return(round(min_PSI, 3))
}

calculate_average_count_sum <- function(I, S){
  I_values <- as.numeric(unlist(strsplit(I, ",")))
  S_values <- as.numeric(unlist(strsplit(S, ",")))
  total_sum <- I_values + S_values
  average_count_sum <- mean(total_sum)
  return(round(average_count_sum, 0))
}

# Apply the function to the data frame
combined_psi <- combined_psi %>%
  # Filter out folders that contain the word "A172" or Kelly or T47D
  filter(!grepl("Kelly", Folder)) %>%
  filter(!grepl("T47D", Folder)) %>%
  filter(!grepl("A172", Folder)) %>%
  filter(!grepl("HCC38", Folder)) %>%
  filter(!grepl("IPC298", Folder)) %>%
  filter(!grepl("KMRC20", Folder)) %>%
  filter(!grepl("KMRC1", Folder)) %>%
  filter(!grepl("HEK", Folder)) %>%
  mutate(
    PSI1 = mapply(calculate_ratio, I1, S1),
    PSI2 = mapply(calculate_ratio, I2, S2)
  ) %>% 
  mutate(
    PSI1_average = mapply(calculate_average, PSI1),
    PSI2_average = mapply(calculate_average, PSI2)
  ) %>%
  mutate(PSI_diff = PSI1_average - PSI2_average) %>% 
  mutate(
    count_sum_average1 = mapply(calculate_average_count_sum, I1, S1),
    count_sum_average2 = mapply(calculate_average_count_sum, I2, S2)
  ) %>% mutate(PSI_ratio = PSI1_average / PSI2_average) %>% 
  mutate(PSI_reverse_ratio = (1-PSI1_average)/(1-PSI2_average)) %>% 
  mutate(max_PSI1 = mapply(get_max_PSI, PSI1), max_PSI2 = mapply(get_max_PSI, PSI2)) %>% 
  mutate(min_PSI1 = mapply(get_min_PSI, PSI1), min_PSI2 = mapply(get_min_PSI, PSI2))


combined_psi_filtered <- combined_psi %>% 
  filter(count_sum_average1 >30) %>% 
  filter(count_sum_average2 > 30) %>% 
  mutate(log2_PSI_ratio = log2(PSI_ratio), log2_PSI_reverse_ratio = log2(PSI_reverse_ratio))
# write_csv(combined_psi_filtered, "combined_psi_filtered.csv")


sgRNA_only <- combined_psi_filtered %>% 
  filter(grepl("MUT-sgCh3-1_MUT-sg", Folder))



# test2 <- combined_psi %>% filter(ExonID %in% c("ENSG00000138821.14;SLC39A8;chr4-102285948-102285973-102267871-102268079-102304316-102304481", "ENSG00000135365.16;PHF21A;chr11-45946075-45946098-45938156-45938312-45948885-45948946"))

top_seq <- sgRNA_only %>% filter(abs(log2_PSI_ratio) > 2 | abs(log2_PSI_reverse_ratio) > 2)
num_obs_per_seq <- top_seq %>% group_by(ExonID) %>% summarise(num_obs = n())
top_seq <- top_seq %>% left_join(num_obs_per_seq, by = "ExonID") %>% filter(num_obs > 1)


# Only look at MUT-sgCh3-1_MUT-sgRUBP1
combined_psi_RUBP1 <- combined_psi_filtered %>% 
  filter(grepl("MUT-sgCh3-1_MUT-sgRUBP1", Folder)) %>% 
  filter(count_sum_average1 >30) %>% 
  filter(count_sum_average2 > 30) %>% 
  mutate(log2_PSI_ratio = log2(PSI_ratio), log2_PSI_reverse_ratio = log2(PSI_reverse_ratio))

# Only look at MUT-sgCh3-1_MUT-sgRBM5
combined_psi_RBM5 <- combined_psi_filtered %>% 
  filter(grepl("MUT-sgCh3-1_MUT-sgRBM5", Folder)) %>% 
  filter(count_sum_average1 >30) %>% 
  filter(count_sum_average2 > 30) %>% 
  mutate(log2_PSI_ratio = log2(PSI_ratio), log2_PSI_reverse_ratio = log2(PSI_reverse_ratio))

# Only look at MUT-sgCh3-1_MUT-sgRBM10
combined_psi_RBM10 <- combined_psi_filtered %>% 
  filter(grepl("MUT-sgCh3-1_MUT-sgRBM10", Folder)) %>% 
  filter(count_sum_average1 >30) %>% 
  filter(count_sum_average2 > 30) %>% 
  mutate(log2_PSI_ratio = log2(PSI_ratio), log2_PSI_reverse_ratio = log2(PSI_reverse_ratio))

# Only look at MUT-plx-LacZ_MUT-sgRUBP1
combined_psi_RUBP1_plx <- combined_psi_filtered %>% 
  filter(grepl("MUT-plx-LacZ_MUT-sgRUBP1", Folder)) %>% 
  filter(count_sum_average1 >30) %>% 
  filter(count_sum_average2 > 30) %>% 
  mutate(log2_PSI_ratio = log2(PSI_ratio), log2_PSI_reverse_ratio = log2(PSI_reverse_ratio))

# Only look at MUT-plx-LacZ_MUT-sgRBM5
combined_psi_RBM5_plx <- combined_psi_filtered %>% 
  filter(grepl("MUT-plx-LacZ_MUT-sgRBM5", Folder)) %>% 
  filter(count_sum_average1 >30) %>% 
  filter(count_sum_average2 > 30) %>% 
  mutate(log2_PSI_ratio = log2(PSI_ratio), log2_PSI_reverse_ratio = log2(PSI_reverse_ratio))

# Only look at MUT-plx-LacZ_MUT-sgRBM10
combined_psi_RBM10_plx <- combined_psi_filtered %>% 
  filter(grepl("MUT-plx-LacZ_MUT-sgRBM10", Folder)) %>% 
  filter(count_sum_average1 >30) %>% 
  filter(count_sum_average2 > 30) %>% 
  mutate(log2_PSI_ratio = log2(PSI_ratio), log2_PSI_reverse_ratio = log2(PSI_reverse_ratio))

# For RBM10, what are the overlap numbers? 
num_unique_RBM10 <- length(unique(combined_psi_RBM10$ExonID))
num_unique_RBM10_plx <- length(unique(combined_psi_RBM10_plx$ExonID))
combined_psi_RBM10_overlap <- intersect(combined_psi_RBM10$ExonID, combined_psi_RBM10_plx$ExonID)
paste("Num unique in RBM10", num_unique_RBM10)
paste("Num unique in RBM10 plx", num_unique_RBM10_plx)
paste("Num overlap in RBM10", length(combined_psi_RBM10_overlap))

# Do the same for RBM5
num_unique_RBM5 <- length(unique(combined_psi_RBM5$ExonID))
num_unique_RBM5_plx <- length(unique(combined_psi_RBM5_plx$ExonID))
combined_psi_RBM5_overlap <- intersect(combined_psi_RBM5$ExonID, combined_psi_RBM5_plx$ExonID)
paste("Num unique in RBM5", num_unique_RBM5)
paste("Num unique in RBM5 plx", num_unique_RBM5_plx)
paste("Num overlap in RBM5", length(combined_psi_RBM5_overlap))

# Do the same for RUBP1
num_unique_RUBP1 <- length(unique(combined_psi_RUBP1$ExonID))
num_unique_RUBP1_plx <- length(unique(combined_psi_RUBP1_plx$ExonID))
combined_psi_RUBP1_overlap <- intersect(combined_psi_RUBP1$ExonID, combined_psi_RUBP1_plx$ExonID)
paste("Num unique in RUBP1", num_unique_RUBP1)
paste("Num unique in RUBP1 plx", num_unique_RUBP1_plx)
paste("Num overlap in RUBP1", length(combined_psi_RUBP1_overlap))


######## Look at just one file and plot a volcano plot.
RBP1_all_df <- fread("/Volumes/broad_dawnccle/processed_data/missplicing_processed_df/novaseq241106/pairadise_indiv_PSI/splicelib_hek_no_dox_rbp1_splicelib_hek_dox_rbp1/splicelib_hek_no_dox_rbp1_splicelib_hek_dox_rbp1_rMATS_Result_P.FDR.txt")
RBP1_all_df <- RBP1_all_df %>% 
  mutate(
    PSI1 = mapply(calculate_ratio, I1, S1),
    PSI2 = mapply(calculate_ratio, I2, S2)
  ) %>% 
  mutate(
    PSI1_average = mapply(calculate_average, PSI1),
    PSI2_average = mapply(calculate_average, PSI2)
  ) %>%
  mutate(PSI_diff = PSI1_average - PSI2_average) %>% 
  mutate(
    count_sum_average1 = mapply(calculate_average_count_sum, I1, S1),
    count_sum_average2 = mapply(calculate_average_count_sum, I2, S2)
  ) %>% mutate(PSI_ratio = PSI1_average / PSI2_average) %>% 
  mutate(PSI_reverse_ratio = (1-PSI1_average)/(1-PSI2_average)) %>% 
  mutate(max_PSI1 = mapply(get_max_PSI, PSI1), max_PSI2 = mapply(get_max_PSI, PSI2)) %>% 
  mutate(min_PSI1 = mapply(get_min_PSI, PSI1), min_PSI2 = mapply(get_min_PSI, PSI2)) %>% 
  mutate(gene_name = map_chr(ExonID, ~str_split(.x, ";")[[1]][2])) %>% 
  filter(count_sum_average1 > 30 & count_sum_average2 > 30)


# Add a column to identify the top genes based on a condition
RBP1_all_df <- RBP1_all_df %>%
  # Set the one where p-value = 0 to be 1e-10.
  mutate(FDR = ifelse(FDR == 0, 1e-10, FDR)) %>%
  # mutate(top_genes = ifelse((log2(PSI_ratio) > 2 | log2(PSI_ratio) < -2) & (-log10(FDR) > -log10(0.01)), gene_name, NA)) %>%
  mutate(top_genes = ifelse(abs(PSI_diff) > 0.2 & (-log10(FDR) > -log10(0.01)), gene_name, NA)) %>%
  mutate(is_top_gene = ifelse(!is.na(top_genes), "top_gene", "other"))

# Plot with ggplot2 and ggrepel
ggplot(RBP1_all_df, aes(x = PSI_diff, y = -log10(FDR), color = is_top_gene)) + 
  geom_point() + 
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed") + 
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") + 
  geom_text_repel(aes(label = top_genes), na.rm = TRUE) + 
  scale_color_manual(values = c("top_gene" = "red", "other" = "grey")) +
  theme_bw() + 
  ggtitle("Control sgRNA vs RBM5 KO")


#################### All RBP ###################
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)

# Define the RBPs and file base path
rbps <- paste0("rbp", 1:12)
base_path <- "/Volumes/broad_dawnccle/processed_data/missplicing_processed_df/novaseq241106/pairadise_indiv_PSI"

# Define the save path base
save_path_base <- "~/Dropbox (Harvard University)/02Splicing/SplicingManuscript/figure_outputs"

# Loop through each RBP
for (rbp in rbps) {
  # Construct the file path dynamically
  file_path <- file.path(base_path, paste0("splicelib_hek_no_dox_", rbp, "_splicelib_hek_dox_", rbp),
                         paste0("splicelib_hek_no_dox_", rbp, "_splicelib_hek_dox_", rbp, "_rMATS_Result_P.FDR.txt"))
  
  # Read and process the data
  all_df <- fread(file_path) %>%
    mutate(
      PSI1 = mapply(calculate_ratio, I1, S1),
      PSI2 = mapply(calculate_ratio, I2, S2),
      PSI1_average = mapply(calculate_average, PSI1),
      PSI2_average = mapply(calculate_average, PSI2),
      PSI_diff = PSI1_average - PSI2_average,
      count_sum_average1 = mapply(calculate_average_count_sum, I1, S1),
      count_sum_average2 = mapply(calculate_average_count_sum, I2, S2),
      PSI_ratio = PSI1_average / PSI2_average,
      PSI_reverse_ratio = (1 - PSI1_average) / (1 - PSI2_average),
      max_PSI1 = mapply(get_max_PSI, PSI1),
      max_PSI2 = mapply(get_max_PSI, PSI2),
      min_PSI1 = mapply(get_min_PSI, PSI1),
      min_PSI2 = mapply(get_min_PSI, PSI2),
      gene_name = map_chr(ExonID, ~str_split(.x, ";")[[1]][2])
    ) %>%
    filter(count_sum_average1 > 30 & count_sum_average2 > 30) %>%
    mutate(
      FDR = ifelse(FDR == 0, 1e-10, FDR),
      top_genes = ifelse(abs(PSI_diff) > 0.2 & (-log10(FDR) > -log10(0.01)), gene_name, NA),
      is_top_gene = ifelse(!is.na(top_genes), "top_gene", "other")
    )
  
  # Generate the plot
  plot <- ggplot(all_df, aes(x = log2(PSI_ratio), y = -log10(FDR), color = is_top_gene)) +
    geom_point() +
    geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
    geom_text_repel(aes(label = top_genes), na.rm = TRUE) +
    scale_color_manual(values = c("top_gene" = "red", "other" = "grey")) +
    theme_bw() +
    xlab("log2(PSI ratio)") +
    ggtitle(paste("No Dox vs Dox for", toupper(rbp)))
  
  # Construct the save path dynamically
  save_path <- file.path(save_path_base, paste0(toupper(rbp), "_volcano_plot_ratio.png"))
  
  # Save the plot
  ggsave(save_path, plot, width = 10, height = 6, units = "in", dpi = 300)
}



##### Plot for U2AF1 WT vs U2AF1 S34F
U2AF1_S34F_all_df <- fread("/Volumes/broad_dawnccle/processed_data/missplicing_processed_df/novaseq241106/pairadise_indiv_PSI/splicelib_U2AF1_WT_splicelib_U2AF1_S34F/splicelib_U2AF1_WT_splicelib_U2AF1_S34F_rMATS_Result_P.FDR.txt")

U2AF1_S34F_all_df <- U2AF1_S34F_all_df %>% 
  mutate(
    PSI1 = mapply(calculate_ratio, I1, S1),
    PSI2 = mapply(calculate_ratio, I2, S2)
  ) %>% 
  mutate(
    PSI1_average = mapply(calculate_average, PSI1),
    PSI2_average = mapply(calculate_average, PSI2)
  ) %>%
  mutate(PSI_diff = PSI1_average - PSI2_average) %>% 
  mutate(
    count_sum_average1 = mapply(calculate_average_count_sum, I1, S1),
    count_sum_average2 = mapply(calculate_average_count_sum, I2, S2)
  ) %>% mutate(PSI_ratio = PSI1_average / PSI2_average) %>% 
  mutate(PSI_reverse_ratio = (1-PSI1_average)/(1-PSI2_average)) %>% 
  mutate(max_PSI1 = mapply(get_max_PSI, PSI1), max_PSI2 = mapply(get_max_PSI, PSI2)) %>% 
  mutate(min_PSI1 = mapply(get_min_PSI, PSI1), min_PSI2 = mapply(get_min_PSI, PSI2)) %>% 
  mutate(gene_name = map_chr(ExonID, ~str_split(.x, ";")[[1]][2])) %>% 
  filter(count_sum_average1 > 30 & count_sum_average2 > 30) %>%
  # Set the one where p-value = 0 to be 1e-10.
  mutate(FDR = ifelse(FDR == 0, 1e-10, FDR)) %>%
  # mutate(top_genes = ifelse((log2(PSI_ratio) > 2 | log2(PSI_ratio) < -2) & (-log10(FDR) > -log10(0.01)), gene_name, NA)) %>%
  mutate(top_genes = ifelse(abs(PSI_diff) > 0.2 & (-log10(FDR) > -log10(0.01)), gene_name, NA)) %>%
  mutate(is_top_gene = ifelse(!is.na(top_genes), "top_gene", "other"))

# Plot with ggplot2 and ggrepel
p1 <- ggplot(U2AF1_S34F_all_df, aes(x = PSI_diff, y = -log10(FDR), color = is_top_gene)) + 
  geom_point() + 
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed") + 
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") + 
  geom_text_repel(aes(label = top_genes), na.rm = TRUE) + 
  scale_color_manual(values = c("top_gene" = "red", "other" = "grey")) +
  theme_bw() + 
  # xlab("log2(PSI ratio)") +
  ggtitle("U2AF1 WT vs U2AF1 S34F")

ggsave("~/Dropbox (Harvard University)/02Splicing/SplicingManuscript/figure_outputs/U2AF1_S34F_volcano_plot.png", p1, width = 10, height = 6, units = "in", dpi = 300)

# Also do for ZRSR2
ZRSR2_all_df <- fread("/Volumes/broad_dawnccle/processed_data/missplicing_processed_df/novaseq241106/pairadise_indiv_PSI/splicelib_sgCh3_splicelib_ZRSR2/splicelib_sgCh3_splicelib_ZRSR2_rMATS_Result_P.FDR.txt")

ZRSR2_all_df <- ZRSR2_all_df %>% 
  mutate(
    PSI1 = mapply(calculate_ratio, I1, S1),
    PSI2 = mapply(calculate_ratio, I2, S2)
  ) %>% 
  mutate(
    PSI1_average = mapply(calculate_average, PSI1),
    PSI2_average = mapply(calculate_average, PSI2)
  ) %>%
  mutate(PSI_diff = PSI1_average - PSI2_average) %>% 
  mutate(
    count_sum_average1 = mapply(calculate_average_count_sum, I1, S1),
    count_sum_average2 = mapply(calculate_average_count_sum, I2, S2)
  ) %>% mutate(PSI_ratio = PSI1_average / PSI2_average) %>% 
  mutate(PSI_reverse_ratio = (1-PSI1_average)/(1-PSI2_average)) %>% 
  mutate(max_PSI1 = mapply(get_max_PSI, PSI1), max_PSI2 = mapply(get_max_PSI, PSI2)) %>% 
  mutate(min_PSI1 = mapply(get_min_PSI, PSI1), min_PSI2 = mapply(get_min_PSI, PSI2)) %>% 
  mutate(gene_name = map_chr(ExonID, ~str_split(.x, ";")[[1]][2])) %>% 
  filter(count_sum_average1 > 30 & count_sum_average2 > 30) %>%
  # Set the one where p-value = 0 to be 1e-10.
  mutate(FDR = ifelse(FDR == 0, 1e-10, FDR)) %>%
  # mutate(top_genes = ifelse((log2(PSI_ratio) > 2 | log2(PSI_ratio) < -2) & (-log10(FDR) > -log10(0.01)), gene_name, NA)) %>%
  mutate(top_genes = ifelse(abs(PSI_diff) > 0.2 & (-log10(FDR) > -log10(0.01)), gene_name, NA)) %>%
  mutate(is_top_gene = ifelse(!is.na(top_genes), "top_gene", "other"))

# Plot with ggplot2 and ggrepel
p2 <- ggplot(ZRSR2_all_df, aes(x = PSI_diff, y = -log10(FDR), color = is_top_gene)) + 
  geom_point() + 
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed") + 
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") + 
  geom_text_repel(aes(label = top_genes), na.rm = TRUE) + 
  scale_color_manual(values = c("top_gene" = "red", "other" = "grey")) +
  theme_bw() + 
  ggtitle("Control sgRNA vs ZRSR2 KO")
ggsave("~/Dropbox (Harvard University)/02Splicing/SplicingManuscript/figure_outputs/ZRSR2_volcano_plot.png", p2, width = 10, height = 6, units = "in", dpi = 300)



#### Also look at the 3ss data. #####

# Define the RBPs and file base path
rbps <- paste0("rbp", 1:12)
base_path <- "/Volumes/broad_dawnccle/processed_data/missplicing_processed_df/novaseq241106/pairadise_indiv_3ss/"

# Define the save path base
save_path_base <- "~/Dropbox (Harvard University)/02Splicing/SplicingManuscript/figure_outputs"

# Loop through each RBP
for (rbp in rbps) {
  # Construct the file path dynamically
  file_path <- file.path(base_path, paste0("splicelib_hek_no_dox_", rbp, "_splicelib_hek_dox_", rbp),
                         paste0("splicelib_hek_no_dox_", rbp, "_splicelib_hek_dox_", rbp, "_rMATS_Result_P.FDR.txt"))
  
  # Read and process the data
  all_df <- fread(file_path) %>%
    mutate(
      PSI1 = mapply(calculate_ratio, I1, S1),
      PSI2 = mapply(calculate_ratio, I2, S2),
      PSI1_average = mapply(calculate_average, PSI1),
      PSI2_average = mapply(calculate_average, PSI2),
      PSI_diff = PSI1_average - PSI2_average,
      count_sum_average1 = mapply(calculate_average_count_sum, I1, S1),
      count_sum_average2 = mapply(calculate_average_count_sum, I2, S2),
      PSI_ratio = PSI1_average / PSI2_average,
      PSI_reverse_ratio = (1 - PSI1_average) / (1 - PSI2_average),
      max_PSI1 = mapply(get_max_PSI, PSI1),
      max_PSI2 = mapply(get_max_PSI, PSI2),
      min_PSI1 = mapply(get_min_PSI, PSI1),
      min_PSI2 = mapply(get_min_PSI, PSI2),
      gene_name = map_chr(ExonID, ~str_split(.x, ";")[[1]][2])
    ) %>%
    filter(count_sum_average1 > 30 & count_sum_average2 > 30) %>%
    mutate(
      FDR = ifelse(FDR == 0, 1e-10, FDR),
      top_genes = ifelse(abs(PSI_diff) > 0.2 & (-log10(FDR) > -log10(0.01)), gene_name, NA),
      is_top_gene = ifelse(!is.na(top_genes), "top_gene", "other")
    )
  
  # Generate the plot
  plot <- ggplot(all_df, aes(x = log2(PSI_ratio), y = -log10(FDR), color = is_top_gene)) +
    geom_point() +
    geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
    geom_text_repel(aes(label = top_genes), na.rm = TRUE) +
    scale_color_manual(values = c("top_gene" = "red", "other" = "grey")) +
    theme_bw() +
    xlab("log2(PSI ratio)") +
    ggtitle(paste("No Dox vs Dox for", toupper(rbp)))
  
  # Construct the save path dynamically
  save_path <- file.path(save_path_base, paste0(toupper(rbp), "_volcano_plot_3ss_ratio.png"))
  
  # Save the plot
  ggsave(save_path, plot, width = 10, height = 6, units = "in", dpi = 300)
}

