library(tidyverse)
library(vroom)
library(data.table)
library(Biostrings)
library(ggpointdensity)
library(pheatmap)
combined_psi <- fread("~/Dropbox (Harvard University)/02Splicing/latest/Nova240826_PSI_combined_output_indiv.tsv")
# combined_psi <- fread("~/Dropbox (Harvard University)/02Splicing/latest/Nova240826_3ss_combined_output_indiv.tsv")
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

get_max_PSI <- function(PSI){
  PSI_values <- as.numeric(unlist(strsplit(PSI, ",")))
  max_PSI <- max(PSI_values)
  return(round(max_PSI, 3))
}

get_min_PSI <- function(PSI){
  PSI_values <- as.numeric(unlist(strsplit(PSI, ",")))
  min_PSI <- min(PSI_values)
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
RBM5_all_df <- fread("/Volumes/broad_dawnccle/processed_data/missplicing_processed_df/novaseq240826/pairadise_indiv_PSI/MUT-sgCh3-1_MUT-sgRBM5/MUT-sgCh3-1_MUT-sgRBM5_rMATS_Result_P.FDR.txt")

RBM5_all_df <- RBM5_all_df %>% 
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

library(ggrepel)
# Add a column to identify the top genes based on a condition
RBM5_all_df <- RBM5_all_df %>%
  # Set the one where p-value = 0 to be 1e-10.
  mutate(FDR = ifelse(FDR == 0, 1e-10, FDR)) %>%
  mutate(top_genes = ifelse((log2(PSI_ratio) > 2 | log2(PSI_ratio) < -2) & (-log10(FDR) > -log10(0.01)), gene_name, NA)) %>%
  mutate(is_top_gene = ifelse(!is.na(top_genes), "top_gene", "other"))

# Plot with ggplot2 and ggrepel
ggplot(RBM5_all_df, aes(x = log2(PSI_reverse_ratio), y = -log10(FDR), color = is_top_gene)) + 
  geom_point() + 
  # geom_vline(xintercept = c(-2, 2), linetype = "dashed") + 
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") + 
  geom_text_repel(aes(label = top_genes), na.rm = TRUE) + 
  scale_color_manual(values = c("top_gene" = "red", "other" = "grey")) +
  theme_minimal()
