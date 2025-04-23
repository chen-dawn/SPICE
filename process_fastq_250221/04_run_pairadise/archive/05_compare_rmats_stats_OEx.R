library(tidyverse)
library(pheatmap)
library(data.table)

######## Look at rmats stats ########
combined_psi <- read_tsv("U:/processed_data/reprocess_250221/pairadise_indiv_PSI/OEx_PSI_combined_output_indiv.tsv")
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

# Apply the function to the data frame
combined_psi <- combined_psi %>%
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
  mutate(PSI_reverse_ratio = (1-PSI1_average)/(1-PSI2_average))

combined_psi_filtered <- combined_psi %>% 
  filter(count_sum_average1 >30) %>% 
  filter(count_sum_average2 > 30) %>% 
  mutate(log2_PSI_ratio = log2(PSI_ratio), log2_PSI_reverse_ratio = log2(PSI_reverse_ratio)) %>% 
  separate(ExonID, sep = "__", into =c("index", "offset"), remove = FALSE) %>%
  separate(offset, into = c("skipped_exon_start", "skipped_exon_end", "downstream_exon_start"), sep = ":", remove = FALSE) %>%
  filter(abs(as.integer(skipped_exon_start)) != 1 & abs(as.integer(skipped_exon_end)) != 1) 

# Plot the volcano plot
ggplot(combined_psi_filtered, aes(PSI_diff, -log10(FDR))) + 
  geom_point() + 
  facet_wrap(~Folder)

# Look at the non rbp 7 and 8
non_rbp_7_8 <- combined_psi_filtered %>% 
  filter(!grepl("rbp7|rbp8", Folder)) 
  
