library(tidyverse)
library(data.table)
library(pheatmap)

# This file is for PSI values grouped by tissue type
combined_psi <- read.table("/Volumes/broad_dawnccle/processed_data/missplicing_processed_df/V5/tissue_type_PSI/PSI_combined_output.tsv", header = T)
# This file is for combined between tissue types of "Misssplicing" (or alterative 3' splice site usage)
combined_psi <- read.table("/Volumes/broad_dawnccle/processed_data/missplicing_processed_df/V5/tissue_type/combined_output.tsv", header = T)
# This file is for individual pairwise comparisons between individual cell lines. 
combined_psi <- read.table("/Volumes/broad_dawnccle/processed_data/missplicing_processed_df/V5/pairadise_out/missplicing_indiv_combined_output.tsv", header = T)                   
# This file is for individual pairwise comparisons between individual cell lines.
combined_psi <- read.table("/Volumes/broad_dawnccle/processed_data/missplicing_processed_df/V5/pairadise_out_PSI/PSI_indiv_combined_output.tsv", header = T)
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
    

ggplot(combined_psi %>% filter(!grepl("blood", Folder)), aes(log2(PSI_ratio), count_sum_average1)) + geom_point() 

not_blood <- combined_psi %>% # filter(!grepl("K562", Folder)) %>% 
  filter(count_sum_average1 >30) %>% 
  filter(count_sum_average2 > 30) %>% 
  mutate(log2_PSI_ratio = log2(PSI_ratio), log2_PSI_reverse_ratio = log2(PSI_reverse_ratio))


# Filter out the top 10 exons with the most number of cell lines
top_seq <- not_blood %>% filter(abs(log2_PSI_ratio) > 5 | abs(log2_PSI_reverse_ratio) > 5)
top_seq_index <- top_seq %>% pull(ExonID)
test_sum <- top_seq %>% group_by(ExonID) %>% summarise(n = n()) %>% arrange(desc(n)) 
test <- top_seq %>% group_by(ExonID) %>% summarise(n = n()) %>% arrange(desc(n)) %>% filter(n>20) %>% pull(ExonID)# %>% head(100) %>% pull(ExonID)

# top_seq %>% filter(ExonID %in% test) %>% ggplot(aes(log2_PSI_ratio, log2_PSI_reverse_ratio)) + geom_point() + facet_wrap(~ExonID)

# Look at across all cell lines.
# all_sample_reps <- fread("/Volumes/broad_dawnccle/processed_data/missplicing_processed_df/V5/all_sample_reps_PSI.csv")
# all_samples_wide <- all_sample_reps %>%
#   mutate(PSI = count/(count + skipped)) %>% 
#   filter((count + skipped) > 30) %>%
#   select(-count, -skipped, - mode, -offset, -condition) %>%
#   pivot_wider(names_from = sample, values_from = PSI, values_fill = -1) 

all_samples_mat <- as.matrix(all_samples_wide[, -1])
rownames(all_samples_mat) <- all_samples_wide$index
# Pull out the top 10 exons with the most number of cell lines
top_seq_mat <- all_samples_mat[unique(test),]

# Load the CSV file
prism_cluster <- read_csv("/Volumes/broad_dawnccle/melange/data/PRISM_celltype_cluster.csv")

# Remove duplicates from prism_cluster
prism_cluster_unique <- prism_cluster %>%
  distinct(display_name, .keep_all = TRUE)

get_cluster_num_from_rep <- function(rep_name, prism_cluster_unique) {
  tmp <- str_extract(rep_name, "^[^_-]+") %>% toupper()
  cluster_num <- prism_cluster_unique %>%
    filter(display_name == tmp) %>%
    pull(cluster)
  
  # Handle cases where there is no match
  if (length(cluster_num) == 0) {
    return(-1)
  } else {
    return(cluster_num)
  }
}
get_lineage_from_rep <- function(rep_name, prism_cluster_unique) {
  tmp <- str_extract(rep_name, "^[^_-]+") %>% toupper()
  lineage <- prism_cluster_unique %>%
    filter(display_name == tmp) %>%
    pull(lineage)
  
  # Handle cases where there is no match
  if (length(lineage) == 0) {
    return("Unknown")
  } else {
    return(lineage)
  }
}

get_subtype_from_rep <- function(rep_name, prism_cluster_unique) {
  tmp <- str_extract(rep_name, "^[^_-]+") %>% toupper()
  subtype <- prism_cluster_unique %>%
    filter(display_name == tmp) %>%
    pull(subtype)
  
  # Handle cases where there is no match
  if (length(subtype) == 0) {
    return("Unknown")
  } else {
    return(subtype)
  }
}

annotation_col <- data.frame(
  cell_line = colnames(top_seq_mat)
)
annotation_col$cluster <- sapply(annotation_col$cell_line, get_cluster_num_from_rep, prism_cluster_unique)
annotation_col$lineage <- sapply(annotation_col$cell_line, get_lineage_from_rep, prism_cluster_unique)
annotation_col$subtype <- sapply(annotation_col$cell_line, get_subtype_from_rep, prism_cluster_unique)

# Convert to a data frame and set rownames
annotation_col <- as.data.frame(annotation_col)
rownames(annotation_col) <- annotation_col$cell_line

# Remove cell_line column and convert cluster to factor
annotation_col <- annotation_col %>% 
  select(-cell_line) %>% 
  mutate(cluster = as.factor(cluster)) 

# Order by annotation_col
annotation_col <- annotation_col %>% arrange(lineage, subtype)
ordered_columns <- rownames(annotation_col)
top_seq_mat <- top_seq_mat[, ordered_columns]

# Find the positions where the lineage change to add gaps. Lineage is not numeric.
lineage_changes <- which(annotation_col$lineage != lag(annotation_col$lineage, default = first(annotation_col$lineage)))
lineage_changes <- lineage_changes - 1  # Shift to get the correct column positions

# Plot heatmap with annotations and gaps between clusters
pheatmap(top_seq_mat,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = T,
         show_colnames = TRUE,
         fontsize = 6,
         color = colorRampPalette(c("black", "black","black","black", "#10194d","#5496b6", "#ffffe0", "#d96c5d","#4a001e"))(100),
         annotation_col = annotation_col,
         gaps_col = lineage_changes)

###### Look at top seq where folder has Kelly. Kelly has a specific block.
# Get sequences from top_seq_mat where Kelly > 0.7 and COGN278 < 0.3. The names are like Kelly-rep1, COGN278-rep1 etc
# Extract column indices for Kelly and COGN278
kelly_cols <- grep("^Kelly", colnames(top_seq_mat))
cogn278_cols <- grep("^COGN278", colnames(top_seq_mat))

# Create logical vectors for the conditions
kelly_condition <- rowSums(top_seq_mat[, kelly_cols] > 0.7) == 3
cogn278_condition <- rowSums(top_seq_mat[, cogn278_cols] < 0.3) == 3

# Get the rows that meet both conditions
filtered_rows <- which(kelly_condition & cogn278_condition)

# Extract the sequences that meet the conditions
filtered_seq_mat <- top_seq_mat[filtered_rows, ]

# Plot the heatmap
pheatmap(filtered_seq_mat,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = TRUE,
         fontsize = 6,
         color = colorRampPalette(c("black", "black", "#45abd7", "white", "#cf1c1a"))(100),
         annotation_col = annotation_col,
         gaps_col = cluster_changes)


# SLC39A8 subset
SLC39A8 <- read_csv("/Volumes/broad_dawnccle/processed_data/missplicing_processed_df/V5/SLC39A8_subset.csv", col_names = F)
colnames(SLC39A8) <- c("index", "mode", "offset","count" , "sample", "condition")
SLC39A8 <- SLC39A8 %>% group_by(sample, condition) %>% 
  filter(mode != "UNSPLICED") %>% 
  mutate(count_sum = sum(count)) %>%
  ungroup() %>% 
  mutate(fraction = count / count_sum) 
  # group_by(condition, mode, offset) %>%
  # summarise(fraction = mean(fraction)) 
  
  
ggplot(SLC39A8 %>% filter(condition %in% c("HEK", "Kelly", "IPC298", "MEWO")), aes(offset, fraction)) + 
    geom_bar(stat = "identity") + 
    facet_wrap(~sample, scales = "free_x") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    # no legend
    theme(legend.position = "none") 

#### Look at kidney and liver
kidney_celllines <- c("786O", "769P", "ACHN", "KMRC1", "KMRC20", "ACHN", "VMRCRCZ")
liver_celllines <- c("SNU449", "SNU423", "SNU398", "PLCPRF5", "JHH6")
breast_celllines <- c("T47D", "MCF7", "HCC1428", "MDAMB231", "CAL120", "HCC38")
top_seq <- not_blood %>% filter(abs(log2_PSI_ratio) > 4.5 | abs(log2_PSI_reverse_ratio) > 4.5)
top_seq_kidney_liver <- top_seq %>% 
  filter(grepl(paste(kidney_celllines, collapse = "|"), Folder)) %>%
  filter(grepl(paste(liver_celllines, collapse = "|"), Folder)) %>%
  filter(count_sum_average1 > 30) %>%
  filter(count_sum_average2 > 30) %>% 
  mutate(log2_PSI_ratio = log2(PSI_ratio), log2_PSI_reverse_ratio = log2(PSI_reverse_ratio))

# Merge with the barcodes 
top_seq_kidney_liver_bc <- merge(top_seq_kidney_liver, barcodes, by.x = "ExonID", by.y = "ID")
write_csv(top_seq_kidney_liver_bc, "/Volumes/broad_dawnccle/processed_data/missplicing_processed_df/V5/top_seq_kidney_liver.csv")

# Heatmap of the top seq kidney and liver
top_seq_mat_kidney_liver <- all_samples_mat[unique(top_seq_kidney_liver$ExonID),]

# Order by annotation_col
top_seq_mat_kidney_liver <- top_seq_mat_kidney_liver[, ordered_columns]

# Plot heatmap with annotations and gaps between clusters
pheatmap(top_seq_mat_kidney_liver,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = T,
         show_colnames = TRUE,
         fontsize = 6,
         color = colorRampPalette(c("black", "black", "#45abd7", "white", "#cf1c1a"))(100),
         # color = colorRampPalette(c("black", "black","black","black", "#10194d","#5496b6", "#ffffe0", "#d96c5d","#4a001e"))(100),
         annotation_col = annotation_col,
         gaps_col = lineage_changes,
         main = "Top seq in kidney and liver")

# Also look at liver vs breast
top_seq_liver_breast <- top_seq %>% 
  filter(grepl(paste(liver_celllines, collapse = "|"), Folder)) %>%
  filter(grepl(paste(breast_celllines, collapse = "|"), Folder)) %>%
  filter(count_sum_average1 > 30) %>%
  filter(count_sum_average2 > 30) %>% 
  mutate(log2_PSI_ratio = log2(PSI_ratio), log2_PSI_reverse_ratio = log2(PSI_reverse_ratio))

# Merge with the barcodes
top_seq_liver_breast_bc <- merge(top_seq_liver_breast, barcodes, by.x = "ExonID", by.y = "ID")
write_csv(top_seq_liver_breast_bc, "/Volumes/broad_dawnccle/processed_data/missplicing_processed_df/V5/top_seq_liver_breast.csv")

# Heatmap of the top seq liver and breast
top_seq_mat_liver_breast <- all_samples_mat[unique(top_seq_liver_breast$ExonID),]

# Order by annotation_col
top_seq_mat_liver_breast <- top_seq_mat_liver_breast[, ordered_columns]

# Plot heatmap with annotations and gaps between clusters
pheatmap(top_seq_mat_liver_breast,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = T,
         show_colnames = TRUE,
         fontsize = 6,
         color = colorRampPalette(c("black", "black", "#45abd7", "white", "#cf1c1a"))(100),
         # color = colorRampPalette(c("black", "black","black","black", "#10194d","#5496b6", "#ffffe0", "#d96c5d","#4a001e"))(100),
         annotation_col = annotation_col,
         gaps_col = lineage_changes,
         main = "Top seq in liver and breast")
