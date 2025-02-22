library(tidyverse)
library(data.table)
library(pheatmap)

# # This file is for PSI values grouped by tissue type
# combined_psi <- read.table("/Volumes/broad_dawnccle/processed_data/missplicing_processed_df/V5/tissue_type_PSI/PSI_combined_output.tsv", header = T)
# # This file is for combined between tissue types of "Misssplicing" (or alterative 3' splice site usage)
# combined_psi <- read.table("/Volumes/broad_dawnccle/processed_data/missplicing_processed_df/V5/tissue_type/combined_output.tsv", header = T)
# # This file is for individual pairwise comparisons between individual cell lines. 
# combined_psi <- read.table("/Volumes/broad_dawnccle/processed_data/missplicing_processed_df/V5/pairadise_out/missplicing_indiv_combined_output.tsv", header = T)                   
# # This file is for individual pairwise comparisons between individual cell lines.
# combined_psi <- read.table("/Volumes/broad_dawnccle/processed_data/missplicing_processed_df/V5/pairadise_out_PSI/PSI_indiv_combined_output.tsv", header = T)
# # This file is for individual one vs all comparisons between individual cell lines.
# combined_psi <- read.table("/Volumes/broad_dawnccle/processed_data/missplicing_processed_df/V5/pairadise_out_PSI_indiv/PSI_combined_output_indiv.tsv", header = T)
# Updated PSI.
combined_psi <- read.table("~/Dropbox (Harvard University)/02Splicing/latest/rmats_one_vs_all_combined_output_PSI.tsv", header = T)
all_sample_reps <- fread("~/Dropbox (Harvard University)/02Splicing/latest/all_sample_reps_PSI.csv")
# Load the CSV file
prism_cluster <- read_csv("~/melange/data/PRISM_celltype_cluster.csv")

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
  

combined_psi_filtered <- combined_psi %>% # filter(!grepl("K562", Folder)) %>% 
  filter(count_sum_average1 >30) %>% 
  filter(count_sum_average2 > 30) %>% 
  mutate(log2_PSI_ratio = log2(PSI_ratio), log2_PSI_reverse_ratio = log2(PSI_reverse_ratio))


# Filter out the top 10 exons with the most number of cell lines
top_seq <- combined_psi_filtered %>% filter(abs(log2_PSI_ratio) > 3 | abs(log2_PSI_reverse_ratio) > 2)
top_seq_index <- top_seq %>% pull(ExonID)
test_sum <- top_seq %>% group_by(ExonID) %>% summarise(n = n()) %>% arrange(desc(n)) 
test <- top_seq %>% 
  filter(!(Folder %in% c("HEK"))) %>% 
  filter(abs(log2_PSI_ratio) > 4 | abs(log2_PSI_reverse_ratio) > 4) %>% pull(ExonID) %>% unique() # group_by(ExonID) %>% summarise(n = n()) %>% arrange(desc(n)) %>% head(100) %>% pull(ExonID)

# top_seq %>% filter(ExonID %in% test) %>% ggplot(aes(log2_PSI_ratio, log2_PSI_reverse_ratio)) + geom_point() + facet_wrap(~ExonID)

######## Look at across all cell lines. ########
all_samples_wide <- all_sample_reps %>%
  mutate(PSI = count/(count + skipped)) %>%
  filter((count + skipped) > 30) %>%
  select(-count, -skipped, - mode, -offset, -condition) %>%
  pivot_wider(names_from = sample, values_from = PSI, values_fill = -1)

all_samples_mat <- as.matrix(all_samples_wide[, -1])
rownames(all_samples_mat) <- all_samples_wide$index
# Pull out the top 10 exons with the most number of cell lines
top_seq_mat <- all_samples_mat[unique(test),]


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
# pheatmap(top_seq_mat,
#          cluster_rows = TRUE,
#          cluster_cols = FALSE,
#          show_rownames = T,
#          show_colnames = TRUE,
#          fontsize = 6,
#          color = colorRampPalette(c("black", "black","black","black", "#10194d","#5496b6", "#ffffe0", "#d96c5d","#4a001e"))(100),
#          annotation_col = annotation_col,
#          gaps_col = lineage_changes)


plot_heatmap_for_ids <- function(exon_ids, all_sample_reps, prism_cluster, title = "") {
  # Prepare the wide format of all samples
  all_samples_wide <- all_sample_reps %>%
    mutate(PSI = count / (count + skipped)) %>%
    filter((count + skipped) > 30) %>%
    select(-count, -skipped, -mode, -offset, -condition) %>%
    pivot_wider(names_from = sample, values_from = PSI, values_fill = -1)
  
  all_samples_mat <- as.matrix(all_samples_wide[, -1])
  rownames(all_samples_mat) <- all_samples_wide$index
  
  # Remove duplicates from prism_cluster
  prism_cluster_unique <- prism_cluster %>%
    distinct(display_name, .keep_all = TRUE)
  
  # Loop through the list of ExonIDs and plot heatmaps
    top_seq_mat <- all_samples_mat[exon_ids,]
    
    annotation_col <- data.frame(
      cell_line = colnames(top_seq_mat)
    )
    annotation_col$cluster <- sapply(annotation_col$cell_line, get_cluster_num_from_rep, prism_cluster_unique)
    annotation_col$lineage <- sapply(annotation_col$cell_line, get_lineage_from_rep, prism_cluster_unique)
    annotation_col$subtype <- sapply(annotation_col$cell_line, get_subtype_from_rep, prism_cluster_unique)
    
    annotation_col <- as.data.frame(annotation_col)
    rownames(annotation_col) <- annotation_col$cell_line
    
    annotation_col <- annotation_col %>%
      select(-cell_line) %>%
      mutate(cluster = as.factor(cluster)) %>%
      arrange(lineage, subtype)
    
    ordered_columns <- rownames(annotation_col)
    top_seq_mat <- top_seq_mat[, ordered_columns]
    
    lineage_changes <- which(annotation_col$lineage != lag(annotation_col$lineage, default = first(annotation_col$lineage)))
    lineage_changes <- lineage_changes - 1
    
    pheatmap(top_seq_mat,
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             show_rownames = TRUE,
             show_colnames = TRUE,
             fontsize = 6,
             color = colorRampPalette(c("black", "black", "black", "black", "#10194d", "#5496b6", "#ffffe0", "#d96c5d", "#4a001e"))(100),
             annotation_col = annotation_col,
             gaps_col = lineage_changes,
             main = title)
    
  return(top_seq_mat)
}

# plot_heatmap_for_ids(test[1:2], all_sample_reps, prism_cluster)

############################################
##### Look at top sequences for Kelly ######
############################################
Kelly <- combined_psi_filtered %>% filter(grepl("Kelly", Folder))
Kelly_top <- Kelly %>% filter(abs(log2_PSI_ratio) > 3 | abs(log2_PSI_reverse_ratio) > 3 | abs(PSI_diff) > 0.6)

Kelly_top_mat <- plot_heatmap_for_ids(Kelly_top %>% pull(ExonID), all_sample_reps, prism_cluster)
# Get num of -1 for each row. 
num_missing <- rowSums(Kelly_top_mat == -1)
# Get rows with at most 75 missing.
Kelly_top_mat_filtered <- Kelly_top_mat[num_missing <= 75,]
# Plot the heatmap
final_kelly_top <- plot_heatmap_for_ids(rownames(Kelly_top_mat_filtered), all_sample_reps, prism_cluster, title = "Top seq in Kelly")
# Write the final top seq to a CSV file.
write.table(as.data.frame(final_kelly_top), "/Volumes/broad_dawnccle/processed_data/latest/mutagenesis_shortlist/kelly_top_mat.tsv", row.names = T, sep = "\t")
write_csv(Kelly_top %>% filter(ExonID %in% rownames(final_kelly_top)), "/Volumes/broad_dawnccle/processed_data/latest/mutagenesis_shortlist/kelly_top_stats.csv")

############################################
##### Look at top sequences for JHH6 ######
############################################
JHH6 <- combined_psi_filtered %>% filter(grepl("JHH6", Folder))
JHH6_top <- JHH6 %>% filter(abs(log2_PSI_ratio) > 2 | abs(log2_PSI_reverse_ratio) > 2 | abs(PSI_diff) > 0.3)

JHH6_top_mat <- plot_heatmap_for_ids(JHH6_top %>% pull(ExonID), all_sample_reps, prism_cluster)
# Get num of -1 for each row.
num_missing <- rowSums(JHH6_top_mat == -1)
# Get rows with at most 75 missing.
JHH6_top_mat_filtered <- JHH6_top_mat[num_missing <= 75,]
# Plot the heatmap
final_JHH6_top <- plot_heatmap_for_ids(rownames(JHH6_top_mat_filtered), all_sample_reps, prism_cluster, title = "Top seq in JHH6")
# Write the final top seq to a CSV file.
write.table(as.data.frame(final_JHH6_top), "/Volumes/broad_dawnccle/processed_data/latest/mutagenesis_shortlist/JHH6_top_mat.tsv", row.names = T, sep = "\t")
write_csv(JHH6_top %>% filter(ExonID %in% rownames(final_JHH6_top)), "/Volumes/broad_dawnccle/processed_data/latest/mutagenesis_shortlist/JHH6_top_stats.csv")

############################################
##### Look at top sequences for MCF7 ######
############################################
MCF7 <- combined_psi_filtered %>% filter(grepl("MCF7", Folder))
MCF7_top <- MCF7 %>% filter(abs(log2_PSI_ratio) > 2 | abs(log2_PSI_reverse_ratio) > 2 | abs(PSI_diff) > 0.3)

MCF7_top_mat <- plot_heatmap_for_ids(MCF7_top %>% pull(ExonID), all_sample_reps, prism_cluster)
# Get num of -1 for each row.
num_missing <- rowSums(MCF7_top_mat == -1)
# Get rows with at most 75 missing.
MCF7_top_mat_filtered <- MCF7_top_mat[num_missing <= 75,]
# Plot the heatmap
final_MCF7_top <- plot_heatmap_for_ids(rownames(MCF7_top_mat_filtered), all_sample_reps, prism_cluster, title = "Top seq in MCF7")
# Write the final top seq to a CSV file.
write.table(as.data.frame(final_MCF7_top), "/Volumes/broad_dawnccle/processed_data/latest/mutagenesis_shortlist/MCF7_top_mat.tsv", row.names = T, sep = "\t")
write_csv(MCF7_top %>% filter(ExonID %in% rownames(final_MCF7_top)), "/Volumes/broad_dawnccle/processed_data/latest/mutagenesis_shortlist/MCF7_top_stats.csv")

############################################
##### Look at top sequences for IPC298 ######
############################################
IPC298 <- combined_psi_filtered %>% filter(grepl("IPC298", Folder))
IPC298_top <- IPC298 %>% filter(abs(log2_PSI_ratio) > 2 | abs(log2_PSI_reverse_ratio) > 2 | abs(PSI_diff) > 0.3)

IPC298_top_mat <- plot_heatmap_for_ids(IPC298_top %>% pull(ExonID), all_sample_reps, prism_cluster)
# Get num of -1 for each row.
num_missing <- rowSums(IPC298_top_mat == -1)
# Get rows with at most 75 missing.
IPC298_top_mat_filtered <- IPC298_top_mat[num_missing <= 75,]
# Plot the heatmap
final_IPC298_top <- plot_heatmap_for_ids(rownames(IPC298_top_mat_filtered), all_sample_reps, prism_cluster, title = "Top seq in IPC298")
# Write the final top seq to a CSV file.
write.table(as.data.frame(final_IPC298_top), "/Volumes/broad_dawnccle/processed_data/latest/mutagenesis_shortlist/IPC298_top_mat.tsv", row.names = T, sep = "\t")
write_csv(IPC298_top %>% filter(ExonID %in% rownames(final_IPC298_top)), "/Volumes/broad_dawnccle/processed_data/latest/mutagenesis_shortlist/IPC298_top_stats.csv")

############################################
##### Look at top sequences for IGR37 ######
############################################
IGR37 <- combined_psi_filtered %>% filter(grepl("IGR37", Folder))
IGR37_top <- IGR37 %>% filter( abs(log2_PSI_ratio) > 2 | abs(log2_PSI_reverse_ratio) > 2 | abs(PSI_diff) > 0.3)

IGR37_top_mat <- plot_heatmap_for_ids(IGR37_top %>% pull(ExonID), all_sample_reps, prism_cluster)
# Get num of -1 for each row.
num_missing <- rowSums(IGR37_top_mat == -1)
# Get rows with at most 75 missing.
IGR37_top_mat_filtered <- IGR37_top_mat[num_missing <= 75,]
# Plot the heatmap
final_IGR37_top <- plot_heatmap_for_ids(rownames(IGR37_top_mat_filtered), all_sample_reps, prism_cluster, title = "Top seq in IGR37")
# Write the final top seq to a CSV file.
write.table(as.data.frame(final_IGR37_top), "/Volumes/broad_dawnccle/processed_data/latest/mutagenesis_shortlist/IGR37_top_mat.tsv", row.names = T, sep = "\t")
write_csv(IGR37_top %>% filter(ExonID %in% rownames(final_IGR37_top)), "/Volumes/broad_dawnccle/processed_data/latest/mutagenesis_shortlist/IGR37_top_stats.csv")

############################################
##### Look at top sequences for KMRC20 ######
############################################
KMRC20 <- combined_psi_filtered %>% filter(grepl("KMRC20", Folder))
KMRC20_top <- KMRC20 %>% filter(abs(log2_PSI_ratio) > 3 | abs(log2_PSI_reverse_ratio) > 3 | abs(PSI_diff) > 0.5)

KMRC20_top_mat <- plot_heatmap_for_ids(KMRC20_top %>% pull(ExonID), all_sample_reps, prism_cluster)
# Get num of -1 for each row.
num_missing <- rowSums(KMRC20_top_mat == -1)
# Get rows with at most 75 missing.
KMRC20_top_mat_filtered <- KMRC20_top_mat[num_missing <= 75,]
# Plot the heatmap
final_KMRC20_top <- plot_heatmap_for_ids(rownames(KMRC20_top_mat_filtered), all_sample_reps, prism_cluster, title = "Top seq in KMRC20")
# Write the final top seq to a CSV file.
write.table(as.data.frame(final_KMRC20_top), "/Volumes/broad_dawnccle/processed_data/latest/mutagenesis_shortlist/KMRC20_top_mat.tsv", row.names = T, sep = "\t")
write_csv(KMRC20_top %>% filter(ExonID %in% rownames(final_KMRC20_top)), "/Volumes/broad_dawnccle/processed_data/latest/mutagenesis_shortlist/KMRC20_top_stats.csv")

############################################
##### Look at top sequences for HEK ######
############################################
HEK <- combined_psi_filtered %>% filter(grepl("HEK", Folder))
HEK_top <- HEK %>% filter( abs(PSI_diff) > 0.6)

HEK_top_mat <- plot_heatmap_for_ids(HEK_top %>% pull(ExonID), all_sample_reps, prism_cluster)
# Get num of -1 for each row.
num_missing <- rowSums(HEK_top_mat == -1)
# Get rows with at most 75 missing.
HEK_top_mat_filtered <- HEK_top_mat[num_missing <= 75,]
# Plot the heatmap
final_HEK_top <- plot_heatmap_for_ids(rownames(HEK_top_mat_filtered), all_sample_reps, prism_cluster, title = "Top seq in HEK")
# Write the final top seq to a CSV file.
write.table(as.data.frame(final_HEK_top), "/Volumes/broad_dawnccle/processed_data/latest/mutagenesis_shortlist/HEK_top_mat.tsv", row.names = T, sep = "\t")
write_csv(HEK_top %>% filter(ExonID %in% rownames(final_HEK_top)), "/Volumes/broad_dawnccle/processed_data/latest/mutagenesis_shortlist/HEK_top_stats.csv")


###### Retrying the filter #####
constant_in_most <- combined_psi_filtered %>% 
  filter(PSI2_average > 0.93 | PSI2_average < 0.05) %>% 
  filter(!(Folder %in% c("K562", "HEK")))

shorter_shortlist <- constant_in_most %>% 
  filter(Folder %in% c("Kelly", "T47D", "HCC38", "KMRC20"))

best_hits_top_mat <- plot_heatmap_for_ids(constant_in_most %>% pull(ExonID)%>%unique(), all_sample_reps, prism_cluster)
# Get num of -1 for each row.
num_missing <- rowSums(best_hits_top_mat == -1)
# Get rows with at most 75 missing.
best_hits_top_mat_filtered <- best_hits_top_mat[num_missing <= 75,]

# Plot the heatmap
final_best_hits <- plot_heatmap_for_ids(unique(rownames(best_hits_top_mat_filtered)), all_sample_reps, prism_cluster, title = "Constant in most cell lines")
# Write the final top seq to a CSV file.
write.table(as.data.frame(final_best_hits), "/Volumes/broad_dawnccle/processed_data/latest/mutagenesis_shortlist/best_hits_top_mat.tsv", row.names = T, sep = "\t")
write_csv(constant_in_most %>% filter(ExonID %in% rownames(final_best_hits)), "/Volumes/broad_dawnccle/processed_data/latest/mutagenesis_shortlist/best_hits_top_stats.csv")


####### Get the sequences for target cell lines ######
high_ratio_and_fold_change <- combined_psi_filtered %>% 
  mutate(max_log2_abs_ratio = pmax(abs(log2_PSI_ratio), abs(log2_PSI_reverse_ratio))) %>%
  arrange(desc(max_log2_abs_ratio)) %>% 
  filter(abs(PSI_diff) > 0.3)

get_top_ID_for_cellline <- function(df, cellline, top_n = 4){
  topseq <- df %>% filter(Folder == cellline) %>% head(top_n) 
  return(topseq)
}

JHH6_top_seq <- get_top_ID_for_cellline(high_ratio_and_fold_change, "JHH6")
MCF7_top_seq <- get_top_ID_for_cellline(high_ratio_and_fold_change, "MCF7")
IPC298_top_seq <- get_top_ID_for_cellline(high_ratio_and_fold_change, "IPC298")
IGR37_top_seq <- get_top_ID_for_cellline(high_ratio_and_fold_change, "IGR37")
KMRC20_top_seq <- get_top_ID_for_cellline(high_ratio_and_fold_change, "KMRC20")
HEK_top_seq <- get_top_ID_for_cellline(high_ratio_and_fold_change, "HEK")
Kelly_top_seq <- get_top_ID_for_cellline(high_ratio_and_fold_change, "Kelly")
T47D_top_seq <- get_top_ID_for_cellline(high_ratio_and_fold_change, "T47D")
HCC38_top_seq <- get_top_ID_for_cellline(high_ratio_and_fold_change, "HCC38")
all_top_seq <- bind_rows(JHH6_top_seq, MCF7_top_seq, IPC298_top_seq, IGR37_top_seq, KMRC20_top_seq, HEK_top_seq, Kelly_top_seq, T47D_top_seq, HCC38_top_seq)

# Plot heatmap for top constant
final_top_seq <- plot_heatmap_for_ids(all_top_seq$ExonID, all_sample_reps, prism_cluster, title = "Top variable in target cell lines")
# Get num of -1 for each row.
num_missing <- rowSums(final_top_seq == -1)
# Get rows with at most 75 missing.
final_top_seq_filtered <- final_top_seq[num_missing <= 75,]
# Plot the heatmap
final_top_seq_filtered <- plot_heatmap_for_ids(rownames(final_top_seq_filtered), all_sample_reps, prism_cluster, title = "Top variable in target cell lines")
all_top_variable <- all_top_seq %>% filter(ExonID %in% rownames(final_top_seq_filtered)) 
# Write the final top seq to a CSV file.
write_csv(all_top_variable, "/Volumes/broad_dawnccle/processed_data/latest/validation_shortlist/all_top_variable.csv")
write.table(as.data.frame(final_top_seq_filtered), "/Volumes/broad_dawnccle/processed_data/latest/validation_shortlist/all_top_variable_mat.tsv", row.names = T, sep = "\t")

all_samples_wide <- all_sample_reps %>%
  mutate(PSI = count / (count + skipped)) %>%
  filter((count + skipped) > 30) %>%
  select(-count, -skipped, -mode, -offset, -condition) %>%
  pivot_wider(names_from = sample, values_from = PSI, values_fill = NA)

all_samples_mat <- as.matrix(all_samples_wide[, -1])
rownames(all_samples_mat) <- all_samples_wide$index

min_each_row <- apply(all_samples_mat, 1, min, na.rm = T)
max_each_row <- apply(all_samples_mat, 1, max, na.rm = T)
num_NA_each_row <- apply(all_samples_mat, 1, function(x) sum(is.na(x)))

min_max_df <- data.frame(ExonID = rownames(all_samples_mat), min = min_each_row, max = max_each_row, num_NA = num_NA_each_row)
# Clear row names.
rownames(min_max_df) <- NULL
const_max <- min_max_df %>% filter(max < 0.04) %>% filter(num_NA <= 45)
const_min <- min_max_df %>% filter(min > 0.875) %>% filter(num_NA <= 45)
top_constant <- rbind(const_max, const_min)
# Write the final top seq to a CSV file.
write_csv(top_constant, "/Volumes/broad_dawnccle/processed_data/latest/validation_shortlist/top_constant.csv")
# Plot heatmap for top constant
final_top_constant <- plot_heatmap_for_ids(top_constant$ExonID, all_sample_reps, prism_cluster, title = "Top constant in most cell lines")
# Write table to a file.
write.table(as.data.frame(final_top_constant), "/Volumes/broad_dawnccle/processed_data/latest/validation_shortlist/top_constant_mat.tsv", row.names = T, sep = "\t")

ggplot(all_top_variable, aes(PSI1_average, PSI2_average, color = Folder, label = Folder)) + geom_point() + 
  geom_label() + theme_bw() + 
  xlab("Average PSI in target cell line") + 
  ylab("Average PSI in non-target cell line")


# Get the raw sequences for these seuqences. 
revcomp <- function(seq) {
  # Create a named vector for complement bases
  complement <- c(A = "T", T = "A", G = "C", C = "G")
  
  # Reverse the sequence and find complement
  rev_seq <- rev(strsplit(seq, NULL)[[1]])
  rev_comp_seq <- sapply(rev_seq, function(base) complement[base])
  
  # Combine into a single string
  paste(rev_comp_seq, collapse = "")
}
twist_seq <- read_csv("~/melange/data/guide_library_cleaned/20240605_twist_library_v3_filtered.csv")

# Get the sequences in the Kelly_top_seq
top_seq_sequences <- high_ratio_and_fold_change %>% left_join(twist_seq, by = c("ExonID" = "ID")) %>% 
  filter(Folder == "Kelly")

write_csv(top_seq_sequences, "/Volumes/broad_dawnccle/processed_data/latest/validation_shortlist/KELLY_ONLY_top_nucleotide_sequence.csv")






############## Look at tissue type level data ##############
# Load the CSV file
prism_cluster <- read_csv("~/melange/data/PRISM_celltype_cluster.csv")
combined_psi <- read.table("~/Dropbox (Harvard University)/02Splicing/latest/rmats_tissue_type_combined_output_PSI.tsv", header = T)

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
  mutate(PSI_reverse_ratio = (1-PSI1_average)/(1-PSI2_average)) %>% 
  mutate(max_PSI1 = mapply(get_max_PSI, PSI1), max_PSI2 = mapply(get_max_PSI, PSI2)) %>% 
  mutate(min_PSI1 = mapply(get_min_PSI, PSI1), min_PSI2 = mapply(get_min_PSI, PSI2))


combined_psi_filtered <- combined_psi %>% 
  filter(count_sum_average1 >30) %>% 
  filter(count_sum_average2 > 30) %>% 
  mutate(log2_PSI_ratio = log2(PSI_ratio), log2_PSI_reverse_ratio = log2(PSI_reverse_ratio))


############################################
##### Look at top sequences for breast_kidney ######
############################################
breast_kidney <- combined_psi_filtered %>% filter(grepl("breast_kidney", Folder))
breast_kidney_top <- breast_kidney %>% filter(abs(log2_PSI_ratio) > 2 | abs(log2_PSI_reverse_ratio) > 2 | abs(PSI_diff) > 0.3)

breast_kidney_top_mat <- plot_heatmap_for_ids(breast_kidney_top %>% pull(ExonID), all_sample_reps, prism_cluster)
# Get num of -1 for each row.
num_missing <- rowSums(breast_kidney_top_mat == -1)
# Get rows with at most 75 missing.
breast_kidney_top_mat_filtered <- breast_kidney_top_mat[num_missing <= 75,]
# Plot the heatmap
final_breast_kidney_top <- plot_heatmap_for_ids(rownames(breast_kidney_top_mat_filtered), all_sample_reps, prism_cluster, title = "Top seq in breast_kidney")
# Write the final top seq to a CSV file.
write.table(as.data.frame(final_breast_kidney_top), "/Volumes/broad_dawnccle/processed_data/latest/mutagenesis_shortlist/breast_kidney_top_mat.tsv", row.names = T, sep = "\t")
write_csv(breast_kidney_top %>% filter(ExonID %in% rownames(final_breast_kidney_top)), "/Volumes/broad_dawnccle/processed_data/latest/mutagenesis_shortlist/breast_kidney_top_stats.csv")

############################################
##### Look at top sequences for breast_skin ######
############################################
breast_skin <- combined_psi_filtered %>% filter(grepl("breast_skin", Folder))
breast_skin_top <- breast_skin %>% filter(abs(log2_PSI_ratio) > 2 | abs(log2_PSI_reverse_ratio) > 2 | abs(PSI_diff) > 0.3)

breast_skin_top_mat <- plot_heatmap_for_ids(breast_skin_top %>% pull(ExonID), all_sample_reps, prism_cluster)
# Get num of -1 for each row.
num_missing <- rowSums(breast_skin_top_mat == -1)
# Get rows with at most 75 missing.
breast_skin_top_mat_filtered <- breast_skin_top_mat[num_missing <= 75,]
# Plot the heatmap
final_breast_skin_top <- plot_heatmap_for_ids(rownames(breast_skin_top_mat_filtered), all_sample_reps, prism_cluster, title = "Top seq in breast_skin")
# Write the final top seq to a CSV file.
write.table(as.data.frame(final_breast_skin_top), "/Volumes/broad_dawnccle/processed_data/latest/mutagenesis_shortlist/breast_skin_top_mat.tsv", row.names = T, sep = "\t")
write_csv(breast_skin_top %>% filter(ExonID %in% rownames(final_breast_skin_top)), "/Volumes/broad_dawnccle/processed_data/latest/mutagenesis_shortlist/breast_skin_top_stats.csv")

############################################
##### Look at top sequences for kidney_skin ######
############################################
kidney_skin <- combined_psi_filtered %>% filter(grepl("kidney_skin", Folder))
kidney_skin_top <- kidney_skin %>% filter(abs(log2_PSI_ratio) > 2 | abs(log2_PSI_reverse_ratio) > 2 | abs(PSI_diff) > 0.3)

kidney_skin_top_mat <- plot_heatmap_for_ids(kidney_skin_top %>% pull(ExonID), all_sample_reps, prism_cluster)
# Get num of -1 for each row.
num_missing <- rowSums(kidney_skin_top_mat == -1)
# Get rows with at most 75 missing.
kidney_skin_top_mat_filtered <- kidney_skin_top_mat[num_missing <= 75,]
# Plot the heatmap
final_kidney_skin_top <- plot_heatmap_for_ids(rownames(kidney_skin_top_mat_filtered), all_sample_reps, prism_cluster, title = "Top seq in kidney_skin")
# Write the final top seq to a CSV file.
write.table(as.data.frame(final_kidney_skin_top), "/Volumes/broad_dawnccle/processed_data/latest/mutagenesis_shortlist/kidney_skin_top_mat.tsv", row.names = T, sep = "\t")
write_csv(kidney_skin_top %>% filter(ExonID %in% rownames(final_kidney_skin_top)), "/Volumes/broad_dawnccle/processed_data/latest/mutagenesis_shortlist/kidney_skin_top_stats.csv")






