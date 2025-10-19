# Compare RBP gene expression with PSI.

library(tidyverse)
library(vroom)
library(data.table)
library(Biostrings)
library(ggpointdensity)
library(pheatmap)

reverse_complement <- function(dna_seq) {
  complement <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
  nucleotides <- unlist(strsplit(dna_seq, ""))
  complement_nucleotides <- complement[nucleotides]
  reverse_complement_seq <- paste(rev(complement_nucleotides), collapse = "")
  return(reverse_complement_seq)
}


###### Read in all CCLE gene expression ######
all_gex <- fread("~/Dropbox (Harvard University)/02Splicing/latest/CCLE_expression.csv", sep = ",")
# For each variable name it looks like TSPAN6 (7105). Only keep the gene name not in the brackets. 
all_gex <- all_gex %>% 
  rename_with(~str_extract(., "^[^\\s]+")) %>% 
  rename(V1 = "DepMap_ID") 

all_predicted_RBP <- fread("~/Dropbox (Harvard University)/02Splicing/latest/RBP_possible_Gerstberger_NatRev2014.csv")

# Get Gene names that are in colnames of all_gex and all_predicted_RBP
all_available_RBP <- intersect(colnames(all_gex), all_predicted_RBP$`gene name`)

# Filter to only the predicted RBPs
all_gex_RBP <- all_gex %>% select(DepMap_ID, all_of(all_available_RBP))


all_sample_reps <- fread("~/Dropbox (Harvard University)/02Splicing/latest/all_sample_reps_PSI.csv")
# all_sample_reps <- fread("/Volumes/broad_dawnccle/processed_data/latest/all_sample_reps_PSI.csv")
all_samples_wide <- all_sample_reps %>%
  # mutate(PSI = count/(count + skipped)) %>%
  mutate(PSI = log2((count+1)/(skipped+1))) %>%
  filter((count + skipped) > 30) %>%
  mutate(condition = toupper(str_extract(condition, "^[^_-]+"))) %>% 
  group_by(condition, index) %>%
  summarise(PSI = mean(PSI)) %>%
  ungroup() %>%
  pivot_wider(names_from = condition, values_from = PSI, values_fill = NA)

all_samples_mat <- as.matrix(all_samples_wide[, -1])
rownames(all_samples_mat) <- all_samples_wide$index

# Get gene expression
gex <- fread("/Volumes/broad_dawnccle/for_anisha/CCLE_expression.RBPs.DEDUPLICATED.csv") %>% 
  rename(V1 = "DepMap_ID") 
# Arrange gex colnames roughly by all_gex_RBP colnames
existing_RBPs <- intersect(colnames(all_gex_RBP), colnames(gex))
gex <- gex %>% select(DepMap_ID, all_of(existing_RBPs)) %>% 
  rename_with(~str_extract(., "^[^\\s]+")) %>% 
  rename(DepMap_ID = "V1")

cellline_metadata <- fread("/Volumes/broad_dawnccle/for_anisha/cellline_data_full_metadata.csv") %>% 
  select(DepMap_ID, StrippedName) %>% distinct()
# Rename V1 in gex based on cellline_metadata
gex_formatted <- all_gex_RBP %>% 
  left_join(cellline_metadata, by = "DepMap_ID") %>% 
  select(-DepMap_ID) %>% 
  rename(StrippedName = "condition") %>%
  filter(condition %in% c(colnames(all_samples_mat), "K562", "8MGBA", "A375", "SKNAS"))


# Get common cell lines in gex and all_samples_mat
common_cell_lines <- intersect(colnames(all_samples_mat), gex_formatted$condition)
#convert to matrix no condition column
gex_mat <- as.matrix(gex_formatted %>% select(-condition))
rownames(gex_mat) <- gex_formatted$condition

# Order the all_samples_mat based on the names
all_samples_mat_aligned <- all_samples_mat[, common_cell_lines]

# Subset matrices based on the common names
gex_mat_aligned <- gex_mat[common_cell_lines,]

correlation_matrix <- cor(t(all_samples_mat_aligned), gex_mat_aligned, use = "pairwise.complete.obs")

# # Initialize a matrix to store the correlation values
# correlation_matrix <- matrix(NA, nrow = nrow(all_samples_mat_aligned), ncol = ncol(gex_mat_aligned))
# rownames(correlation_matrix) <- rownames(all_samples_mat_aligned)
# colnames(correlation_matrix) <- colnames(gex_mat_aligned)
# 
# # Calculate the correlation for each element-gene pair
# all_samples_mat_aligned_t <- t(all_samples_mat_aligned)
# for (i in 1:ncol(all_samples_mat_aligned_t)) {
#   for (j in 1:ncol(gex_mat_aligned)) {
#     element_values <- all_samples_mat_aligned_t[, i]
#     gene_values <- gex_mat_aligned[, j]
#     
#     # df <- data.frame(element_values, gene_values)
#     # Check if there are enough complete pairs to compute the correlation
#     if (sum(complete.cases(element_values, gene_values)) > 5) {
#       # Check for zero standard deviation
#       if (sd(element_values, na.rm = TRUE) != 0 && sd(gene_values, na.rm = TRUE) != 0) {
#         # Compute correlation, use method = "pearson" by default
#         correlation_matrix[i, j] <- cor(element_values, gene_values, method = "pearson", use = "complete.obs")
#       } else {
#         correlation_matrix[i, j] <- NA  # Set to NA if standard deviation is zero
#       }
#     } else {
#       correlation_matrix[i, j] <- NA  # Set to NA if not enough complete pairs
#     }
#   }
# }
# 
# # Remove rows with >80% NA values
# correlation_matrix <- correlation_matrix[rowSums(is.na(correlation_matrix)) <= 0.8 * ncol(correlation_matrix), ]
# 
# 
# # Create a heatmap of the correlation matrix
# # pheatmap(correlation_matrix, cluster_rows = TRUE, 
# #          cluster_cols = TRUE, show_rownames = F, 
# #          show_colnames = T, fontsize = 6)
# 
# # Filter for rows and columns where correlation > 0.8.
# indices <- which(abs(correlation_matrix) > 0.8, arr.ind = TRUE)
# # Extract unique row and column names
# rows_to_keep <- unique(rownames(correlation_matrix)[indices[, 1]])
# cols_to_keep <- unique(colnames(correlation_matrix)[indices[, 2]])
# 
# # Subset the correlation matrix
# subset_correlation_matrix <- correlation_matrix[rows_to_keep, cols_to_keep, drop = FALSE]
# pheatmap(subset_correlation_matrix, cluster_rows = TRUE, 
#          cluster_cols = TRUE, show_rownames = F, 
#          show_colnames = T, fontsize = 6, main = "RBP to PSI correlation > 0.8")
# 
# # now only subset to the shortlist.
# seq_shortlist <- c("ENSG00000138821.14;SLC39A8;chr4-102285948-102285973-102267871-102268079-102304316-102304481", 
#                    "ENSG00000135365.16;PHF21A;chr11-45946075-45946098-45938156-45938312-45948885-45948946")
# subset_correlation_matrix <- correlation_matrix[rownames(correlation_matrix) %in% seq_shortlist, ]
# #subset values >abs 0.6. for values
# 
# subset_correlation_matrix <- subset_correlation_matrix[, apply(abs(subset_correlation_matrix), 2, function(x) any(x > 0.3))]
# color_palette <- colorRampPalette(c("#45abd7", "white", "#cf1c1a"))(100)
# 
# # Set the breaks for the color scale
# breaks <- seq(-1, 1, length.out = 101)
# 
# # Generate the heatmap with specified parameters
# pheatmap(subset_correlation_matrix,
#          cluster_rows = FALSE, 
#          cluster_cols = TRUE, 
#          show_rownames = FALSE, 
#          show_colnames = TRUE, 
#          fontsize = 7, 
#          main = "High in Kelly",
#          breaks = breaks)


# Plot scatter plot for KHSRP gene.
gex_KHSRP <- gex_mat_aligned[, "RBFOX2"]
psi_KHSRP <- all_samples_mat_aligned["ENSG00000138821.14;SLC39A8;chr4-102285948-102285973-102267871-102268079-102304316-102304481", ]
plot_df <- data.frame(gex_KHSRP, psi_KHSRP, cell_line = rownames(gex_mat_aligned))
ggplot(plot_df, aes(x = gex_KHSRP, y = psi_KHSRP, color = cell_line)) +
  geom_point() + 
  theme_minimal() +
  labs(x = "RFCOL5 gene expression", y = "PSI")

# Plot scatter for RBFOX2
gene_name <- "RBFOX2"
seq_name <- "ENSG00000135365.16;PHF21A;chr11-45946075-45946098-45938156-45938312-45948885-45948946"
gex_RBFOX2 <- gex_mat_aligned[, gene_name]
psi_RBFOX2 <- all_samples_mat_aligned[seq_name, ]
corr_val <- cor(gex_RBFOX2, psi_RBFOX2, use = "pairwise.complete.obs", method = "pearson")
plot_df <- data.frame(gex_RBFOX2, psi_RBFOX2, cell_line = rownames(gex_mat_aligned))
ggplot(plot_df, aes(x = gex_RBFOX2, y = psi_RBFOX2, color = cell_line)) +
  geom_point() + 
  theme_minimal() +
  labs(x = "Gene expression", y = "PSI") + 
  ggtitle(paste(gene_name, seq_name)) + 
  # Add subtitle with the gene name and sequence name
  theme(plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, face = "italic")) + 
  # Add correlation as subtitle text
  labs(subtitle = paste("Correlation:", round(corr_val, 3)))

# Plot scatter for RBFOX2
gex_RBFOX2 <- gex_mat_aligned[, "RBFOX2"]
psi_RBFOX2 <- all_samples_mat_aligned["ENSG00000135365.16;PHF21A;chr11-45946075-45946098-45938156-45938312-45948885-45948946", ]
plot_df <- data.frame(gex_RBFOX2, psi_RBFOX2, cell_line = rownames(gex_mat_aligned))
ggplot(plot_df, aes(x = gex_RBFOX2, y = psi_RBFOX2, color = cell_line)) +
  geom_point() + 
  theme_minimal() +
  labs(x = "RBFOX2 gene expression", y = "PSI")


############# Look at Tau and PSI #################
# Calculate Tau for each gene.
calculate_tau <- function(row){
  # Remove NA values from the row
  non_na_row <- row[!is.na(row)]
  # If the row is empty after removing NAs, return NA
  if (length(non_na_row) < 20) {
    return(NA)
  }
  # Normalize the row by dividing by the max value of the non-NA row
  norm_row <- non_na_row / max(non_na_row)
  # Calculate tau using the number of non-NA values
  tau <- sum(1 - norm_row) / (length(non_na_row) - 1)
  return(tau)
}
  
# upsilon is the metric for PSI. We add 1 so that values close to 0 will not be inflated.
calculate_upsilon <- function(row) {
  # Remove NA values from the row
  non_na_row <- row[!is.na(row)]
  # If the row is empty after removing NAs, return NA
  if (length(non_na_row) < 20) {
    return(NA)
  }
  # Add 1 to every value in the row.
  non_na_row <- non_na_row + 1
  # Normalize the row by dividing by the max value of the non-NA row
  norm_row <- non_na_row / max(non_na_row)
  # Calculate tau using the number of non-NA values
  tau <- sum(1 - norm_row) / (length(non_na_row) - 1)
  return(tau * 2)
}

# Calculate tau for the gex values. 
tau_values <- apply(gex_mat_aligned, 2, calculate_tau)
# calculate dynamic range.
dynamic_range <- apply(gex_mat_aligned, 2, function(x) max(x, na.rm =T) - min(x, na.rm = T))
# convert the values to a df.
tau_df <- data.frame(gene = names(tau_values), tau = tau_values, dynamic_range = dynamic_range)
# Shortlist to tau >0.5 in the gex_mat_aligned.
shortlist_genes <- tau_df %>% filter(tau > 0.5 & dynamic_range > 2) %>% pull(gene)
# Subset the gex_mat_aligned to the shortlist_genes.
gex_mat_aligned_shortlist <- gex_mat_aligned[, c(shortlist_genes)]
# Plot heatmap for the shortlist_genes and add the tau values.
pheatmap(gex_mat_aligned_shortlist, cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = TRUE, show_colnames = TRUE, fontsize = 6,)

# Calculate upsilon for the psi values.
# Scale the all_samples_mat_aligned so that all values are between 0 and 1
upsilon_values <- apply(all_samples_mat_aligned + 11.85651, 1, calculate_upsilon)
# convert the values to a df.
upsilon_df <- data.frame(gene = names(upsilon_values), upsilon = upsilon_values)
# Plot histogram of upsilon values.
ggplot(upsilon_df, aes(x = upsilon)) +
  geom_histogram(bins = 30, fill = "blue", color = "black") +
  labs(title = "Histogram of Upsilon values for PSI", x = "Upsilon", y = "Frequency")

# Also calculate reverse upsilon.
reverse_upsilon_values <- apply(-(all_samples_mat_aligned + 11.85651) + 20.575, 1, calculate_upsilon)
# convert the values to a df.
reverse_upsilon_df <- data.frame(gene = names(reverse_upsilon_values), upsilon = reverse_upsilon_values)
# Plot histogram of upsilon values.
ggplot(reverse_upsilon_df, aes(x = upsilon)) +
  geom_histogram(bins = 30, fill = "blue", color = "black") +
  labs(title = "Histogram of Upsilon values for PSI", x = "Upsilon", y = "Frequency")

# Shortlist to >0.6
shortlist_psi <- upsilon_df %>% filter(upsilon > 0.7) %>% pull(gene)
shortlist_pis_reverse <- reverse_upsilon_df %>% filter(upsilon > 0.9) %>% pull(gene)
# merge
shortlist_psi <- c(shortlist_psi, shortlist_pis_reverse)
# Subset the all_samples_mat_aligned to the shortlist_psi.
all_samples_mat_aligned_shortlist <- all_samples_mat_aligned[shortlist_psi, ]
# Plot heatmap for the shortlist_psi.
pheatmap(all_samples_mat_aligned_shortlist, cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = TRUE, show_colnames = TRUE, fontsize = 6,)

# For each of these sequences we will plot the gene expression against the PSI values. And write to output.
outdir <- "~/Dropbox (Harvard University)/02Splicing/test_psi_log_scaled_high_corr/"
# Loop through the shortlisted genes and create scatter plots for each gene's expression against PSI.
# Do the same for the shortlisted PSI values
for (seq in shortlist_psi) {
  for (gene in colnames(gex_mat_aligned_shortlist)) {
    gex_values <- gex_mat_aligned_shortlist[, gene]
    psi_values <- all_samples_mat_aligned_shortlist[seq, ]
    
    # Create a dataframe for plotting
    plot_df <- data.frame(gex_values, psi_values, cell_line = rownames(gex_mat_aligned_shortlist))
    corr_val <- cor(gex_values, psi_values, use = "pairwise.complete.obs", method = "pearson")
    if (is.na(corr_val)) {
      # skip
      next
    }
    if (abs(corr_val) < 0.5){
      # skip
      next
    }
    
    # Generate the scatter plot
    p <- ggplot(plot_df, aes(x = gex_values, y = psi_values, color = cell_line)) +
      geom_point() +
      theme_minimal() +
      labs(x = "Gene expression", y = paste(seq, "PSI")) +
      ggtitle(paste(gene, seq)) + 
      # Add subtitle with the gene name and sequence name
      theme(plot.title = element_text(size = 10, face = "bold"),
            plot.subtitle = element_text(size = 8, face = "italic")) + 
      # Add correlation as subtitle text
      labs(subtitle = paste("Correlation:", round(corr_val, 3)))
    
    # Save the plot to a file
    ggsave(filename = paste0(outdir, gene, "_", seq, "_scatter_plot.png"), plot = p, width = 8, height = 6)
  }
}

outdir <- "~/Dropbox (Harvard University)/02Splicing/test_psi_log_scaled_high_corr_seq_first/"
# Loop through the shortlisted genes and create scatter plots for each gene's expression against PSI.
# Do the same for the shortlisted PSI values
for (seq in shortlist_psi) {
  for (gene in colnames(gex_mat_aligned_shortlist)) {
    gex_values <- gex_mat_aligned_shortlist[, gene]
    psi_values <- all_samples_mat_aligned_shortlist[seq, ]
    
    # Create a dataframe for plotting
    plot_df <- data.frame(gex_values, psi_values, cell_line = rownames(gex_mat_aligned_shortlist))
    corr_val <- cor(gex_values, psi_values, use = "pairwise.complete.obs", method = "pearson")
    if (is.na(corr_val)) {
      # skip
      next
    }
    if (abs(corr_val) < 0.5){
      # skip
      next
    }
    
    # Generate the scatter plot
    p <- ggplot(plot_df, aes(x = gex_values, y = psi_values, color = cell_line)) +
      geom_point() +
      theme_minimal() +
      labs(x = "Gene expression", y = paste(seq, "PSI")) +
      ggtitle(paste(gene, seq)) + 
      # Add subtitle with the gene name and sequence name
      theme(plot.title = element_text(size = 10, face = "bold"),
            plot.subtitle = element_text(size = 8, face = "italic")) + 
      # Add correlation as subtitle text
      labs(subtitle = paste("Correlation:", round(corr_val, 3)))
    
    # Save the plot to a file
    ggsave(filename = paste0(outdir, seq, "_", gene, "_scatter_plot.png"), plot = p, width = 8, height = 6)
  }
}

outdir <- "~/Dropbox (Harvard University)/02Splicing/test_psi_log_scaled_dumb_max_metric/"
# Get only the "shortlist".
# for (seq in shortlist_psi) {
#   for (gene in colnames(gex_mat_aligned_shortlist)) {
#     gex_values <- gex_mat_aligned_shortlist[, gene]
#     psi_values <- all_samples_mat_aligned_shortlist[seq, ]
#     
#     # Get cell line name with the highest PSI value.
#     cell_line_max_psi <- names(psi_values)[which.max(psi_values)]
#     
#     # Get the cell line with the highest gene expression value
#     cell_line_max_gex <- names(gex_values)[which.max(gex_values)]
#     
#     # Check if the cell line with the highest PSI value also has the highest gene expression
#     if (cell_line_max_psi == cell_line_max_gex) {
#       plot_df <- data.frame(gex_values, psi_values, cell_line = rownames(gex_mat_aligned_shortlist))
#       corr_val <- cor(gex_values, psi_values, use = "pairwise.complete.obs", method = "pearson")
#       
#       # Generate the scatter plot
#       p <- ggplot(plot_df, aes(x = gex_values, y = psi_values, color = cell_line)) +
#         geom_point() +
#         theme_minimal() +
#         labs(x = "Gene expression", y = paste(seq, "PSI")) +
#         ggtitle(paste(gene, seq)) + 
#         # Add subtitle with the gene name and sequence name
#         theme(plot.title = element_text(size = 10, face = "bold"),
#               plot.subtitle = element_text(size = 8, face = "italic")) + 
#         # Add correlation as subtitle text
#         labs(subtitle = paste("Correlation:", round(corr_val, 3)))
#       # Save the plot to a file
#       ggsave(filename = paste0(outdir, gene, "_", seq, "_scatter_plot_shortlist.png"), plot = p, width = 8, height = 6)
#     }
#   }
# }

# Also get for reverse upsilon.
shortlist_psi <- reverse_upsilon_df %>% filter(upsilon > 0.6) %>% pull(gene)
# Subset the all_samples_mat_aligned to the shortlist_psi.
all_samples_mat_aligned_shortlist <- all_samples_mat_aligned[shortlist_psi, ]
# Plot heatmap for the shortlist_psi.
pheatmap(all_samples_mat_aligned_shortlist, cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = TRUE, show_colnames = TRUE, fontsize = 6,)

# for (seq in shortlist_psi) {
#   for (gene in colnames(gex_mat_aligned_shortlist)) {
#     gex_values <- gex_mat_aligned_shortlist[, gene]
#     psi_values <- all_samples_mat_aligned_shortlist[seq, ]
#     
#     # Get cell line name with the highest PSI value.
#     cell_line_max_psi <- names(psi_values)[which.min(psi_values)]
#     
#     # Get the cell line with the highest gene expression value
#     cell_line_max_gex <- names(gex_values)[which.max(gex_values)]
#     
#     # Check if the cell line with the highest PSI value also has the highest gene expression
#     if (cell_line_max_psi == cell_line_max_gex) {
#       print(paste("Gene:", gene, "Sequence:", seq, "Cell line:", cell_line_max_psi, "has the highest PSI and GEX"))
#       
#       # Create a dataframe for plotting
#       plot_df <- data.frame(gex_values, psi_values, cell_line = rownames(gex_mat_aligned_shortlist))
#       
#       # Generate the scatter plot
#       p <- ggplot(plot_df, aes(x = gex_values, y = psi_values, color = cell_line)) +
#         geom_point() +
#         theme_minimal() +
#         labs(x = "Gene expression", y = paste(seq, "PSI")) +
#         ggtitle(paste("Scatter plot for", seq, "and", gene))
#       
#       # Save the plot to a file
#       ggsave(filename = paste0(outdir, gene, "_", seq, "_scatter_plot_shortlist.png"), plot = p, width = 8, height = 6)
#     }
#   }
# }


##### We try to look at the sequences that are significant based on one-vs-all #####
# Read in the significant sequences
combined_psi <- fread("~/Dropbox (Harvard University)/02Splicing/latest/rmats_one_vs_all_combined_output_PSI.tsv")
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


combined_psi_filtered <- combined_psi %>% 
  filter(count_sum_average1 >30) %>% 
  filter(count_sum_average2 > 30) %>% 
  mutate(log2_PSI_ratio = log2(PSI_ratio), log2_PSI_reverse_ratio = log2(PSI_reverse_ratio))

test2 <- combined_psi %>% filter(ExonID %in% c("ENSG00000138821.14;SLC39A8;chr4-102285948-102285973-102267871-102268079-102304316-102304481", "ENSG00000135365.16;PHF21A;chr11-45946075-45946098-45938156-45938312-45948885-45948946"))

top_seq <- combined_psi_filtered %>% filter(abs(log2_PSI_ratio) > 2 | abs(log2_PSI_reverse_ratio) > 2)
num_obs_per_seq <- top_seq %>% group_by(ExonID) %>% summarise(num_obs = n())
top_seq <- top_seq %>% left_join(num_obs_per_seq, by = "ExonID") %>% filter(num_obs > 1)


outdir <- "~/Dropbox (Harvard University)/02Splicing/test_psi_high_corr_rmats_sig_seq_first/"
# Loop through the shortlisted genes and create scatter plots for each gene's expression against PSI.
# Do the same for the shortlisted PSI values
# for (seq in unique(top_seq$ExonID)) {
#   for (gene in colnames(gex_mat_aligned_shortlist)) {
#     gex_values <- gex_mat_aligned_shortlist[, gene]
#     psi_values <- all_samples_mat_aligned[seq, ]
#     
#     # Create a dataframe for plotting
#     plot_df <- data.frame(gex_values, psi_values, cell_line = rownames(gex_mat_aligned_shortlist))
#     corr_val <- cor(gex_values, psi_values, use = "pairwise.complete.obs", method = "pearson")
#     if (is.na(corr_val)) {
#       # skip
#       next
#     }
#     if (abs(corr_val) < 0.5){
#       # skip
#       next
#     }
#     
#     # Generate the scatter plot
#     p <- ggplot(plot_df, aes(x = gex_values, y = psi_values, color = cell_line)) +
#       geom_point() +
#       theme_minimal() +
#       labs(x = "Gene expression", y = paste(seq, "PSI")) +
#       ggtitle(paste(gene, seq)) + 
#       # Add subtitle with the gene name and sequence name
#       theme(plot.title = element_text(size = 10, face = "bold"),
#             plot.subtitle = element_text(size = 8, face = "italic")) + 
#       # Add correlation as subtitle text
#       labs(subtitle = paste("Correlation:", round(corr_val, 3))) + 
#       # Set the X and Y limits
#       xlim(c(0, 8)) + ylim(c(-5, 5))
#     
#     
#     # Save the plot to a file
#     ggsave(filename = paste0(outdir, gene, "_", seq, "_scatter_plot.png"), plot = p, width = 8, height = 6)
#   }
# }

# Extract gene and exon IDs
exon_list <- unique(top_seq$ExonID)
gene_list <- colnames(gex_mat_aligned_shortlist)

# Function to calculate correlations for a given exon
calculate_correlation_for_exon <- function(seq) {
  # Get PSI values for the given exon
  psi_values <- all_samples_mat_aligned[seq, ]
  
  # Calculate correlations with all genes
  correlations <- map_dfr(gene_list, function(gene) {
    gex_values <- gex_mat_aligned_shortlist[, gene]
    corr_val <- cor(gex_values, psi_values, use = "pairwise.complete.obs", method = "pearson")
    tibble(Gene = gene, Exon = seq, Correlation = corr_val)
  })
  
  return(correlations)
}

# Calculate correlations for all exons and combine results
correlation_df <- map_dfr(exon_list, calculate_correlation_for_exon)
# get top 20 absolute values
top_correlations <- correlation_df %>% arrange(desc(abs(Correlation))) %>% head(20)
high_correlations <- correlation_df %>% filter(abs(Correlation) > 0.5) %>% group_by(Gene) %>% summarise(num_high_correlations = n()) %>% arrange(desc(num_high_correlations))

# Plot scatter plot for the top 20 correlations
outdir <- "~/Dropbox (Harvard University)/02Splicing/test_psi_high_corr_rmats_sig_seq_first/"
for (i in 1:nrow(top_correlations)) {
  gene <- top_correlations$Gene[i]
  seq <- top_correlations$Exon[i]
  
  gex_values <- gex_mat_aligned_shortlist[, gene]
  psi_values <- all_samples_mat_aligned[seq, ]
  
  # Create a dataframe for plotting
  plot_df <- data.frame(gex_values, psi_values, cell_line = rownames(gex_mat_aligned_shortlist))
  corr_val <- cor(gex_values, psi_values, use = "pairwise.complete.obs", method = "pearson")
  
  # Generate the scatter plot
  p <- ggplot(plot_df, aes(x = gex_values, y = psi_values, color = cell_line)) +
    geom_point() +
    theme_minimal() +
    labs(x = "Gene expression", y = paste(seq, "PSI")) +
    ggtitle(paste(gene, seq)) + 
    # Add subtitle with the gene name and sequence name
    theme(plot.title = element_text(size = 10, face = "bold"),
          plot.subtitle = element_text(size = 8, face = "italic")) + 
    # Add correlation as subtitle text
    labs(subtitle = paste("Correlation:", round(corr_val, 3))) + 
    # Set the X and Y limits
    xlim(c(0, 8)) + ylim(c(-5, 5))
  
  # Save the plot to a file
  ggsave(filename = paste0(outdir,"00TOP_", gene, "_", seq, "_scatter_plot.png"), plot = p, width = 8, height = 6)
}


# # What if I calculate tau from the non-normalized values? We will do a reverse log2 transform.
# # Not going to do this anymore looks interesting...
# gex_TPM <- 2^gex_mat_aligned - 1
# tau_values_TPM <- apply(gex_TPM, 2, calculate_tau)
# # convert the values to a df.
# tau_df_TPM <- data.frame(gene = names(tau_values_TPM), tau = tau_values_TPM)
# # Plot histogram of tau values.
# ggplot(tau_df_TPM, aes(x = tau)) +
#   geom_histogram(bins = 30, fill = "blue", color = "black") +
#   labs(title = "Histogram of Tau values for TPM", x = "Tau", y = "Frequency")
# # Shortlist to tau >0.5 in the gex_mat_aligned.
# shortlist_genes_TPM <- tau_df_TPM %>% filter(tau > 0.8) %>% pull(gene)
# # Subset the gex_mat_aligned to the shortlist_genes.
# gex_mat_aligned_shortlist_TPM <- gex_mat_aligned[, c(shortlist_genes_TPM)]
# # Plot heatmap for the shortlist_genes and add the tau values.
# pheatmap(gex_mat_aligned_shortlist_TPM, cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = TRUE, show_colnames = TRUE, fontsize = 6,)

# # Filter to complete rows of all_samples_mat_aligned.
# impute_mean <- function(x) {
#   x[is.na(x)] <- mean(x, na.rm = TRUE)
#   return(x)
# }
# 
# all_samples_mat_aligned_clean <- apply(all_samples_mat_aligned, 2, impute_mean)
# 
# all_samples_mat_aligned_clean <- all_samples_mat_aligned[complete.cases(all_samples_mat_aligned), ]
# 
# # Perform Canonical Correlation Analysis with cleaned data
# require(CCA)
# require(CCP)
# cca_result <- cc(gex_mat_aligned, t(all_samples_mat_aligned_clean)+0.01)
# 
# # View results
# print(cca_result)
# 
# # Heatmap of cca_result$xcoef
# pheatmap(cca_result$xcoef, cluster_rows = T, cluster_cols = T, show_rownames = T, show_colnames = T, fontsize = 6)
