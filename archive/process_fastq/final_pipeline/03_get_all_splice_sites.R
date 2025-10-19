library(tidyverse)
library(vroom)
library(data.table)
library(pheatmap)
library(preprocessCore)

all_files_df_path <- "U:/processed_data/latest/250131_merged_v3/WT_all_samples_raw_counts.csv"

all_files_df <- vroom(all_files_df_path)

# First get all perfect sequences.
perfect_PSI <- all_files_df %>% 
  filter((mode == "INCLUDED" & offset == "0:0:0")|(mode == "SKIPPED" & offset == "0")) %>% 
  group_by(index, sample) %>% 
  mutate(total_count = sum(count_scaled)) %>% 
  ungroup() %>% 
  filter(mode == "INCLUDED") %>% 
  mutate(skipped = total_count - count_scaled) %>%
  mutate(included = count_scaled) %>%
  select(-count, -count_scaled) 

perfect_PSI_pass_filter <- perfect_PSI %>% 
  filter(total_count > 30) %>% 
  mutate(PSI = included/total_count) 

# Get number of pass filter per sample:
num_pass_filter_by_sample <- perfect_PSI_pass_filter %>% 
  group_by(sample) %>% 
  summarise(num_pass = n())

# Plot heatmap of values:
PSI_mat <- perfect_PSI_pass_filter %>% 
  group_by(index, condition) %>% 
  summarise(PSI = mean(PSI)) %>% 
  pivot_wider(names_from = condition, values_from = PSI) %>% 
  column_to_rownames(var = "index") %>% 
  as.matrix()

# Quantile normalization:

PSI_mat <- normalize.quantiles(PSI_mat)
# Add back row and col names.
rownames(PSI_mat) <- rownames(PSI_mat)
colnames(PSI_mat) <- colnames(PSI_mat)
# Scale all columns between 0 and 1
PSI_mat <- apply(PSI_mat, 2, function(x) {
  min_x <- min(x, na.rm = TRUE)
  max_x <- max(x, na.rm = TRUE)
  if (max_x == min_x) {
    return(rep(0, length(x)))  # Avoid division by zero
  }
  (x - min_x) / (max_x - min_x)
})

# Filter to rows that has < 50% NA.
PSI_mat <- PSI_mat[rowSums(is.na(PSI_mat)) < ncol(PSI_mat)/2,]

# pheatmap(PSI_mat,
#          cluster_rows = T,
#          cluster_cols = F,
#          show_rownames = F,
#          show_colnames = T,
#          fontsize_col = 8,
#          main = "Perfect PSI values")

# > head(all_files_df)
# # A tibble: 6 × 7
# index       mode  offset count count_scaled sample
# <chr>       <chr> <chr>  <dbl>        <dbl> <chr> 
#   1 ENSG000000… INCL… 0:-1:0     1            1 769P-…
# 2 ENSG000000… INCL… 0:0:0     29           29 769P-…
# 3 ENSG000000… SKIP… 0          3            2 769P-…
# 4 ENSG000000… UNSP… 113        3            3 769P-…
# 5 ENSG000000… INCL… 0:-3:3     1            1 769P-…
# 6 ENSG000000… INCL… 0:0:0     91           91 769P-…
# # ℹ 1 more variable: condition <chr>
# Let's just look at 1 cell line for testing. 
HEK <- all_files_df %>% filter(condition == "HEK") %>% 
  group_by(index, sample) %>%
  mutate(total_count = sum(count_scaled)) %>% 
  filter(total_count > 30) %>% 
  mutate(percentage = count_scaled/total_count) %>%
  ungroup()

unspliced <- HEK %>% 
  filter(mode == "UNSPLICED")

# Filter for unspliced events
unspliced_df <- all_files_df %>% 
  filter(mode == "UNSPLICED") %>% 
  group_by(index, sample) %>%
  mutate(total_count = sum(count_scaled)) %>% 
  filter(total_count > 30) %>% 
  mutate(percentage = count_scaled / total_count) %>% 
  ungroup() 

# Compute max percentage per cell line and offset, keeping the index of the max sequence
unspliced_max_df <- unspliced_df %>% 
  group_by(sample, offset) %>% 
  slice_max(order_by = percentage, with_ties = FALSE) %>%  # Keep only the max
  ungroup() %>% 
  mutate(max_percentage = percentage * 100) %>%  # Convert to percentage
  filter(!(offset %in% c("0", "113", "144")))  

# Convert to wide format for heatmap (matrix of max % unspliced)
unspliced_matrix <- unspliced_max_df %>% 
  select(sample, offset, max_percentage, index) %>% 
  select(-index) %>% 
  pivot_wider(names_from = offset, values_from = max_percentage, values_fill = 0) %>% 
  column_to_rownames(var = "sample") %>% 
  as.matrix()

# order the columns by increasing offset
unspliced_matrix <- unspliced_matrix[, order(as.numeric(colnames(unspliced_matrix))), drop = FALSE]

# Plot heatmap of max percentage unspliced
pheatmap(unspliced_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 8,
         fontsize_col = 8,
         main = "Max % Unspliced by Offset")

# Compute the mean percentage per cell line and offset
unspliced_mean_df <- unspliced_df %>% 
  group_by(sample, offset) %>% 
  summarise(mean_percentage = mean(percentage * 100, na.rm = TRUE)) %>% 
  ungroup() %>% 
  filter(!(offset %in% c("0", "113", "144")))

# Convert to wide format for heatmap (matrix of mean % unspliced)
unspliced_mean_matrix <- unspliced_mean_df %>% 
  select(sample, offset, mean_percentage) %>% 
  pivot_wider(names_from = offset, values_from = mean_percentage, values_fill = 0) %>% 
  column_to_rownames(var = "sample") %>% 
  as.matrix()

# order the columns by increasing offset
unspliced_mean_matrix <- unspliced_mean_matrix[, order(as.numeric(colnames(unspliced_mean_matrix))), drop = FALSE]
# # Identify row indices where gaps should be inserted
# row_indices <- seq(3, nrow(unspliced_mean_matrix), by = 3)
# 
# # Create an empty row template
# empty_row <- matrix(NA, ncol = ncol(unspliced_mean_matrix), nrow = 1)
# colnames(empty_row) <- colnames(unspliced_mean_matrix)
# 
# # Insert empty rows into the matrix
# new_matrix <- unspliced_mean_matrix  # Start with original matrix
# 
# for (i in rev(row_indices)) {  # Reverse order to prevent index shifting
#   if (i == nrow(new_matrix)) {
#     new_matrix <- rbind(new_matrix, empty_row)  # Append empty row at the end
#   } else {
#     new_matrix <- rbind(new_matrix[1:i, ], empty_row, new_matrix[(i+1):nrow(new_matrix), , drop=FALSE])
#   }
# }
# 
# # Plot heatmap with row gaps
# pheatmap(new_matrix,
#          cluster_rows = FALSE,
#          cluster_cols = FALSE,
#          show_rownames = TRUE,
#          show_colnames = TRUE,
#          fontsize_row = 8,
#          fontsize_col = 8,
#          main = "Mean % Unspliced by Offset",
#          na_col = "white")  # Make NA values appear as gaps

pheatmap(unspliced_mean_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 8,
         fontsize_col = 8,
         main = "Mean % Unspliced by Offset")

HEK_included <- HEK %>% 
  filter(mode == "INCLUDED") %>% 
  filter(offset != "0:0:0") %>% 
  filter(count_scaled > 5 & percentage > 0.05) %>% 
  # offset doesn't contain 9999
  filter(!grepl("9999", offset)) 

# Get all these sequences included and 0:0:0 or skipped 0
HEK_included_controls <- HEK %>% 
  filter((mode == "INCLUDED" & offset == "0:0:0")|(mode == "SKIPPED" & offset == "0")) %>% 
  filter(index %in% HEK_included$index) 

# Concate the two dataframes
HEK_included_controls <- rbind(HEK_included_controls, HEK_included)
