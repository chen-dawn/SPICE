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

all_sample_reps <- fread("/Volumes/broad_dawnccle/processed_data/missplicing_processed_df/V5/all_sample_reps_PSI.csv")
all_samples_wide <- all_sample_reps %>%
  mutate(PSI = count/(count + skipped)) %>%
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
cellline_metadata <- fread("/Volumes/broad_dawnccle/for_anisha/cellline_data_full_metadata.csv") %>% 
  select(DepMap_ID, StrippedName) %>% distinct()
# Rename V1 in gex based on cellline_metadata
gex_formatted <- gex %>% 
  left_join(cellline_metadata, by = "DepMap_ID") %>% 
  select(-DepMap_ID) %>% 
  rename(StrippedName = "condition") %>%
  filter(condition %in% c(colnames(all_samples_mat), "K562", "8MGBA", "A375", "SKNAS"))

#convert to matrix no condition column
gex_mat <- as.matrix(gex_formatted %>% select(-condition))
rownames(gex_mat) <- gex_formatted$condition

# Get rownames order of gex_mat
gex_mat_order <- rownames(gex_mat)

# Order the all_samples_mat based on the gex_mat_order
all_samples_mat_aligned <- all_samples_mat[,gex_mat_order]

# Subset matrices based on the common names
gex_mat_aligned <- gex_mat

# Initialize a matrix to store the correlation values
correlation_matrix <- matrix(NA, nrow = nrow(all_samples_mat_aligned), ncol = ncol(gex_mat_aligned))
rownames(correlation_matrix) <- rownames(all_samples_mat_aligned)
colnames(correlation_matrix) <- colnames(gex_mat_aligned)

# Calculate the correlation for each element-gene pair
all_samples_mat_aligned_t <- t(all_samples_mat_aligned)
for (i in 1:ncol(all_samples_mat_aligned_t)) {
  for (j in 1:ncol(gex_mat_aligned)) {
    element_values <- all_samples_mat_aligned_t[, i]
    gene_values <- gex_mat_aligned[, j]
    
    # df <- data.frame(element_values, gene_values)
    # Check if there are enough complete pairs to compute the correlation
    if (sum(complete.cases(element_values, gene_values)) > 5) {
      # Check for zero standard deviation
      if (sd(element_values, na.rm = TRUE) != 0 && sd(gene_values, na.rm = TRUE) != 0) {
        # Compute correlation, use method = "pearson" by default
        correlation_matrix[i, j] <- cor(element_values, gene_values, method = "pearson", use = "complete.obs")
      } else {
        correlation_matrix[i, j] <- NA  # Set to NA if standard deviation is zero
      }
    } else {
      correlation_matrix[i, j] <- NA  # Set to NA if not enough complete pairs
    }
  }
}

# Remove rows with >80% NA values
correlation_matrix <- correlation_matrix[rowSums(is.na(correlation_matrix)) <= 0.8 * ncol(correlation_matrix), ]


# Create a heatmap of the correlation matrix
pheatmap(correlation_matrix, cluster_rows = TRUE, 
         cluster_cols = TRUE, show_rownames = F, 
         show_colnames = T, fontsize = 6)

# Filter for rows and columns where correlation > 0.8.
indices <- which(correlation_matrix > 0.8, arr.ind = TRUE)
# Extract unique row and column names
rows_to_keep <- unique(rownames(correlation_matrix)[indices[, 1]])
cols_to_keep <- unique(colnames(correlation_matrix)[indices[, 2]])

# Subset the correlation matrix
subset_correlation_matrix <- correlation_matrix[rows_to_keep, cols_to_keep, drop = FALSE]
pheatmap(subset_correlation_matrix, cluster_rows = TRUE, 
         cluster_cols = TRUE, show_rownames = F, 
         show_colnames = T, fontsize = 4, main = "RBP to PSI correlation > 0.8")

# now only subset to the shortlist.
seq_shortlist <- c("ENSG00000138821.14;SLC39A8;chr4-102285948-102285973-102267871-102268079-102304316-102304481", 
                   "ENSG00000135365.16;PHF21A;chr11-45946075-45946098-45938156-45938312-45948885-45948946")
subset_correlation_matrix <- correlation_matrix[rownames(correlation_matrix) %in% seq_shortlist, ]
color_palette <- colorRampPalette(c("#45abd7", "white", "#cf1c1a"))(100)

# Set the breaks for the color scale
breaks <- seq(-1, 1, length.out = 101)

# Generate the heatmap with specified parameters
pheatmap(subset_correlation_matrix,
         cluster_rows = FALSE, 
         cluster_cols = TRUE, 
         show_rownames = FALSE, 
         show_colnames = TRUE, 
         fontsize = 6, 
         main = "High in Kelly",
         breaks = breaks)

# Plot scatter plot for KHSRP gene.
gex_KHSRP <- gex_mat_aligned[, "KHSRP"]
psi_KHSRP <- all_samples_mat_aligned["ENSG00000138821.14;SLC39A8;chr4-102285948-102285973-102267871-102268079-102304316-102304481", ]
plot_df <- data.frame(gex_KHSRP, psi_KHSRP, cell_line = rownames(gex_mat_aligned))
ggplot(plot_df, aes(x = gex_KHSRP, y = psi_KHSRP, color = cell_line)) +
  geom_point() + 
  theme_minimal() +
  labs(x = "KHSRP gene expression", y = "PSI")


element <- c("ENSG00000138821.14;SLC39A8;chr4-102285948-102285973-102267871-102268079-102304316-102304481", 
             "ENSG00000135365.16;PHF21A;chr11-45946075-45946098-45938156-45938312-45948885-45948946")
test <- all_samples_wide %>% filter(index %in% element)

# Extract PSI values for the specific element
psi_values <- all_samples_mat[element, ]
pheatmap(psi_values, cluster_rows = F, cluster_cols = F, show_rownames = T, show_colnames = T, main = "PSI values for KHSRP gene")

# Subset to the first 50 genes
subset_genes <- colnames(gex_mat_aligned)[1:20]
gex_subset <- gex_mat_aligned[, subset_genes]

# Create a data frame for plotting
plot_df <- data.frame(
  psi_values = rep(psi_values, times = ncol(gex_subset)),
  gene_expression = as.vector(as.matrix(gex_subset)),
  gene = rep(subset_genes, each = nrow(gex_subset)),
  cell_line = rep(rownames(gex_subset), times = ncol(gex_subset))
)
# Plot using ggplot with facets
ggplot(plot_df, aes(x = gene_expression, y = psi_values)) +
  geom_point(aes(color = cell_line)) + 
  geom_label(aes(label = cell_line), box.padding = 0.5, point.padding = 0.5) +
  facet_wrap(~ gene, scales = "free") +
  theme_minimal() +
  labs(x = "Gene Expression", y = "PSI", title = "Scatter Plots of PSI vs Gene Expression")

