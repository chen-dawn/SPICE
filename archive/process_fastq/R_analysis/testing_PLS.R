
gex <- matrix(c(2,2,2,7,2,3,4,5,5,2,2,2,1,1,1,1,2,2,2,2), nrow=5, byrow=TRUE)
rownames(gex) <- c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5")
colnames(gex) <- c("Cell1", "Cell2", "Cell3", "Cell4")

PSI <- matrix(c(0.1,0.2,0.1,0.7,0.2,0.3,0.4,0.5,0.1,0.1,0.1,0.1), nrow=3, byrow=TRUE)
rownames(PSI) <- c("Seq1", "Seq2", "Seq3")
colnames(PSI) <- c("Cell1", "Cell2", "Cell3", "Cell4")

# Create a larger gex matrix with 10 genes and 30 cells
set.seed(1)  # For reproducibility
gex <- matrix(sample(1:10, 300, replace = TRUE), nrow = 10)
rownames(gex) <- paste0("Gene", 1:10)
colnames(gex) <- paste0("Cell", 1:30)

# Create a PSI matrix with 3 rows (sequences) and 30 columns (cells)
PSI <- matrix(runif(90, 0, 1), nrow = 3)  # Random values between 0 and 1
rownames(PSI) <- paste0("Seq", 1:3)
colnames(PSI) <- paste0("Cell", 1:30)

# Add the first additional row
gex <- rbind(gex, c(10, rep(1, 29)))
rownames(gex)[nrow(gex)] <- "Gene11"

PSI <- rbind(PSI, c(1, rep(0, 29)))
rownames(PSI)[nrow(PSI)] <- "Seq4"

# Add the second additional row
gex <- rbind(gex, c(rep(2, 15), rep(1, 15)))
rownames(gex)[nrow(gex)] <- "Gene12"

PSI <- rbind(PSI, rep(0.1, 30))
rownames(PSI)[nrow(PSI)] <- "Seq5"

# Load the pls library
library(pls)
set.seed(1)

# Fit the PLSR model using the fourth row of PSI and transposed gex
model <- plsr(PSI[4, ] ~ t(gex), scale = TRUE, validation = "CV", segments = 10)

# Extract and sort the loading weights for the first component
loadings <- loading.weights(model)
first_component_loadings <- loadings[, 1]
important_variables <- sort(first_component_loadings, decreasing = TRUE)

# Print the important variables for the first component
print(important_variables)

# Plot the important variables for the first component
barplot(important_variables, main = "Important Variables for the First Component", xlab = "Variable", ylab = "Loading Weight")


###### Now try it for the actual matrices #####
gex <- gex_mat_aligned
psi_shortlist <- all_samples_mat_aligned_shortlist
psi_one_seq <- all_samples_mat_aligned["ENSG00000138821.14;SLC39A8;chr4-102285948-102285973-102267871-102268079-102304316-102304481", ]

# Fit the PLSR model using the psi_one_seq and gex.
model <- plsr(psi_one_seq ~ gex, scale = TRUE, validation = "CV")
loadings <- loading.weights(model)
first_component_loadings <- loadings[, 1]
important_variables <- sort((first_component_loadings), decreasing = TRUE)

# Print the important variables for the first component
print(important_variables)

# Plot the important variables for the first component
barplot(important_variables, main = "Important Variables for the First Component", xlab = "Variable", ylab = "Loading Weight")

# Output all plots for this sequence:
outdir <- "~/Dropbox (Harvard University)/02Splicing/test2/"
for (seq in c("ENSG00000138821.14;SLC39A8;chr4-102285948-102285973-102267871-102268079-102304316-102304481")) {
  for (gene in colnames(gex_mat_aligned)) {
    gex_values <- gex_mat_aligned[, gene]
    psi_values <- all_samples_mat_aligned_shortlist[seq, ]
    
    # Get cell line name with the highest PSI value.
    cell_line_max_psi <- names(psi_values)[which.max(psi_values)]
    
    # Get the cell line with the highest gene expression value
    cell_line_max_gex <- names(gex_values)[which.max(gex_values)]
    
    # Check if the cell line with the highest PSI value also has the highest gene expression
      # Create a dataframe for plotting
      plot_df <- data.frame(gex_values, psi_values, cell_line = rownames(gex_mat_aligned))
      
      # Generate the scatter plot
      p <- ggplot(plot_df, aes(x = gex_values, y = psi_values, color = cell_line)) +
        geom_point() +
        theme_minimal() +
        labs(x = "Gene expression", y = paste(seq, "PSI")) +
        ggtitle(paste(important_variables[which(names(important_variables) == gene)], seq, "and", gene))
      
      # Save the plot to a file
      ggsave(filename = paste0(outdir, "oneseq_",gene, "_", seq, "_scatter_plot_shortlist.png"), plot = p, width = 8, height = 6)
  }
}

# Plot heatmap of the psi_shortlist
pheatmap::pheatmap(all_samples_mat_aligned_shortlist, cluster_rows = T, cluster_cols = T, main = "PSI values for all samples")
