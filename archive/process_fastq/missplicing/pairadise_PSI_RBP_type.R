library(tidyverse)
library(vroom)
library(data.table)
library(future)
library(future.apply)

# if pairadise is not installed, install it
if(!requireNamespace("PAIRADISE", quietly = TRUE)) {
  BiocManager::install("PAIRADISE")
}
# Also for BiocParallel
if(!requireNamespace("BiocParallel", quietly = TRUE)) {
  BiocManager::install("BiocParallel")
}
library(PAIRADISE)

reverse_complement <- function(dna_seq) {
  complement <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
  nucleotides <- unlist(strsplit(dna_seq, ""))
  complement_nucleotides <- complement[nucleotides]
  reverse_complement_seq <- paste(rev(complement_nucleotides), collapse = "")
  return(reverse_complement_seq)
}

args <- commandArgs(trailingOnly = TRUE)
genename <- args[1]
output_filename <- args[2]

# fwrite(all_sample_reps, file = "/Volumes/broad_dawnccle/processed_data/missplicing_processed_df/V5/all_sample_reps_PSI.csv")
all_sample_reps <- fread("/broad/dawnccle/processed_data/missplicing_processed_df/V5/all_sample_reps_PSI.csv")
all_sample_reps <- all_sample_reps %>% 
  # Change the DBTR05MG to DBTRG05MG
  mutate(condition = str_replace(condition, "DBTR05MG", "DBTRG05MG")) %>%
  # Change MEWO to MeWo
  mutate(condition = str_replace(condition, "MEWO", "MeWo")) %>%
  # Change JHOM to JHOM1
  mutate(condition = str_replace(condition, "JHOM", "JHOM1")) %>%
  # Remove condition == HEK
  filter(condition != "HEK") %>%
  mutate(celltype = if_else(str_detect(condition, "_"), 
                            str_extract(condition, ".*(?=_)"), 
                            condition) %>% str_to_upper())

print("Finished reading data")

get_paradise_output <- function(df, condition1, condition2) {
  # Filter data frames based on conditions
  condition1_df <- df %>% filter(celltype %in% condition1) # %>% filter(offset == "0:0:0")
  condition2_df <- df %>% filter(celltype %in% condition2) # %>% filter(offset == "0:0:0")
  
  # Create a list of the first 10 unique indices
  unique_indices <- unique(condition1_df$index)
  
  # Function to pad vectors with zeros
  pad_with_zeros <- function(vec, target_len = 3) {
    if (length(vec) < target_len) {
      vec <- c(vec, rep(0, target_len - length(vec)))
    }
    return(vec)
  }
  
  # Setup future plan to use multiple threads
  plan(multisession, workers = 8)
  
  # Create the paradise data frame using dplyr operations and future.apply for parallel processing
  paradise_df <- future_lapply(unique_indices, function(idx) {
    condition1_index_df <- condition1_df %>% filter(index == idx) %>% arrange(sample)
    condition2_index_df <- condition2_df %>% filter(index == idx) %>% arrange(sample)
    
    I1 <- pad_with_zeros(condition1_index_df$count)
    S1 <- pad_with_zeros(condition1_index_df$skipped)
    I2 <- pad_with_zeros(condition2_index_df$count)
    S2 <- pad_with_zeros(condition2_index_df$skipped)
    
    data.frame(
      ExonID = idx, 
      I1 = paste(I1, collapse = ","), 
      S1 = paste(S1, collapse = ","), 
      I2 = paste(I2, collapse = ","), 
      S2 = paste(S2, collapse = ","), 
      I_len = 1, 
      S_len = 1
    )
  }) %>% bind_rows()
  
  return(paradise_df)
}

rbp_clusters <- read_csv("/broad/dawnccle/melange/data/RBP_cluster_matrix.csv") 
# rbp_clusters <- read_csv("~/Dropbox (Harvard University)/02Splicing/SplicingManuscript/RBP_cluster_matrix.csv") 
# rename column 1 to cell line
colnames(rbp_clusters)[1] <- "celltype"

# Extract columns celltype and genename. Genename is the gene name of the RBP, which is a colname
# genename <- "RBFOX1"
gene_shortlist <- rbp_clusters[, c("celltype", genename)]
colnames(gene_shortlist)[2] <- "cluster"

cluster1_condition <- gene_shortlist %>% filter(cluster == 1) %>% mutate(celltype = str_to_upper(celltype)) %>% pull(celltype)
cluster2_condition <- gene_shortlist %>% filter(cluster == 2) %>% mutate(celltype = str_to_upper(celltype)) %>% pull(celltype)

formatted_df <- get_paradise_output(all_sample_reps, cluster1_condition, cluster2_condition)
print("Finished formatting data")
write_tsv(formatted_df, file = output_filename)
