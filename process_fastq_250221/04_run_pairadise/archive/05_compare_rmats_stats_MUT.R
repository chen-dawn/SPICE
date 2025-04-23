library(tidyverse)
library(pheatmap)
library(data.table)

reverse_complement <- function(dna_seq) {
  complement <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
  nucleotides <- unlist(strsplit(dna_seq, ""))
  complement_nucleotides <- complement[nucleotides]
  reverse_complement_seq <- paste(rev(complement_nucleotides), collapse = "")
  return(reverse_complement_seq)
}

twist_barcodes <- read_csv("U:/melange/data/guide_library/20230130_twist_library_v3.csv") %>%
  mutate(barcodeRevcomp = sapply(barcode, reverse_complement))

######## Look at rmats stats ########
combined_psi <- read_tsv("U:/processed_data/reprocess_250221/pairadise_indiv_PSI/MUT_PSI_combined_output_indiv.tsv")
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

wanted_pairs <- c("MUT_MUT-plx317_U2AF1_Q157A_MUT-plx317_U2AF1_WT",
                  "MUT_MUT-sgCh3-1_MUT-sgRBM10",
                  "MUT_MUT-sgCh3-1_MUT-sgRBM5",
                  "MUT_MUT-sgCh3-1_MUT-sgRUBP1",
                  "MUT_splicelib_U2AF1_S34F_splicelib_U2AF1_WT",
                  "MUT_splicelib_sgCh3_splicelib_ZRSR2")

combined_psi_filtered <- combined_psi %>% 
  filter(Folder %in% wanted_pairs) %>%
  filter(count_sum_average1 >30) %>% 
  filter(count_sum_average2 > 30) %>% 
  mutate(log2_PSI_ratio = log2(PSI_ratio), log2_PSI_reverse_ratio = log2(PSI_reverse_ratio)) %>% 
  separate(ExonID, sep = "__", into =c("index", "offset"), remove = FALSE) %>%
  separate(offset, into = c("skipped_exon_start", "skipped_exon_end", "downstream_exon_start"), sep = ":", remove = FALSE) %>%
  filter(abs(as.integer(skipped_exon_start)) != 1 & abs(as.integer(skipped_exon_end)) != 1) 

# Plot the volcano plot
ggplot(combined_psi_filtered, aes(log2(PSI_ratio), -log10(FDR))) + 
  geom_point() + 
  facet_wrap(~Folder)
  
# Look at S34F pair only.
combined_psi_filtered_S34F <- combined_psi_filtered %>% 
  filter(Folder == "MUT_splicelib_U2AF1_S34F_splicelib_U2AF1_WT")


######## Process and Save the Individual Sequences as Fasta ######
# Define a function to process and save sequences
process_and_save <- function(df, name) {
  # Define file paths
  upstream_fasta <- paste0("U:/processed_data/reprocess_250221/pairadise_indiv_PSI/seq_output/", name, "_upstreamIntronSeq_adj.fasta")
  skipped_fasta <- paste0("U:/processed_data/reprocess_250221/pairadise_indiv_PSI/seq_output/", name, "_skippedExonSeq_adj.fasta")
  downstream_fasta <- paste0("U:/processed_data/reprocess_250221/pairadise_indiv_PSI/seq_output/", name, "_downstreamIntronSeq_adj.fasta")
  
  # Write upstream intron sequences
  fileConn <- file(upstream_fasta, "w")
  apply(df, 1, function(row) {
    cat(paste0(">", paste0(row["barcode"], "_", row["offset"]), "\n", row["upstreamIntronSeq_adj"], "\n"), file = fileConn)
  })
  close(fileConn)
  
  # Write skipped exon sequences
  fileConn <- file(skipped_fasta, "w")
  apply(df, 1, function(row) {
    cat(paste0(">", paste0(row["barcode"], "_", row["offset"]), "\n", row["skippedExonSeq_adj"], "\n"), file = fileConn)
  })
  close(fileConn)
  
  # Write downstream intron sequences
  fileConn <- file(downstream_fasta, "w")
  apply(df, 1, function(row) {
    cat(paste0(">", paste0(row["barcode"], "_", row["offset"]), "\n", row["downstreamIntronSeq_adj"], "\n"), file = fileConn)
  })
  close(fileConn)
  
  # Save CSV file
  csv_path <- paste0("U:/processed_data/reprocess_250221/pairadise_indiv_PSI/seq_output/MUT_", name, "_seq.csv")
  write_csv(df, csv_path)
  
  # Message to indicate completion
  cat("FASTA and CSV files saved for:", name, "\n")
}

# Process and save each dataset
datasets <- list(
  "RBM10" = combined_psi_filtered %>% filter(Folder == "MUT_MUT-sgCh3-1_MUT-sgRBM10"),
  "RBM5" = combined_psi_filtered %>% filter(Folder == "MUT_MUT-sgCh3-1_MUT-sgRBM5"),
  "FUBP1" = combined_psi_filtered %>% filter(Folder == "MUT_MUT-sgCh3-1_MUT-sgRUBP1"),
  "U2AF1_S34F" = combined_psi_filtered %>% filter(Folder == "MUT_splicelib_U2AF1_S34F_splicelib_U2AF1_WT")
)

for (name in names(datasets)) {
  df <- datasets[[name]] %>%
    left_join(twist_barcodes, by = c("index" = "ID")) %>%
    mutate(
      upstream_offset = as.integer(skipped_exon_start),
      downstream_offset = as.integer(skipped_exon_end),
      const_offset = as.integer(downstream_exon_start),
      upstreamIntron_len = nchar(upstreamIntronSeq),
      downstreamIntron_len = nchar(downstreamIntronSeq),
      skippedExon_len = nchar(skippedExonSeq),
      upstreamIntron_len_adj = upstreamIntron_len + upstream_offset,
      downstreamIntron_len_adj = downstreamIntron_len - downstream_offset,
      skippedExon_len_adj = skippedExon_len - upstream_offset + downstream_offset,
      upstreamIntronSeq_adj = substr(librarySequence, 1, upstreamIntron_len_adj),
      skippedExonSeq_adj = substr(librarySequence, upstreamIntron_len_adj + 1, upstreamIntron_len_adj + skippedExon_len_adj),
      downstreamIntronSeq_adj = substr(librarySequence, upstreamIntron_len_adj + skippedExon_len_adj + 1, upstreamIntron_len_adj + skippedExon_len_adj + downstreamIntron_len_adj)
    )
  
  process_and_save(df, name)
}

cat("Processing complete for all datasets.\n")

# Make the "background sequence" file 
upstream_fasta <- paste0("U:/processed_data/reprocess_250221/pairadise_indiv_PSI/seq_output/", "background_upstreamIntronSeq_adj.fasta")
skipped_fasta <- paste0("U:/processed_data/reprocess_250221/pairadise_indiv_PSI/seq_output/", "background_skippedExonSeq_adj.fasta")
downstream_fasta <- paste0("U:/processed_data/reprocess_250221/pairadise_indiv_PSI/seq_output/", "background_downstreamIntronSeq_adj.fasta")

# Write upstream intron sequences
fileConn <- file(upstream_fasta, "w")
apply(twist_barcodes, 1, function(row) {
  cat(paste0(">", paste0(row["barcode"], "_", row["offset"]), "\n", row["upstreamIntronSeq"], "\n"), file = fileConn)
})
close(fileConn)

# Write skipped exon sequences
fileConn <- file(skipped_fasta, "w")
apply(twist_barcodes, 1, function(row) {
  cat(paste0(">", paste0(row["barcode"], "_", row["offset"]), "\n", row["skippedExonSeq"], "\n"), file = fileConn)
})
close(fileConn)

# Write downstream intron sequences
fileConn <- file(downstream_fasta, "w")
apply(twist_barcodes, 1, function(row) {
  cat(paste0(">", paste0(row["barcode"], "_", row["offset"]), "\n", row["downstreamIntronSeq"], "\n"), file = fileConn)
})
close(fileConn)