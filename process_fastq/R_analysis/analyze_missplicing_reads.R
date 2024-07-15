library(tidyverse)
library(data.table)
library(pheatmap)

reverse_complement <- function(dna_seq) {
  complement <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
  nucleotides <- unlist(strsplit(dna_seq, ""))
  complement_nucleotides <- complement[nucleotides]
  reverse_complement_seq <- paste(rev(complement_nucleotides), collapse = "")
  return(reverse_complement_seq)
}

out_dir <- "~/Dropbox (Harvard University)/02Splicing/library_47k_missplicing/data_analysis/"
dir.create(out_dir, showWarnings = FALSE)
all_files_df <- fread("~/Dropbox (Harvard University)/02Splicing/library_47k_missplicing/V5_results/umi_count_merged_to_ref_normalized.csv")
barcodes <- read_csv("~/melange/data/guide_library/20230130_twist_library_v3.csv") %>%
  mutate(barcodeRevcomp = sapply(barcode, reverse_complement))
all_sample_reps <- fread("/Volumes/broad_dawnccle/processed_data/missplicing_processed_df/V5/all_sample_reps_PSI.csv")


all_samples_wide <- all_sample_reps %>% 
  mutate(PSI = count/(count + skipped)) %>% 
  filter((count + skipped) > 30) %>%
  select(-count, -skipped, - mode, -offset, -condition) %>%
  pivot_wider(names_from = sample, values_from = PSI, values_fill = -1) 

all_samples_mat <- as.matrix(all_samples_wide[, -1])
rownames(all_samples_mat) <- all_samples_wide$index
# pheatmap(all_samples_mat, cluster_rows = T, cluster_cols = F, show_rownames = F,
#          show_colnames = T,
#          color = colorRampPalette(c("black", "black", "#45abd7", "white", "#cf1c1a"))(100),
#          main = "PSI values for all samples")

calculate_tau <- function(row) {
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

# Apply the function to each row of the matrix
all_samples_wide <- all_sample_reps %>%
  filter(!condition %in% c("KMRC1", "A172", "OSRC2", "OVTOKO", "Kelly_old", "SKNAS_Nuc", "K562_1ugNuc")) %>% 
  filter((count + skipped) > 30) %>%
  mutate(PSI = count/(count + skipped)) %>% 
  select(-count, -skipped, - mode, -offset, -condition) %>%
  pivot_wider(names_from = sample, values_from = PSI, values_fill = NA) 


all_samples_mat <- as.matrix(all_samples_wide[, -1])
rownames(all_samples_mat) <- all_samples_wide$index

all_samples_condition_wide <- all_sample_reps %>%
  filter(!condition %in% c("KMRC1", "A172", "OSRC2", "OVTOKO", "Kelly_old", "SKNAS_Nuc", "K562_1ugNuc")) %>% 
  mutate(PSI = count/(count + skipped)) %>% 
  group_by(index, condition) %>%
  summarise(PSI = mean(PSI)) %>%
  ungroup() %>%
  pivot_wider(names_from = condition, values_from = PSI, values_fill = NA) 
all_samples_condition_mat <- as.matrix(all_samples_condition_wide[, -1])
rownames(all_samples_condition_mat) <- all_samples_condition_wide$index
tau_score <- apply(1- all_samples_condition_mat, 1, calculate_tau)

ggplot(data.frame(tau = tau_score), aes(tau)) +
  geom_histogram(bins = 50) +
  labs(title = "Distribution of tau scores", x = "Tau score", y = "Frequency")

# Get sequences with tau > 0.95. We want the row names
high_tau_rows <- rownames(all_samples_condition_mat)[which(tau_score > 0.6)]
high_tau_indices <- which(rownames(all_samples_mat) %in% high_tau_rows)
high_tau_sequences <-all_samples_mat[high_tau_indices, ]
pheatmap(high_tau_sequences, cluster_rows = F, cluster_cols = F, show_rownames = F,
         show_colnames = T,
         color = colorRampPalette(c("#45abd7", "white", "#cf1c1a"))(100),
         main = "PSI values for high tau sequences", na_col = "black")

# Plot the histogram of all PSI values.
ggplot(all_sample_reps %>% mutate(PSI = count/(count + skipped)), aes(PSI)) +
  geom_histogram(bins = 100) +
  labs(title = "Distribution of PSI values", x = "PSI", y = "Frequency")


# Get the sequences with PSI == 1
PSI_1 <- all_sample_reps %>%
  filter(!condition %in% c("KMRC1", "A172", "OSRC2", "OVTOKO", "Kelly_old", "SKNAS_Nuc", "K562_1ugNuc")) %>% 
  mutate(PSI = count/(count + skipped)) %>% 
  filter(PSI == 1) 


## Plot with the old:
good_umis <- read_csv("~/Dropbox (Harvard University)/02Splicing/library_47k_novaseq/240217_rerun_umi_dedup_all.csv")
colnames(good_umis) <- c("cb", "unspliced", "included", "skipped", "sample")


# ggplot(good_umis %>% filter(condition == "KELLY") %>% filter((included + skipped) > 30) %>% 
#          mutate(PSI = included/(included + skipped)), aes(PSI)) + geom_histogram(bins = 100) + theme_classic()
# Calculate PSI and convert data to a matrix format. 
good_umis_reps <- good_umis %>% 
  mutate(PSI = included/(included + skipped)) %>% 
  filter((included + skipped) > 30) %>% 
  select(cb, PSI, sample)

new_PSI_reps <- all_sample_reps %>% 
  mutate(PSI_new = count/(count + skipped)) %>% 
  filter((count + skipped) > 30) %>%
  select(index, PSI_new, sample) %>% 
  # merge with barcodes to get the cell line
  left_join(barcodes, by = c("index" = "ID")) %>% 
  select(barcodeRevcomp, PSI_new, sample)

# Merge the good umis with the new umis
all_psi_reps <- merge(good_umis_reps, new_PSI_reps, by.x = c("sample", "cb"), 
                      by.y = c("sample", "barcodeRevcomp"))

ggplot(all_psi_reps[1:50000, ], aes(PSI, PSI_new)) + geom_pointdensity() + theme_classic() + 
  labs(title = "Comparison of PSI values between old and new data", x = "Old PSI", y = "New PSI")

# Plot sample to sample replicate correlation.
sample_correlation <- new_PSI_reps %>% 
  pivot_wider(names_from = sample, values_from = PSI_new, values_fill = NA) %>% 
  select(-barcodeRevcomp) %>% 
  as.matrix() %>% 
  cor(use = "pairwise.complete.obs")

# Plot on heatmap
pheatmap(sample_correlation, cluster_rows = T, cluster_cols = T, show_rownames = T,
         show_colnames = T,
         main = "Sample to sample correlation")

# Look at andrea sequences. 
andrea_sequences <- read_csv("~/Downloads/Andrea_Sequences.csv")
# Merge with barcodes.
andrea_sequences <- andrea_sequences %>% 
  left_join(barcodes, by = "ID") %>% 
  select(barcodeRevcomp, ID, LibNum) 

# Filter from all_psi_reps
andrea_sequences_subset <- all_psi_reps %>% 
  filter(cb %in% andrea_sequences$barcodeRevcomp) %>% 
  left_join(andrea_sequences, by = c("cb" = "barcodeRevcomp")) 

# Plot heatmap for PSI
andrea_sequences_subset_wide <- andrea_sequences_subset %>% 
  select(sample, PSI, ID) %>% 
  filter(grepl("HCC38|IPC298|HEK|T47D|Kelly", sample))%>% 
  pivot_wider(names_from = sample, values_from = PSI, values_fill = -1) 

andrea_mat <- as.matrix(andrea_sequences_subset_wide[, -1])
rownames(andrea_mat) <- andrea_sequences_subset_wide$ID

pheatmap(andrea_mat, cluster_rows = F, cluster_cols = F, show_rownames = T,
         color = colorRampPalette(c("black", "black", "#45abd7", "white", "#cf1c1a"))(100),
         main = "PSI (old) values for Andrea sequences", display_numbers = TRUE)

# Plot heatmap for new PSI values
andrea_sequences_subset_wide <- andrea_sequences_subset %>% 
  select(sample, PSI_new, ID) %>% 
  filter(grepl("HCC38|IPC298|HEK|T47D|Kelly", sample))%>% 
  pivot_wider(names_from = sample, values_from = PSI_new, values_fill = -1)

andrea_mat <- as.matrix(andrea_sequences_subset_wide[, -1])
rownames(andrea_mat) <- andrea_sequences_subset_wide$ID

pheatmap(andrea_mat, cluster_rows = F, cluster_cols = F, show_rownames = T,
         color = colorRampPalette(c("black", "black", "#45abd7", "white", "#cf1c1a"))(100),
         main = "PSI (new) values for Andrea sequences", display_numbers = TRUE)

write.csv(andrea_sequences_subset, file = "~/Downloads/Andrea_Sequences_PSI.csv")

# Sequences: 
seq_shortlist <- c("ENSG00000156860.16;FBRS;chr16-30661179-30661215-30660262-30660442-30661303-30661333", 
                   "ENSG00000138821.14;SLC39A8;chr4-102285948-102285973-102267871-102268079-102304316-102304481")

# Get the barcoderevcomp for the sequences
seq_shortlist_barcodes <- barcodes %>% 
  filter(ID %in% seq_shortlist) %>% 
  pull(barcodeRevcomp)

# Filter from old PSI
seq_shortlist_subset <- all_psi_reps %>% 
  filter(cb %in% seq_shortlist_barcodes) 

# Plot heatmap for PSI
seq_shortlist_subset_wide <- seq_shortlist_subset %>% 
  select(sample, PSI, cb) %>% 
  # filter(grepl("HCC38|IPC298|HEK|T47D|Kelly", sample))%>% 
  pivot_wider(names_from = sample, values_from = PSI, values_fill = -1)

seq_shortlist_mat <- as.matrix(seq_shortlist_subset_wide[, -1])
rownames(seq_shortlist_mat) <- seq_shortlist_subset_wide$cb

pheatmap(seq_shortlist_mat, cluster_rows = F, cluster_cols = F, show_rownames = T,
         color = colorRampPalette(c("black", "black", "#45abd7", "white", "#cf1c1a"))(100),
         main = "PSI (old) values for shortlisted sequences", display_numbers = TRUE)


# Plot heatmap for new PSI values
seq_shortlist_subset_wide <- seq_shortlist_subset %>% 
  select(sample, PSI_new, cb) %>% 
  # filter(grepl("HCC38|IPC298|HEK|T47D|Kelly", sample))%>% 
  pivot_wider(names_from = sample, values_from = PSI_new, values_fill = -1)

seq_shortlist_mat <- as.matrix(seq_shortlist_subset_wide[, -1])
rownames(seq_shortlist_mat) <- seq_shortlist_subset_wide$cb

pheatmap(seq_shortlist_mat, cluster_rows = F, cluster_cols = F, show_rownames = T,
         color = colorRampPalette(c("black", "black", "#45abd7", "white", "#cf1c1a"))(100),
         main = "PSI (new) values for shortlisted sequences", display_numbers = TRUE)


high_tau_cb <- c("GGTCATTCAGTATT", "AACGAAGTAGAACT", # High in kelly
                 "CATGCACATACAAG", "TGCTAGGTAGCCTT", # High in HCC38 and T47D
                 "TTTGGGCATTTCCT", # HCC38 AND T47D AND SNU423 (LIVER?)
                 "GGCGAAAGGCCTTT", "TCGTTTTCGGGATG", "ACCAACGTTACCGT", "GCCTAAGTATGTTT", # T47D ONLY
                 "ATTGCGTCATGTTA") #T47D AND GAMG ALSO?

high_tau_cb_subset <- all_psi_reps %>%
  filter(cb %in% high_tau_cb) %>%
  left_join(andrea_sequences, by = c("cb" = "barcodeRevcomp")) 

# Plot heatmap. 
high_tau_cb_subset_wide <- high_tau_cb_subset %>% 
  select(sample, PSI_new, ID) %>% 
  filter(grepl("HCC38|IPC298|HEK|T47D|Kelly", sample))%>% 
  pivot_wider(names_from = sample, values_from = PSI_new, values_fill = -1)

high_tau_cb_mat <- as.matrix(high_tau_cb_subset_wide[, -1])
rownames(high_tau_cb_mat) <- high_tau_cb_subset_wide$ID

pheatmap(high_tau_cb_mat, cluster_rows = F, cluster_cols = F, show_rownames = T,
         color = colorRampPalette(c("black", "black", "#45abd7", "white", "#cf1c1a"))(100),
         main = "PSI (new) values for high tau sequences", display_numbers = TRUE)
