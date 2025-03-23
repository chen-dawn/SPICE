library(tidyverse)
library(vroom)
library(data.table)
library(pheatmap)
library(preprocessCore)
library(purrr)
library(RColorBrewer)
library(ggpubr)
library(viridis)

reverse_complement <- function(dna_seq) {
  complement <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
  nucleotides <- unlist(strsplit(dna_seq, ""))
  complement_nucleotides <- complement[nucleotides]
  reverse_complement_seq <- paste(rev(complement_nucleotides), collapse = "")
  return(reverse_complement_seq)
}

output_filepath <- "C:/Users/dawnxi/Dropbox (Harvard University)/02Splicing/SplicingManuscript/figure_outputs"

# final_psi_table_filtered <- fread("U:/processed_data/reprocess_250221/count_normalized_chimeric_rate_considering_included/WT_all_samples_PSI_count_table.csv")
# final_psi_table_filtered <- final_psi_table_filtered %>% 
#   filter(!(condition %in% c("K562WT", "K562K700E"))) %>% 
#   filter(!(condition %in% c("JHOM1", "RVH421", "KNS60", "OVTOKO"))) %>% 
#   mutate(total_count = included_count + skipped_count) %>%
#   filter(total_count >= 20) %>%
#   mutate(index_offset = paste(index, offset, sep = "__")) %>% 
#   separate(offset, into = c("upstream_offset", "downstream_offset", "const_offset"), sep = ":") %>% 
#   mutate(upstream_offset = as.integer(upstream_offset)) %>% 
#   mutate(downstream_offset = as.integer(downstream_offset)) %>%
#   mutate(const_offset = as.integer(const_offset)) %>% 
#   filter(abs(upstream_offset) != 1 & abs(downstream_offset)!= 1) %>% 
#   mutate(PSI = included_count/(included_count + skipped_count))

# fwrite(final_psi_table_filtered, "U:/processed_data/reprocess_250221/count_normalized_chimeric_rate_considering_included/WT_all_samples_PSI_count_table_filtered.csv")
final_psi_table_filtered <- fread("U:/processed_data/reprocess_250221/count_normalized_chimeric_rate_considering_included/WT_all_samples_PSI_count_table_filtered.csv")
twist_barcodes <- read_csv("U:/melange/data/guide_library/20230130_twist_library_v3.csv") %>%
  mutate(barcodeRevcomp = sapply(barcode, reverse_complement))

final_psi_event_sequences <- final_psi_table_filtered %>%
  left_join(twist_barcodes, by = c("index" = "ID"))

# > colnames(high_in_kelly_sequences)
# [1] "ID"                  "offset"              "upsilon"
# [4] "upsilon_reverse"     "tau"                 "tau_reverse"
# [7] "num_na_per_row"      "row_max"             "row_min"
# [10] "max_sample"          "min_sample"          "target_cell_type"
# [13] "barcode"             "upstreamIntronSeq"   "skippedExonSeq"
# [16] "downstreamIntronSeq" "librarySequence"     "twistSequence"
# [19] "barcodeRevcomp"


# Adjust the sequences upstreamIntronSeq, skippedExonSeq, downstreamIntronSeq based on offset.
final_psi_event_sequences <- final_psi_event_sequences %>%
  mutate(upstreamIntron_len = nchar(upstreamIntronSeq)) %>%
  mutate(downstreamIntron_len = nchar(downstreamIntronSeq)) %>%
  mutate(skippedExon_len = nchar(skippedExonSeq)) %>%
  mutate(upstreamIntron_len_adj = upstreamIntron_len + upstream_offset) %>%
  mutate(downstreamIntron_len_adj = downstreamIntron_len - downstream_offset) %>%
  mutate(skippedExon_len_adj = skippedExon_len - upstream_offset + downstream_offset) %>%
  mutate(upstreamIntronSeq_adj = substr(librarySequence, 1, upstreamIntron_len_adj)) %>%
  mutate(skippedExonSeq_adj = substr(librarySequence, upstreamIntron_len_adj + 1, upstreamIntron_len_adj + skippedExon_len_adj)) %>%
  mutate(downstreamIntronSeq_adj = substr(librarySequence, upstreamIntron_len_adj + skippedExon_len_adj + 1, upstreamIntron_len_adj + skippedExon_len_adj + downstreamIntron_len_adj))

########################################################
###### Look at how many seuqences have exon <20bp ######
########################################################
# Plot exon len distribution
exon_lens <- final_psi_event_sequences %>% 
  select(index_offset, skippedExon_len_adj) %>%
  distinct()

# Create histogram with enhanced aesthetics
p1 <- ggplot(exon_lens, aes(x = skippedExon_len_adj)) +
  geom_histogram(binwidth = 1, fill = "#0073C2", color = "black", alpha = 0.8) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  labs(
    title = "Distribution of Exon Lengths",
    x = "Exon Length (bp)",
    y = "Frequency"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 12)
  )

ggsave(file.path(output_filepath, "micro_exon_length_distribution.pdf"), 
       plot = p1, width = 8, height = 6, dpi = 300)  # Publication quality

exon_lens_perfect <- exon_lens %>%
  filter(grepl("__0:0:0", index_offset))

p2 <- ggplot(exon_lens_perfect, aes(x = skippedExon_len_adj)) +
  geom_histogram(binwidth = 1, fill = "#0073C2", color = "black", alpha = 0.8) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  labs(
    title = "Distribution of Exon Lengths (reference only)",
    x = "Exon Length (bp)",
    y = "Frequency"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 12)
  )

ggsave(file.path(output_filepath, "micro_exon_length_distribution_perfect_only.pdf"), 
       plot = p2, width = 8, height = 6, dpi = 300)  # Publication quality

# How many sequences have exon < 20bp
small_exon_index <- exon_lens %>% 
  filter(skippedExon_len_adj < 20) %>% 
  pull(index_offset)

small_exon_seq <- final_psi_event_sequences %>% 
  filter(index_offset %in% small_exon_index) 

small_exon_PSI <- small_exon_seq %>% 
  group_by(index_offset, condition) %>% 
  summarise(PSI = mean(PSI, na.rm = T)) %>% 
  pivot_wider(names_from = condition, values_from = PSI) %>% 
  ungroup()

small_exon_mat <- as.matrix(small_exon_PSI %>% select(-index_offset))
rownames(small_exon_mat) <- small_exon_PSI$index_offset
# Filter out rows with > 20% NA
small_exon_mat <- small_exon_mat[rowMeans(is.na(small_exon_mat)) < 0.2, ]
pdf(file.path(output_filepath, "micro_exon_PSI_heatmap.pdf"), width = 10, height = 8)
pheatmap(small_exon_mat, 
         show_rownames = F, 
         show_colnames = T,
         fontsize_row = 8,
         fontsize_col = 8,
         color = viridis(100),
         cluster_rows = T,
         cluster_cols = F,
         main = "PSI distribution of micro exons (<20bp)")
dev.off()

##################################
###### Plot Sequence Motifs ######
##################################
library(ggseqlogo)

# Reference offset only.
ref_offset <- final_psi_event_sequences %>% 
  select(index_offset, upstreamIntronSeq_adj, skippedExonSeq_adj, downstreamIntronSeq_adj) %>% 
  distinct()

# Function to safely extract substrings
safe_substr <- function(seq, start, stop) {
  ifelse(nchar(seq) >= stop, substr(seq, start, stop), NA)
}

# Extract Upstream Splice Site Motif (Last 20 bp of Upstream Intron + First 10 bp of Skipped Exon)
ref_offset <- ref_offset %>%
  mutate(
    upstream_motif = paste0(
      safe_substr(upstreamIntronSeq_adj, nchar(upstreamIntronSeq_adj) - 19, nchar(upstreamIntronSeq_adj)),
      safe_substr(skippedExonSeq_adj, 1, 3)
    )
  )

# Extract Downstream Splice Site Motif (Last 10 bp of Skipped Exon + First 20 bp of Downstream Intron)
ref_offset <- ref_offset %>%
  mutate(
    downstream_motif = paste0(
      safe_substr(skippedExonSeq_adj, nchar(skippedExonSeq_adj) - 2, nchar(skippedExonSeq_adj)),
      safe_substr(downstreamIntronSeq_adj, 1, 6)
    )
  )

# Remove any NA values
ref_offset <- ref_offset %>%
  filter(!is.na(upstream_motif) & !is.na(downstream_motif)) %>% 
  mutate(upstream_nchar = nchar(upstream_motif)) %>% 
  mutate(downstream_nchar = nchar(downstream_motif)) %>%
  filter(upstream_nchar == 23 & downstream_nchar == 9)

# Plot Upstream Splice Site Motif Logo
p_upstream <- ggseqlogo(ref_offset$upstream_motif, method = "bits") + ylim(0,2)
p_downstream <- ggseqlogo(ref_offset$downstream_motif, method = "bits") + ylim(0,2)

# Save the plots
p <- gridExtra::grid.arrange(p_upstream, p_downstream, ncol = 2)
ggsave(filename = paste0(output_filepath, "/micro_exon_motif_logo.pdf"), 
       plot = p, 
       width = 12, height = 3, dpi = 300)
