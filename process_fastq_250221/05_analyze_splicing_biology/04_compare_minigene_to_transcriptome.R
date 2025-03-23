library(tidyverse)
library(vroom)
library(data.table)
library(pheatmap)
library(preprocessCore)
library(purrr)
library(RColorBrewer)
library(ggpubr)
library(viridis)
library(ggpubr)

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
unique_celllines <- unique(final_psi_table_filtered$condition)


transcriptome_se <- fread("U:/for_anisha/SE_JC_filtered_with_id.txt")
transcriptome_se <- transcriptome_se %>% 
   filter(cellline_name %in% unique_celllines) 

final_psi_table_filtered <- final_psi_table_filtered %>% 
  filter(grepl("__0:0:0", index_offset)) %>% 
  separate(index, into = c("GeneID", "geneSymbol", "ID"), sep = ";")
  
transcriptome_minigene_merged <- merge(transcriptome_se, final_psi_table_filtered, by.x = c("ID", "cellline_name"), by.y = c("ID", "condition"))
transcriptome_minigene_merged <- transcriptome_minigene_merged %>% 
  group_by(ID, cellline_name) %>%
  summarise(PSI_transcriptome = mean(IncLevel1, na.rm = T),
            PSI_minigene = mean(PSI, na.rm = T))

p1 <- ggplot(transcriptome_minigene_merged, aes(x = PSI_transcriptome, y = PSI_minigene)) +
  geom_point( fill = "lightblue", shape = 21, alpha = 0.6, size = 1) +  # High contrast points
  # geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = 1) +  # Reference line
  labs(title = "Comparison of PSI values between Transcriptome and Minigene",
       x = "PSI (Transcriptome)",
       y = "PSI (Minigene)") +
  coord_fixed() +  # Equal aspect ratio for x and y axes
  theme_classic(base_size = 16) +  # Clean, publication-ready theme
  theme(
    axis.text.x = element_text(size = 14),  # Readable x-axis
    axis.text.y = element_text(size = 14),  # Readable y-axis
    axis.title = element_text(size = 16),   # Readable axis titles
    plot.title = element_text(size = 18, hjust = 0.5),  # Centered title
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )
ggsave(file.path(output_filepath, "transcriptome_minigene_comparison.pdf"), p1, width = 8, height = 8)

library(ggpointdensity)
p2 <- ggplot(transcriptome_minigene_merged, aes(x = PSI_transcriptome, y = PSI_minigene)) +
  geom_pointdensity() +
  scale_color_viridis_c() +
  labs(title = "Point Density of PSI Comparison",
       x = "PSI (Transcriptome)",
       y = "PSI (Minigene)",
       color = "Density") +
  theme_classic(base_size = 16) +  # Clean, publication-ready theme
  theme(
    axis.text.x = element_text(size = 14),  # Readable x-axis
    axis.text.y = element_text(size = 14),  # Readable y-axis
    axis.title = element_text(size = 16),   # Readable axis titles
    plot.title = element_text(size = 18, hjust = 0.5)  # Centered title
  )
ggsave(file.path(output_filepath, "transcriptome_minigene_comparison_density.pdf"), p2, width = 10, height = 8)
  
# Compute correlation per cell line
transcriptome_minigene_merged <- transcriptome_minigene_merged %>%
  group_by(cellline_name) %>%
  mutate(correlation = cor(PSI_transcriptome, PSI_minigene, use = "complete.obs"))

# Define output path
output_filepath <- "C:/Users/dawnxi/Dropbox (Harvard University)/02Splicing/SplicingManuscript/figure_outputs"
pdf_filename <- file.path(output_filepath, "transcriptome_minigene_comparison_facet.pdf")

# Create the plot
p1 <- ggplot(transcriptome_minigene_merged, aes(x = PSI_transcriptome, y = PSI_minigene)) +
  geom_pointdensity() +
  scale_color_viridis_c() +
  facet_wrap(~cellline_name) +  # Facet by cell line
  stat_cor(method = "pearson", label.x = 0.1, label.y = 0.9, size = 5) +  # Add correlation coefficient
  labs(title = "Comparison of PSI Values Between Transcriptome and Minigene",
       x = "PSI (Transcriptome)",
       y = "PSI (Minigene)") +
  theme_classic(base_size = 12) +  # Publication-quality theme
  theme(
    axis.text.x = element_text(size = 12),  # Readable x-axis
    axis.text.y = element_text(size = 12),  # Readable y-axis
    axis.title = element_text(size = 14),   # Readable axis titles
    plot.title = element_text(size = 16, hjust = 0.5),  # Centered title
    strip.text = element_text(size = 12, face = "bold"),  # Facet labels
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )

# Save as high-resolution PDF
ggsave(pdf_filename, plot = p1, width = 20, height = 16, dpi = 300)

# twist_barcodes <- read_csv("U:/melange/data/guide_library/20230130_twist_library_v3.csv") %>%
#   mutate(barcodeRevcomp = sapply(barcode, reverse_complement))
# 
# final_psi_event_sequences <- final_psi_table_filtered %>%
#   left_join(twist_barcodes, by = c("index" = "ID"))

#### Look at the things that are in the bottom right corner. Do they have native short introns? 
bottom_right_quad <- transcriptome_minigene_merged %>%
  ungroup() %>% 
  filter(PSI_transcriptome > 0.5 & PSI_minigene < 0.5) %>% 
  select(ID) %>% 
  distinct() %>% 
  pull(ID) 

top_right_quad <- transcriptome_minigene_merged %>%
  ungroup() %>% 
  filter(PSI_transcriptome > 0.5 & PSI_minigene > 0.5) %>% 
  select(ID) %>% 
  distinct() %>% 
  pull(ID)

top_left_quad <- transcriptome_minigene_merged %>%
  ungroup() %>% 
  filter(PSI_transcriptome < 0.5 & PSI_minigene > 0.5) %>% 
  select(ID) %>% 
  distinct() %>% 
  pull(ID)

bottom_left_quad <- transcriptome_minigene_merged %>%
  ungroup() %>% 
  filter(PSI_transcriptome < 0.5 & PSI_minigene < 0.5) %>%
  select(ID) %>%
  distinct() %>% 
  pull(ID)

transcriptome <-  merge(transcriptome_se, final_psi_table_filtered, by.x = c("ID", "cellline_name"), by.y = c("ID", "condition")) %>% 
  select(-cellline_name, -IncLevel1, - IJC_SAMPLE_1, -SJC_SAMPLE_1, -sample, -included_count, -skipped_count, -total_count, - PSI) %>%
  distinct() %>% 
  mutate(upstreamESReal = ifelse(strand == "-", downstreamES, upstreamES)) %>%
  mutate(upstreamEEReal = ifelse(strand == "-", downstreamEE, upstreamEE)) %>%
  mutate(downstreamESReal = ifelse(strand == "-", upstreamES, downstreamES)) %>%
  mutate(downstreamEEReal = ifelse(strand == "-", upstreamEE, downstreamEE)) %>%
  mutate(upstreamISReal = ifelse(strand== "-", exonEnd, upstreamEE)) %>%
  mutate(upstreamIEReal = ifelse(strand== "-", downstreamES, exonStart_0base)) %>%
  mutate(downstreamISReal = ifelse(strand =="-", upstreamEE, exonEnd)) %>%
  mutate(downstreamIEReal = ifelse(strand == "-", exonStart_0base, downstreamES)) %>% 
  # Annotate with quadrant
  mutate(quadrant = case_when(
    ID %in% bottom_right_quad ~ "bottom_right",
    ID %in% top_right_quad ~ "0top_right",
    ID %in% top_left_quad ~ "0top_left",
    ID %in% bottom_left_quad ~ "bottom_left",
    TRUE ~ "NONE"
  )) %>% 
  filter(quadrant != "NONE")

transcriptome <- transcriptome %>%
  mutate(upstream_intron_len = upstreamIEReal - upstreamISReal,
         downstream_intron_len = downstreamIEReal - downstreamISReal,
         exon_len = exonEnd  - exonStart_0base,
         upstream_plus_downstream_intron_len = upstream_intron_len + downstream_intron_len) 

save_density_plot <- function(data, feature, title, filename) {
  p <- ggplot(data, aes_string(feature)) + 
    geom_density(fill = "gray70", color = "black", alpha = 0.7) +  # High contrast density plot
    facet_wrap(~quadrant, scales = "free") +  # Facet by quadrant
    scale_x_log10() +  # Log scale for better visualization
    labs(title = title, x = feature, y = "Density") +
    theme_classic(base_size = 16) +  # Publication-ready theme
    theme(
      axis.text.x = element_text(size = 12),  
      axis.text.y = element_text(size = 12),  
      axis.title = element_text(size = 14),   
      plot.title = element_text(size = 16, hjust = 0.5),  # Centered title
      strip.text = element_text(size = 14, face = "bold")  # Facet labels
    )
  
  # Save the figure as a high-resolution PDF
  ggsave(file.path(output_filepath, filename), plot = p, width = 8, height = 6, dpi = 300)
}

# Generate and save plots
save_density_plot(transcriptome, "upstream_intron_len", "Density of Upstream Intron Length", "upstream_intron_length.pdf")
save_density_plot(transcriptome, "downstream_intron_len", "Density of Downstream Intron Length", "downstream_intron_length.pdf")
save_density_plot(transcriptome, "exon_len", "Density of Exon Length", "exon_length.pdf")
save_density_plot(transcriptome, "upstream_plus_downstream_intron_len", "Density of Combined Intron Length", "combined_intron_length.pdf")

p3 <- ggplot(transcriptome, aes(upstream_intron_len, downstream_intron_len)) +
  geom_pointdensity() + 
  scale_color_viridis_c() +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic(base_size = 16) + 
  labs(x = "Upstream Intron Length", y = "Downstream Intron Length", color = "Quadrant") +
  facet_wrap(~quadrant)
ggsave(file.path(output_filepath, "transcriptome_intron_length_scatter.pdf"), p3, width = 12, height = 10, dpi = 300)
 