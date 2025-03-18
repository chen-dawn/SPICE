library(tidyverse)
library(zoo)
library(furrr)

plan(multicore, workers = availableCores() - 1)

reverse_complement <- function(dna_seq) {
  complement <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
  paste(rev(complement[unlist(strsplit(dna_seq, ""))]), collapse = "")
}

get_kmer_frequency <- function(seq, k){
  chars <- strsplit(seq, "")[[1]]
  if (length(chars) < k) return(data.frame(kmer = character(), frequency = integer()))
  kmers <- rollapply(chars, k, paste, collapse = "", by = 1)
  kmer_freq_table <- table(kmers) %>% 
    as.data.frame(stringsAsFactors = FALSE) %>% 
    setNames(c(paste0("kmer", k), "frequency"))
  return(kmer_freq_table)
}

get_background_distribution_for_kmer <- function(twist_barcodes, k, sample_n, mode, sample_times) {
  target_seqs <- na.omit(twist_barcodes[[mode]])
  if (length(target_seqs) == 0) return(tibble(kmer = character(), frequency = integer(), sample = integer()))
  
  sample_indices <- map(1:sample_times, ~ sample(target_seqs, sample_n, replace = FALSE))
  
  kmer_freq_list <- future_map2(sample_indices, seq_len(sample_times), ~{
    sampled_seqs <- .x
    
    kmer_counts <- map_dfr(sampled_seqs, ~get_kmer_frequency(.x, k)) %>%
      rename_with(~ paste0("kmer", k), starts_with("kmer")) %>%  # Ensure column name matches expected format
      group_by(!!sym(paste0("kmer", k))) %>% 
      summarise(frequency = sum(frequency), .groups = "drop") %>% 
      mutate(sample = .y)  # Assign correct sample number
    
    kmer_counts
  }, .progress = TRUE)
  
  bind_rows(kmer_freq_list)
}

########## Read in and process ##########
twist_barcodes <- read_csv("U:/melange/data/guide_library/20230130_twist_library_v3.csv") %>%
  mutate(barcodeRevcomp = map_chr(barcode, reverse_complement))  # Vectorized

MUT_FUBP1_seq <- read_csv("U:/processed_data/reprocess_250221/pairadise_indiv_PSI/seq_output/MUT_FUBP1_seq.csv")
MUT_RBM10_seq <- read_csv("U:/processed_data/reprocess_250221/pairadise_indiv_PSI/seq_output/MUT_RBM10_seq.csv")
MUT_RBM5_seq <- read_csv("U:/processed_data/reprocess_250221/pairadise_indiv_PSI/seq_output/MUT_RBM5_seq.csv")

###### Bootstrap background for FUBP1 ########
background_1mer_freq_upstreamIntron <- get_background_distribution_for_kmer(twist_barcodes, 1, nrow(MUT_FUBP1_seq), "upstreamIntronSeq", 1000)
background_1mer_freq_skippedExon <- get_background_distribution_for_kmer(twist_barcodes, 1, nrow(MUT_FUBP1_seq), "skippedExonSeq", 1000)
background_1mer_freq_downstreamIntron <- get_background_distribution_for_kmer(twist_barcodes, 1, nrow(MUT_FUBP1_seq), "downstreamIntronSeq", 1000)


# time_taken <- system.time({
#   background_1mer_freq_upstreamIntron <- get_background_distribution_for_kmer(
#     twist_barcodes, 1, nrow(MUT_FUBP1_seq), "upstreamIntronSeq", 10
#   )
# })
# 
# print(time_taken)
# 
# time_taken <- system.time({
#   background_1mer_freq_upstreamIntron <- get_background_distribution_for_kmer(
#     twist_barcodes, 1, nrow(MUT_FUBP1_seq), "upstreamIntronSeq", 100
#   )
# })
# 
# print(time_taken)
# 
# time_taken <- system.time({
#   background_1mer_freq_upstreamIntron <- get_background_distribution_for_kmer(
#     twist_barcodes, 1, nrow(MUT_FUBP1_seq), "upstreamIntronSeq", 1000
#   )
# })
# 
# print(time_taken)



# Combine background frequency data
background_1mer_freq_table <- bind_rows(
  background_1mer_freq_upstreamIntron %>% mutate(region = "upstreamIntron"),
  background_1mer_freq_skippedExon %>% mutate(region = "skippedExon"),
  background_1mer_freq_downstreamIntron %>% mutate(region = "downstreamIntron")
) %>%
  group_by(region, sample) %>%
  mutate(freq_perc = frequency / sum(frequency)) %>%
  ungroup() %>%
  select(kmer1, region, sample, freq_perc)

# Get the values for FUBP1
MUT_FUBP1_seq <- MUT_FUBP1_seq %>% 
  mutate(
    upstreamIntron_kmer1_freq = map(upstreamIntronSeq_adj, ~get_kmer_frequency(.x, 1)),
    skippedExon_kmer1_freq = map(skippedExonSeq_adj, ~get_kmer_frequency(.x, 1)),
    downstreamIntron_kmer1_freq = map(downstreamIntronSeq_adj, ~get_kmer_frequency(.x, 1))
  )

# Aggregate frequency data for each region
MUT_FUBP1_kmer1_freq <- bind_rows(
  MUT_FUBP1_seq %>% select(ExonID, upstreamIntron_kmer1_freq) %>% unnest(upstreamIntron_kmer1_freq) %>% mutate(region = "upstreamIntron"),
  MUT_FUBP1_seq %>% select(ExonID, skippedExon_kmer1_freq) %>% unnest(skippedExon_kmer1_freq) %>% mutate(region = "skippedExon"),
  MUT_FUBP1_seq %>% select(ExonID, downstreamIntron_kmer1_freq) %>% unnest(downstreamIntron_kmer1_freq) %>% mutate(region = "downstreamIntron")
) %>%
  group_by(ExonID, region) %>%
  mutate(freq_perc = frequency / sum(frequency)) %>%
  ungroup() %>%
  group_by(kmer1, region) %>%
  summarise(avg_freq = mean(freq_perc), .groups = "drop")

# Plot the distribution of 1-mer frequencies for each region
ggplot(background_1mer_freq_table, aes(x = freq_perc)) +
  geom_histogram() +
  geom_vline(data = MUT_FUBP1_kmer1_freq, aes(xintercept = avg_freq), color = "red", linetype = "dashed") +
  facet_wrap(~kmer1 + region, nrow = 4) + 
  labs(title = "1-mer frequency distribution for FUBP1", x = "1-mer", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Bootstrap backgroun for RBM5. 
background_1mer_freq_upstreamIntron <- get_background_distribution_for_kmer(twist_barcodes, 1, sample_n = nrow(MUT_RBM5_seq), mode = "upstreamIntronSeq", sample_times = 1000)
background_1mer_freq_skippedExon <- get_background_distribution_for_kmer(twist_barcodes, 1, sample_n = nrow(MUT_RBM5_seq), mode = "skippedExonSeq", sample_times = 1000)
background_1mer_freq_downstreamIntron <- get_background_distribution_for_kmer(twist_barcodes, 1, sample_n = nrow(MUT_RBM5_seq), mode = "downstreamIntronSeq", sample_times = 1000)
# Combine background frequency data
background_1mer_freq_table <- bind_rows(
  background_1mer_freq_upstreamIntron %>% mutate(region = "upstreamIntron"),
  background_1mer_freq_skippedExon %>% mutate(region = "skippedExon"),
  background_1mer_freq_downstreamIntron %>% mutate(region = "downstreamIntron")
) %>%
  group_by(region, sample) %>%
  mutate(freq_perc = frequency / sum(frequency)) %>%
  ungroup() %>%
  select(kmer1, region, sample, freq_perc)

# Get the values for RBM5
MUT_RBM5_seq <- MUT_RBM5_seq %>% 
  mutate(
    upstreamIntron_kmer1_freq = map(upstreamIntronSeq_adj, ~get_kmer_frequency(.x, 1)),
    skippedExon_kmer1_freq = map(skippedExonSeq_adj, ~get_kmer_frequency(.x, 1)),
    downstreamIntron_kmer1_freq = map(downstreamIntronSeq_adj, ~get_kmer_frequency(.x, 1))
  )

# Aggregate frequency data for each region
MUT_RBM5_kmer1_freq <- bind_rows(
  MUT_RBM5_seq %>% select(ExonID, upstreamIntron_kmer1_freq) %>% unnest(upstreamIntron_kmer1_freq) %>% mutate(region = "upstreamIntron"),
  MUT_RBM5_seq %>% select(ExonID, skippedExon_kmer1_freq) %>% unnest(skippedExon_kmer1_freq) %>% mutate(region = "skippedExon"),
  MUT_RBM5_seq %>% select(ExonID, downstreamIntron_kmer1_freq) %>% unnest(downstreamIntron_kmer1_freq) %>% mutate(region = "downstreamIntron")
) %>%
  group_by(ExonID, region) %>%
  mutate(freq_perc = frequency / sum(frequency)) %>%
  ungroup() %>%
  group_by(kmer1, region) %>%
  summarise(avg_freq = mean(freq_perc), .groups = "drop")

# Plot the distribution of 1-mer frequencies for each region
ggplot(background_1mer_freq_table, aes(x = freq_perc)) +
  geom_histogram() +
  geom_vline(data = MUT_RBM5_kmer1_freq, aes(xintercept = avg_freq), color = "red", linetype = "dashed") +
  facet_wrap(~kmer1 + region, nrow = 4) + 
  labs(title = "1-mer frequency distribution for RBM5", x = "1-mer", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Bootstrap background for RBM10.
background_1mer_freq_upstreamIntron <- get_background_distribution_for_kmer(twist_barcodes, 1, sample_n = nrow(MUT_RBM10_seq), mode = "upstreamIntronSeq", sample_times = 1000)
background_1mer_freq_skippedExon <- get_background_distribution_for_kmer(twist_barcodes, 1, sample_n = nrow(MUT_RBM10_seq), mode = "skippedExonSeq", sample_times = 1000)
background_1mer_freq_downstreamIntron <- get_background_distribution_for_kmer(twist_barcodes, 1, sample_n = nrow(MUT_RBM10_seq), mode = "downstreamIntronSeq", sample_times = 1000)
# Combine background frequency data
background_1mer_freq_table <- bind_rows(
  background_1mer_freq_upstreamIntron %>% mutate(region = "upstreamIntron"),
  background_1mer_freq_skippedExon %>% mutate(region = "skippedExon"),
  background_1mer_freq_downstreamIntron %>% mutate(region = "downstreamIntron")
) %>%
  group_by(region, sample) %>%
  mutate(freq_perc = frequency / sum(frequency)) %>%
  ungroup() %>%
  select(kmer1, region, sample, freq_perc)

# Get the values for RBM10
MUT_RBM10_seq <- MUT_RBM10_seq %>% 
  mutate(
    upstreamIntron_kmer1_freq = map(upstreamIntronSeq_adj, ~get_kmer_frequency(.x, 1)),
    skippedExon_kmer1_freq = map(skippedExonSeq_adj, ~get_kmer_frequency(.x, 1)),
    downstreamIntron_kmer1_freq = map(downstreamIntronSeq_adj, ~get_kmer_frequency(.x, 1))
  )

# Aggregate frequency data for each region
MUT_RBM10_kmer1_freq <- bind_rows(
  MUT_RBM10_seq %>% select(ExonID, upstreamIntron_kmer1_freq) %>% unnest(upstreamIntron_kmer1_freq) %>% mutate(region = "upstreamIntron"),
  MUT_RBM10_seq %>% select(ExonID, skippedExon_kmer1_freq) %>% unnest(skippedExon_kmer1_freq) %>% mutate(region = "skippedExon"),
  MUT_RBM10_seq %>% select(ExonID, downstreamIntron_kmer1_freq) %>% unnest(downstreamIntron_kmer1_freq) %>% mutate(region = "downstreamIntron")
) %>%
  group_by(ExonID, region) %>%
  mutate(freq_perc = frequency / sum(frequency)) %>%
  ungroup() %>%
  group_by(kmer1, region) %>%
  summarise(avg_freq = mean(freq_perc), .groups = "drop")

# Plot the distribution of 1-mer frequencies for each region
ggplot(background_1mer_freq_table, aes(x = freq_perc)) +
  geom_histogram() +
  geom_vline(data = MUT_RBM10_kmer1_freq, aes(xintercept = avg_freq), color = "red", linetype = "dashed") +
  facet_wrap(~kmer1 + region, nrow = 4) + 
  labs(title = "1-mer frequency distribution for RBM10", x = "1-mer", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Check the mu and sd for bootstrapped distribution for different n.
sample_size_df_A <- data.frame(sample_size = seq(10, 200, by = 10), 
                             mu = NA, sd = NA)
sample_size_df_C <- data.frame(sample_size = seq(10, 200, by = 10),
                             mu = NA, sd = NA)
sample_size_df_G <- data.frame(sample_size = seq(10, 200, by = 10),
                             mu = NA, sd = NA)
sample_size_df_T <- data.frame(sample_size = seq(10, 200, by = 10),
                             mu = NA, sd = NA)


for (i in 1:nrow(sample_size_df_A)){
  sample_size <- sample_size_df_A$sample_size[i]
  base::print(paste0("Sample size: ", sample_size))
  
  # Compute background distribution for 1-mer in all three regions
  background_1mer_freq_upstreamIntron <- get_background_distribution_for_kmer(twist_barcodes, 1, sample_n = sample_size, mode = "upstreamIntronSeq", sample_times = 100)

  background_1mer_freq_table <- background_1mer_freq_upstreamIntron %>%
    group_by(sample) %>%
    mutate(freq_perc = frequency / sum(frequency)) %>%
    ungroup() %>%
    select(kmer1, sample, freq_perc)
  
  
  # Get the mu and sd for each region
  mu_sd_values <- background_1mer_freq_table %>%
    group_by(kmer1) %>%
    summarise(mu = mean(freq_perc), sd = sd(freq_perc), .groups = "drop")
  # 
  # > mu_sd_values
  # # A tibble: 4 × 3
  # kmer1    mu     sd
  # <chr> <dbl>  <dbl>
  #   1 A     0.225 0.0254
  # 2 C     0.232 0.0309
  # 3 G     0.208 0.0265
  # 4 T     0.334 0.0355
  
  
  # Store the values in the dataframes.
  sample_size_df_A$mu[i] <- mu_sd_values$mu[which(mu_sd_values$kmer1 == "A")]
  sample_size_df_A$sd[i] <- mu_sd_values$sd[which(mu_sd_values$kmer1 == "A")]
  sample_size_df_C$mu[i] <- mu_sd_values$mu[which(mu_sd_values$kmer1 == "C")]
  sample_size_df_C$sd[i] <- mu_sd_values$sd[which(mu_sd_values$kmer1 == "C")]
  sample_size_df_G$mu[i] <- mu_sd_values$mu[which(mu_sd_values$kmer1 == "G")]
  sample_size_df_G$sd[i] <- mu_sd_values$sd[which(mu_sd_values$kmer1 == "G")]
  sample_size_df_T$mu[i] <- mu_sd_values$mu[which(mu_sd_values$kmer1 == "T")]
  sample_size_df_T$sd[i] <- mu_sd_values$sd[which(mu_sd_values$kmer1 == "T")]
  
}

# Plot the mu and sd for each base
sample_size_df_A_longer <- sample_size_df_A %>% 
  pivot_longer(cols = c(mu, sd), names_to = "statistic", values_to = "value") %>% 
  mutate(base = "A")
ggplot(sample_size_df_A_longer, aes(x = sample_size, y = value, color = statistic)) +
  geom_line() +
  geom_point() +
  labs(title = "A base mu and sd for different sample sizes", x = "Sample size", y = "Value") +
  theme_minimal() + 
  facet_wrap(~statistic, scales = "free_y")
