library(tidyverse)
library(vroom)
library(data.table)

# This is the list of genes used for normalization. All of these should have PSI ~ 1 I think. 
# top_normalization_seq_index <- c(
#   "ENSG00000033627.17;ATP6V0A1;chr17-42501196-42501304-42500706-42500923-42507519-42507627",    
#   "ENSG00000134755.18;DSC2;chr18-31093558-31093643-31092100-31092300-31101902-31102421",        
#   "ENSG00000263001.7;GTF2I;chr7-74716893-74716950-74706389-74706433-74728785-74728896",         
#   "ENSG00000101166.16;PRELID3B;chr20-59036689-59036760-59036470-59036573-59038465-59038634",    
#   "ENSG00000121064.13;SCPEP1;chr17-56991098-56991171-56988215-56988290-56994980-56995018",      
#   "ENSG00000175216.15;CKAP5;chr11-46801199-46801304-46798082-46798172-46808030-46808144",       
#   "ENSG00000160190.14;SLC37A1;chr21-42535471-42535550-42525775-42525857-42539511-42539647",     
#   "ENSG00000249242.8;TMEM150C;chr4-82504577-82504667-82503058-82503112-82561905-82561955",      
#   "ENSG00000099204.20;ABLIM1;chr10-114441721-114441786-114441016-114441077-114444028-114444134",
#   "ENSG00000134490.14;TMEM241;chr18-23399565-23399647-23377645-23377671-23421366-23421444",     
#   "ENSG00000155966.14;AFF2;chrX-148842965-148843002-148809875-148809920-148843381-148843415",   
#   "ENSG00000288701.1;PRRC2B;chr9-131492168-131492260-131490764-131491580-131495739-131500197",  
#   "ENSG00000251138.7;LINC02882;chr12-74157789-74157884-74151001-74151050-74285024-74285093",    
#   "ENSG00000158711.14;ELK4;chr1-205618956-205619073-205607942-205616644-205619965-205620838",   
#   "ENSG00000165084.16;C8orf34;chr8-68815885-68815945-68787442-68787536-68818238-68818274",      
#   "ENSG00000262560.1;RP11-296A16.1;chr15-43774475-43774518-43774181-43774353-43774595-43774773",
#   "ENSG00000259256.2;LINC01895;chr18-3365896-3365926-3358401-3358625-3383370-3383439",          
#   "ENSG00000009830.13;POMT2;chr14-77304691-77304800-77291313-77291380-77320433-77320860",       
#   "ENSG00000196305.19;IARS1;chr9-92243215-92243311-92242153-92242330-92244958-92245071",        
#   "ENSG00000101938.15;CHRDL1;chrX-110759660-110759754-110721384-110721530-110762694-110762807"
# )

# Get the top 100 sequences to use as normalization.
top_normalization_seq_index <- read_tsv("U:/processed_data/latest/250131_merged_v2/sequences_used_for_normalization.tsv") %>% 
  pull(index)

geometric_mean <- function(x) {
  exp(mean(log(x[x > 0])))  # Exclude zeros to avoid log(0)
}
normalize_counts_for_sample_PSI <- function(df, anchor_genes_list){
  # Get the sample subset for the cell line. 
  df_subset <- df 
  
  # Get the rows for the anchor genes. 
  anchor_genes <- df_subset %>% filter(index %in% anchor_genes_list) %>% 
    mutate(PSI = count/(count + skipped)) 
  
  # Get the geometric mean of the anchor genes.
  anchor_genes_geom_mean <- anchor_genes %>% 
    summarise(geom_mean_PSI = geometric_mean(PSI)) %>% 
    pull(geom_mean_PSI)
  
  scaling_factor <- 1/anchor_genes_geom_mean
  base::print(paste("scaling factor:", scaling_factor))
  # Scale the PSI based on the scaling factor.
  df_subset_scaled <- df_subset %>% 
    mutate(PSI = count/(count + skipped)) %>% 
    mutate(PSI_scaled = PSI * scaling_factor) %>% 
    # Set PSI above 1 to be 1.
    mutate(PSI_scaled = ifelse(PSI_scaled > 1, 1, PSI_scaled)) %>%
    # Reassign the count and skipped values based on the scaled PSI.
    mutate(count_scaled = PSI_scaled * (count + skipped),
           skipped_scaled = (1 - PSI_scaled) * (count + skipped)) %>% 
    # Make counts integers.
    mutate(count_scaled = as.integer(count_scaled),
           skipped_scaled = as.integer(skipped_scaled)) %>% 
    rename(count_original = count,
           skipped_original = skipped) %>%
    # Rename the scaled counts and skipped to count and skipped.
    rename(count = count_scaled,
           skipped = skipped_scaled)
  
  return(df_subset_scaled)
  
}

normalize_counts_individual <- function(df, anchor_genes_list){
  # This is normalizing the counts by barcode swapping rate. 
  # Get the sample subset for the cell line. 
  df_subset <- df 
  
  # Get the rows for the anchor genes. 
  anchor_genes <- df_subset %>% filter(index %in% anchor_genes_list) %>% 
    mutate(PSI = count/(count + skipped)) 
  
  # Get the geometric mean of the anchor genes.
  anchor_genes_geom_mean <- anchor_genes %>% 
    summarise(geom_mean_PSI = geometric_mean(PSI)) %>% 
    pull(geom_mean_PSI)
  
  # Swapping rate is (1-PSI)/PSI.
  swapping_rate = 1/anchor_genes_geom_mean - 1
  # E_hat = E + swap * I
  
  base::print(paste("swapping_rate:", swapping_rate))
  # Scale the PSI based on the scaling factor.
  df_subset_scaled <- df_subset %>% 
    mutate(skipped_scaled = skipped * (1 - swapping_rate)) %>% 
    mutate(count_scaled = count) %>% 
    # Make all counts integers.
    mutate(count_scaled = as.integer(count_scaled),
           skipped_scaled = as.integer(skipped_scaled)) %>%
    mutate(PSI_scaled = count_scaled/(count_scaled + skipped_scaled))
  
  return(df_subset_scaled)
  
}

scale_all_samples_PSI <- function(df, anchor_genes_list = top_normalization_seq_index){
  unique_samples <- unique(df$sample)
  scaled_list <- lapply(unique_samples, function(x) normalize_counts_for_sample_PSI(df, x, anchor_genes_list))
  scaled_samples_df <- bind_rows(scaled_list)
  return(scaled_samples_df)
}

replace_exact_match <- function(index_str, alt_ref_df) {
  alt_to_ref <- setNames(alt_ref_df$ref, alt_ref_df$alt)
  return(ifelse(index_str %in% names(alt_to_ref), alt_to_ref[index_str], index_str))
}

##### Read in the alt-ref dict. #####
alt_ref_file <- read_tsv("U:/melange/data/guide_library_cleaned/ref_test_alt_ref_dict.tsv")

####### Systematic processing. ########
out_dir <- "U:/processed_data/latest/250131_merged_v4_old_PSI_norm/"
# Create outdir.
dir.create(out_dir, showWarnings = FALSE)


# First we look at all the samples that are available. 
input_dir1 <- "U:/processed_data/latest/raw_47celltype/"
input_dir2 <- "U:/processed_data/latest/raw_novaseq_240826/"
input_dir3 <- "U:/processed_data/latest/raw_novaseq_241106/"
# input_dir4 <- "U:/processed_data/latest/raw_K700E/"

# Get all full filenames in all the input dirs. 
cellline_filenames1 <- list.files(path = input_dir1, pattern = "fine_grained_idx_formatted.tsv$", full.names = TRUE)
cellline_filenames2 <- list.files(path = input_dir2, pattern = "fine_grained_idx_formatted.tsv$", full.names = TRUE)
cellline_filenames3 <- list.files(path = input_dir3, pattern = "fine_grained_idx_formatted.tsv$", full.names = TRUE)
# cellline_filenames4 <- list.files(path = input_dir4, pattern = "fine_grained_idx_formatted.tsv$", full.names = TRUE)

# Make a df of these filenames.
cellline_paths <- data.frame(filename = c(cellline_filenames1, cellline_filenames2, cellline_filenames3))
cellline_paths <- cellline_paths %>% 
  mutate(basename = basename(filename)) %>%
  mutate(sample = str_extract(basename(filename), ".+(?=_umi_dedup)")) %>%
  mutate(condition = str_extract(sample, "^.+(?=-rep\\d)")) %>% 
  # mutate condition to strip the following: _100tfx, _150tfx, _1ugNuc, _2ugNuc.
  mutate(condition = str_replace(condition, "_\\d+tfx", "")) %>%
  mutate(condition = str_replace(condition, "_\\d+ugNuc", "")) %>% 
  # Filter out these conditions: "OSRC2", "Kelly_old", "SKNAS_Nuc"
  filter(!condition %in% c("OSRC2", "Kelly_old", "SKNAS_Nuc")) %>% 
  # Filter out filename that contain raw_47celltype/A172 or raw_47celltype/KMRC1 or 47celltype/K562_1ugNuc
  filter(!str_detect(filename, "raw_47celltype/A172|raw_47celltype/KMRC1|47celltype/K562_1ugNuc")) %>%
  # Strip SKNAS_tfx of the _tfx.
  mutate(condition = str_replace(condition, "_tfx", "")) %>% 
  # ALso extract _Nuc.
  mutate(condition = str_replace(condition, "_Nuc", "")) %>%
  # Change the DBTR05MG to DBTRG05MG
  mutate(condition = str_replace(condition, "DBTR05MG", "DBTRG05MG")) %>%
  # Change MEWO to MeWo
  mutate(condition = str_replace(condition, "MEWO", "MeWo")) %>%
  # Change JHOM to JHOM1
  mutate(condition = str_replace(condition, "JHOM", "JHOM1")) %>%
  mutate(rep_old = str_extract(basename(filename), "rep\\d")) %>%
  # Create rep new, which is rep1 = rep1, rep2 = rep2, rep3 = rep3, rep4 = rep1, rep5 = rep2, rep6 = rep3.
  mutate(rep_new = case_when(
    rep_old == "rep1" ~ "rep1",
    rep_old == "rep2" ~ "rep2",
    rep_old == "rep3" ~ "rep3",
    str_detect(basename(filename), "rep4") ~ "rep1",
    str_detect(basename(filename), "rep5") ~ "rep2",
    str_detect(basename(filename), "rep6") ~ "rep3"
  )) %>% 
  mutate(sample_new = paste0(condition, "-", rep_new))

# Get the umi_sum for each file. And save as a new column in the table. 
cellline_paths$umi_count <- apply(cellline_paths, 1, function(x) {
  file1 <- vroom(x["filename"], id = "filename", delim = "\t")
  umi_sum <- sum(file1$count)
  return(umi_sum)})

# Filter out samples that have < 10k UMI counts.
cellline_paths_filtered <- cellline_paths %>% filter(umi_count >= 10000)

# Write this metadata to file.
fwrite(cellline_paths_filtered, file.path(out_dir, "processed_samples_metadata.csv"))

# Get all unique sample_new
unique_samples <- cellline_paths_filtered %>% group_by(sample_new) %>% summarise(n=n()) %>% ungroup() 

unique_sample_names <- unique(cellline_paths_filtered$sample_new)
for (sample_tmp in unique_sample_names){
  # Get all the filepaths with that sample name.
  sample_filepaths <- cellline_paths_filtered %>% filter(sample_new == sample_tmp) %>% pull(filename)
  tmp_out <- data.frame()
  for (filepath in sample_filepaths) {
    base::print(paste("Processing", filepath))
    # Read in the tsv file.
    tmp <- vroom(filepath, id = "filename", delim = "\t")
    tmp_to_ref <- tmp %>% 
      mutate(index = sapply(index, replace_exact_match, alt_ref_file)) %>%
      group_by(index, mode, offset) %>% 
      summarise(count = sum(count)) %>% 
      ungroup()
    
    tmp_to_ref_PSI <- tmp_to_ref %>% 
      # This one was the original that only included the perfect match. 
      # filter((mode == "INCLUDED" & offset == "0:0:0") | (mode == "SKIPPED" & offset == "0")) %>%
      filter((mode == "INCLUDED" & offset == "0:0:0") | (mode == "SKIPPED" & offset == "0")) %>%
      group_by(index) %>%
      mutate(total_sum = sum(count)) %>%
      ungroup() %>%
      filter(mode == "INCLUDED") %>%
      mutate(skipped = total_sum - count) %>%
      select(-total_sum)

    
    # Normalize this file by the anchor genes.
    tmp_to_ref_PSI <- normalize_counts_for_sample_PSI(tmp_to_ref_PSI, top_normalization_seq_index)
    
    # > names(tmp_to_ref_PSI)
    # [1] "index"            "mode"            
    # [3] "offset"           "count_original"  
    # [5] "skipped_original" "PSI"             
    # [7] "PSI_scaled"       "count"           
    # [9] "skipped"    
    # Make a df of everything else that's not included and perfect. 
    tmp_everything_else <- tmp_to_ref %>% 
      filter(!(mode == "INCLUDED" & offset == "0:0:0")) %>%
      filter(!(mode == "SKIPPED" & offset == "0")) %>% 
      mutate(count_original = count) %>% 
      mutate(skipped_original = NA) %>% 
      mutate(PSI = NA) %>% 
      mutate(PSI_scaled = NA) %>% 
      mutate(skipped = NA) %>% 
      select(index, mode, offset, count_original, skipped_original, PSI, PSI_scaled, count, skipped)
  
    # Merge the dataframes. 
    tmp_final <- bind_rows(tmp_to_ref_PSI, tmp_everything_else)
    # # Make a df of skipped only. 
    # tmp_skipped_only <- tmp_to_ref_PSI %>% 
    #   mutate(count = skipped) %>% 
    #   mutate(count_scaled = skipped_scaled) %>% 
    #   mutate(mode = "SKIPPED") %>% 
    #   mutate(offset = "0") %>%
    #   select(index, mode, offset, count_original, skipped_original, count, skipped)
    # # Make a df of everything except skipped.
    # tmp_everything_else <- tmp_to_ref %>% 
    #   mutate(count_scaled = count ) %>% 
    #   filter(!(mode == "SKIPPED" & offset == "0")) %>%
    #   select(index, mode, offset, count, count_scaled)
    # 
    # tmp_final <- bind_rows(tmp_everything_else, tmp_skipped_only)
    
    tmp_out <- bind_rows(tmp_out, tmp_final)
  }
  # Group the tmp_out.
  tmp_grouped <- tmp_out %>% group_by(index, mode, offset) %>%
    summarise(count = sum(count), 
              count_original = sum(count_original), 
              skipped = sum(skipped),
              skipped_original = sum(skipped_original)) %>% 
    ungroup()
  
  # Write to outdir. 
  base::print(paste("Writing to", file.path(out_dir, paste0(sample_tmp, "_umi_dedup_normalized.tsv"))))
  fwrite(tmp_grouped, file.path(out_dir, paste0(sample_tmp, "_umi_dedup_normalized.tsv")))
  
}
