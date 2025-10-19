library(tidyverse)
library(vroom)
library(data.table)

# This is the list of genes used for normalization. All of these should have PSI ~ 1 I think. 
top_normalization_seq_index <- c(
  "ENSG00000033627.17;ATP6V0A1;chr17-42501196-42501304-42500706-42500923-42507519-42507627",    
  "ENSG00000134755.18;DSC2;chr18-31093558-31093643-31092100-31092300-31101902-31102421",        
  "ENSG00000263001.7;GTF2I;chr7-74716893-74716950-74706389-74706433-74728785-74728896",         
  "ENSG00000101166.16;PRELID3B;chr20-59036689-59036760-59036470-59036573-59038465-59038634",    
  "ENSG00000121064.13;SCPEP1;chr17-56991098-56991171-56988215-56988290-56994980-56995018",      
  "ENSG00000175216.15;CKAP5;chr11-46801199-46801304-46798082-46798172-46808030-46808144",       
  "ENSG00000160190.14;SLC37A1;chr21-42535471-42535550-42525775-42525857-42539511-42539647",     
  "ENSG00000249242.8;TMEM150C;chr4-82504577-82504667-82503058-82503112-82561905-82561955",      
  "ENSG00000099204.20;ABLIM1;chr10-114441721-114441786-114441016-114441077-114444028-114444134",
  "ENSG00000134490.14;TMEM241;chr18-23399565-23399647-23377645-23377671-23421366-23421444",     
  "ENSG00000155966.14;AFF2;chrX-148842965-148843002-148809875-148809920-148843381-148843415",   
  "ENSG00000288701.1;PRRC2B;chr9-131492168-131492260-131490764-131491580-131495739-131500197",  
  "ENSG00000251138.7;LINC02882;chr12-74157789-74157884-74151001-74151050-74285024-74285093",    
  "ENSG00000158711.14;ELK4;chr1-205618956-205619073-205607942-205616644-205619965-205620838",   
  "ENSG00000165084.16;C8orf34;chr8-68815885-68815945-68787442-68787536-68818238-68818274",      
  "ENSG00000262560.1;RP11-296A16.1;chr15-43774475-43774518-43774181-43774353-43774595-43774773",
  "ENSG00000259256.2;LINC01895;chr18-3365896-3365926-3358401-3358625-3383370-3383439",          
  "ENSG00000009830.13;POMT2;chr14-77304691-77304800-77291313-77291380-77320433-77320860",       
  "ENSG00000196305.19;IARS1;chr9-92243215-92243311-92242153-92242330-92244958-92245071",        
  "ENSG00000101938.15;CHRDL1;chrX-110759660-110759754-110721384-110721530-110762694-110762807"
)

geometric_mean <- function(x) {
  exp(mean(log(x[x > 0])))  # Exclude zeros to avoid log(0)
}
normalize_counts_for_sample_PSI <- function(df, sample_name, anchor_genes_list){
  # Get the sample subset for the cell line. 
  df_subset <- df %>% filter(sample == sample_name)
  
  # Get the rows for the anchor genes. 
  anchor_genes <- df_subset %>% filter(index %in% anchor_genes_list) %>% 
    mutate(PSI = count/(count + skipped)) 
  
  # Get the geometric mean of the anchor genes.
  anchor_genes_geom_mean <- anchor_genes %>% 
    summarise(geom_mean_PSI = geometric_mean(PSI)) %>% 
    pull(geom_mean_PSI)
  
  scaling_factor <- 1/anchor_genes_geom_mean
  base::print(paste(sample_name, "scaling factor:", scaling_factor))
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
scale_all_samples_PSI <- function(df, anchor_genes_list = top_normalization_seq_index){
  unique_samples <- unique(df$sample)
  scaled_list <- lapply(unique_samples, function(x) normalize_counts_for_sample_PSI(df, x, anchor_genes_list))
  scaled_samples_df <- bind_rows(scaled_list)
  return(scaled_samples_df)
}


####### OKAY I SHOULD BE MORE SYSTEMATIC. AIYO. ########
out_dir <- "U:/processed_data/latest/updated/"
# Create outdir.
dir.create(out_dir, showWarnings = FALSE)
input_dir <- "U:/processed_data/latest/raw_47celltype/"
alt_ref_file <- read_tsv("U:/melange/data/guide_library_cleaned/ref_test_alt_ref_dict.tsv")

##### Uncomment to process files again. This takes a while. #####
cellline_filenames <- list.files(path = input_dir, pattern = "fine_grained_idx_formatted.tsv$", full.names = TRUE)

all_files <- vroom(cellline_filenames, id = "filename", delim = "\t")
# fwrite(all_files, file.path(out_dir, "umi_count_all_celltypes.csv"))
all_files_df <- all_files %>%  mutate(sample = str_extract(basename(filename), ".+(?=_umi_dedup)")) %>%
  mutate(condition = str_extract(sample, "^.+(?=-rep\\d)")) %>%
  select(-filename)
fwrite(all_files_df, file.path(out_dir, "umi_count_all_celltypes_formatted.csv"))


#### This is processing for K700E samples only. ####
input_dir <- "U:/processed_data/latest/raw_K700E/"
cellline_filenames <- list.files(path = input_dir, pattern = "fine_grained_idx_formatted.tsv$", full.names = TRUE)

all_files <- vroom(cellline_filenames, id = "filename", delim = "\t")
# fwrite(all_files, file.path(out_dir, "K700E_umi_count_all_celltypes.csv"))
all_files_df <- all_files %>%  mutate(sample = str_extract(basename(filename), ".+(?=_umi_dedup)")) %>%
  mutate(condition = str_extract(sample, "^.+(?=-H\\d)")) %>%
  select(-filename)
fwrite(all_files_df, file.path(out_dir, "K700E_umi_count_all_celltypes_formatted.csv"))
##### End processing for K700E samples only. ####

#### This is processing for Novaseq 20240826 ####
input_dir <- "U:/processed_data/latest/raw_novaseq_240826/"
cellline_filenames <- list.files(path = input_dir, pattern = "fine_grained_idx_formatted.tsv$", full.names = TRUE)

all_files <- vroom(cellline_filenames, id = "filename", delim = "\t")
# fwrite(all_files, file.path(out_dir, "K700E_umi_count_all_celltypes.csv"))
all_files_df <- all_files %>%  mutate(sample = str_extract(basename(filename), ".+(?=_umi_dedup)")) %>%
  mutate(condition = str_extract(sample, "^.+(?=-rep\\d)")) %>%
  select(-filename)
fwrite(all_files_df, file.path(out_dir, "Nova240826_umi_count_all_celltypes_formatted.csv"))
##### End processing for Novaseq 20240826. ####



replace_exact_match <- function(index_str, alt_ref_df) {
  alt_to_ref <- setNames(alt_ref_df$ref, alt_ref_df$alt)
  return(ifelse(index_str %in% names(alt_to_ref), alt_to_ref[index_str], index_str))
}

all_files_df <- fread(file.path(out_dir, "umi_count_all_celltypes_formatted.csv"))
K700E_df <- fread(file.path(out_dir, "K700E_umi_count_all_celltypes_formatted.csv"))
Novaseq240826_df <- fread(file.path(out_dir, "Nova240826_umi_count_all_celltypes_formatted.csv"))

# Apply the function to the "index" column of K700E_df
K700E_df_to_ref <- K700E_df %>%
  mutate(index = sapply(index, replace_exact_match, alt_ref_file)) %>%
  group_by(sample, condition, index, mode, offset) %>% 
  summarise(count = sum(count)) %>% 
  ungroup()
fwrite(K700E_df_to_ref, file.path(out_dir, "K700E_umi_count_merged_to_ref_normalized.csv"))

# Also apply it to the all_files_df
all_files_df_to_ref <- all_files_df %>%
  mutate(index = sapply(index, replace_exact_match, alt_ref_file)) %>%
  group_by(sample, condition, index, mode, offset) %>% 
  summarise(count = sum(count)) %>% 
  ungroup()
fwrite(all_files_df_to_ref, file.path(out_dir, "umi_count_merged_to_ref_normalized.csv"))

# Apply to the novaseq240826
Novaseq240826_df_to_ref <- Novaseq240826_df %>%
  mutate(index = sapply(index, replace_exact_match, alt_ref_file)) %>%
  group_by(sample, condition, index, mode, offset) %>% 
  summarise(count = sum(count)) %>% 
  ungroup()
fwrite(Novaseq240826_df_to_ref, file.path(out_dir, "Nova240826_umi_count_merged_to_ref_normalized.csv"))


####### Make the all_samples file that's cleaned. ########
all_files_df <- fread(file.path(out_dir, "umi_count_merged_to_ref_normalized.csv"))
# Look at sample replicates for 3ss
all_sample_reps_3ss <- all_files_df %>%
  filter(mode == "INCLUDED") %>%
  filter(!condition %in% c("OSRC2", "Kelly_old", "SKNAS_Nuc", "A172")) %>%
  separate(offset, c("offset_mid_start", "offset_mid_end", "offset_down_start"), sep = ":", convert = T, remove = F) %>%
  filter(offset_down_start == 0) %>%
  group_by(sample, index) %>%
  mutate(total_sum = sum(count)) %>%
  ungroup() %>%
  mutate(other_splice = total_sum - count)
fwrite(all_sample_reps_3ss, file = file.path(out_dir, "all_sample_reps_3ss.csv"))

# Look at sample replicates for exon skipping PSI
all_sample_reps_PSI <- all_files_df %>%
  filter((mode == "INCLUDED" & offset == "0:0:0") | (mode == "SKIPPED" & offset == "0")) %>%
  filter(!condition %in% c("OSRC2", "Kelly_old", "SKNAS_Nuc", "A172")) %>%
  group_by(sample, index) %>%
  mutate(total_sum = sum(count)) %>%
  ungroup() %>%
  mutate(skipped = total_sum - count) %>%
  filter(mode == "INCLUDED") %>%
  select(-total_sum)
all_sample_reps_PSI <- scale_all_samples_PSI(all_sample_reps_PSI)
fwrite(all_sample_reps_PSI, file = file.path(out_dir, "all_sample_reps_PSI.csv"))

##### Do the same for Nova #####
####### Make the all_samples file that's cleaned. ########
all_files_df <- fread(file.path(out_dir, "Nova240826_umi_count_merged_to_ref_normalized.csv"))
# Look at sample replicates for 5ss
all_sample_reps_3ss <- all_files_df %>%
  filter(mode == "INCLUDED") %>%
  separate(offset, c("offset_mid_start", "offset_mid_end", "offset_down_start"), sep = ":", convert = T, remove = F) %>%
  filter(offset_down_start == 0) %>%
  group_by(sample, index) %>%
  mutate(total_sum = sum(count)) %>%
  ungroup() %>%
  mutate(other_splice = total_sum - count)
# Scale the samples.
fwrite(all_sample_reps_3ss, file = file.path(out_dir, "Nova240826_reps_3ss.csv"))

# Look at sample replicates for exon skipping PSI
all_sample_reps_PSI <- all_files_df %>%
  filter((mode == "INCLUDED" & offset == "0:0:0") | (mode == "SKIPPED" & offset == "0")) %>%
  group_by(sample, index) %>%
  mutate(total_sum = sum(count)) %>%
  ungroup() %>%
  mutate(skipped = total_sum - count) %>%
  filter(mode == "INCLUDED") %>%
  select(-total_sum)
all_sample_reps_PSI <- scale_all_samples_PSI(all_sample_reps_PSI)
fwrite(all_sample_reps_PSI, file = file.path(out_dir, "Nova240826_reps_PSI.csv"))
