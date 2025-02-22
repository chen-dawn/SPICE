#### The data was downloaded from https://cellmodelpassports.sanger.ac.uk/downloads ####
library(data.table)
library(tidyverse)

sanger_gex1 <- fread("~/Downloads/rnaseq_sanger_20210316.csv")
model_id_in_gex1 <- sanger_gex1$model_id
sanger_gex2 <- fread("~/Downloads/rnaseq_broad_20210317.csv")
model_id_in_gex2 <- sanger_gex2$model_id
# remove the model_id that are in both files
overlapping_id <- intersect(model_id_in_gex1, model_id_in_gex2)
# Only remove from gex2 dataset.
sanger_gex2 <- sanger_gex2[!model_id %in% overlapping_id]

sanger_gex <- rbind(sanger_gex1, sanger_gex2)
# sanger_gex <- fread("~/Downloads/rnaseq_all_20220624/rnaseq_tpm_20220624.csv")
sanger_metadata <- fread("~/Downloads/model_list_20240110.csv") %>% 
  select(model_id, sample_id, BROAD_ID)
# > head(sanger_gex)
# model_id model_name data_source   gene_id gene_symbol read_count  fpkm
# <char>     <char>      <char>    <char>      <char>      <int> <num>
#   1: SIDM01240      451Lu      sanger SIDG00001        A1BG        127  0.19
# 2: SIDM01240      451Lu      sanger SIDG00002    A1BG-AS1       1991  4.19
# 3: SIDM01240      451Lu      sanger SIDG00003        A1CF          0  0.00
# 4: SIDM01240      451Lu      sanger SIDG00004         A2M      93979 86.50
# 5: SIDM01240      451Lu      sanger SIDG00005     A2M-AS1         35  0.09

# Change the model_id to Broad_ID
sanger_gex <- merge(sanger_gex, sanger_metadata, by = "model_id")

# For every model ID convert the read_count to TPM
sanger_gex <- sanger_gex %>%
  group_by(BROAD_ID) %>%
  mutate(tpm = read_count / sum(read_count) * 1e6) %>%
  ungroup() %>% 
  # log2(tpm + 1)
  mutate(log2_tpm = log2(tpm + 1))


broad_cellline_metadata <- fread("/Volumes/broad_dawnccle/for_anisha/cellline_data_full_metadata.csv") %>% 
  select(DepMap_ID, StrippedName) %>% distinct()

sanger_gex_wide <- sanger_gex %>% 
  select(BROAD_ID, gene_symbol, log2_tpm) %>% 
  filter(BROAD_ID %in% broad_cellline_metadata$DepMap_ID) %>%
  # Remove the gene SEPTIN4. It has duplicates
  filter(gene_symbol != "SEPTIN4") %>%
  # NA values set to 0.
  pivot_wider(names_from = gene_symbol, values_from = log2_tpm, values_fill = list(log2_tpm = 0)) %>% 
  # rename broad_id to DepMap_ID
  dplyr::rename(DepMap_ID = BROAD_ID) 

# Write this to CSV
write.csv(sanger_gex_wide, "~/Dropbox (Harvard University)/02Splicing/latest/sanger_CCLE_gex_with_Broad_ID.csv", row.names = FALSE)

