library(tidyverse)

# I think taking from the reads is too hard. Now we are going to jsut use the reference instead. 
# ref_alt <- read_tsv("/Volumes/broad_dawnccle/melange/data/guide_library_cleaned/ref_dict_with_alt.tsv", col_names = F)
# colnames(ref_alt) <- c("ref", "alt")
# 
# twist_lib <- read_csv("/Volumes/broad_dawnccle/melange/data/guide_library/20230130_twist_library_v3_ID_barcode_ROUT.csv")
# 
# # Remove the sequences where the ID are the ref_alt$alt.
# twist_lib_filtered <- twist_lib %>% filter(!ID %in% ref_alt$alt)
# write_csv(twist_lib_filtered, "/Volumes/broad_dawnccle/melange/data/guide_library_cleaned/20240605_twist_library_v3_ID_barcode_ROUT_filtered.csv")
# 
# guide_input_csv <- read_csv("/Volumes/broad_thechenlab/Dawn/splicing/library_sequencing_47k/20230130_twist_library_v3.csv")
# guide_input_csv_filtered <- guide_input_csv %>% filter(!ID %in% ref_alt$alt)
# write_csv(guide_input_csv_filtered, "/Volumes/broad_dawnccle/melange/data/guide_library_cleaned/20240605_twist_library_v3_filtered.csv")
# 
# 
# dlst <- ref_alt %>% filter(grepl("DLST", ref) | grepl("DLST", alt))
# dlst <- ref_alt %>% filter(grepl("DLST", alt))


ref_alt <- read_tsv("/Volumes/broad_dawnccle/melange/data/guide_library_cleaned/ref_test_multi_mapped_guides_to_discard.tsv", col_names = F)
colnames(ref_alt) <- c("ID", "count")

twist_lib <- read_csv("/Volumes/broad_dawnccle/melange/data/guide_library/20230130_twist_library_v3_ID_barcode_ROUT.csv")
guide_input_csv <- read_csv("/Volumes/broad_thechenlab/Dawn/splicing/library_sequencing_47k/20230130_twist_library_v3.csv")

# Remove the sequences where the ID are the ref_alt$alt.
twist_lib_filtered <- twist_lib %>% filter(!ID %in% ref_alt$ID)
write_csv(twist_lib_filtered, "/Volumes/broad_dawnccle/melange/data/guide_library_cleaned/20240605_twist_library_v3_ID_barcode_ROUT_filtered.csv")

guide_input_csv_filtered <- guide_input_csv %>% filter(!ID %in% ref_alt$ID)
write_csv(guide_input_csv_filtered, "/Volumes/broad_dawnccle/melange/data/guide_library_cleaned/20240605_twist_library_v3_filtered.csv")

dlst <- ref_alt %>% filter(grepl("DLST", ID))
