library(tidyverse)
library(data.table)
library(magrittr)
library(ggplot2)
library(vroom)
library(ggpointdensity)
library(pheatmap)
library(gridExtra)
library(RColorBrewer)
library(pheatmap)
library(gridExtra)
library(dplyr)
library(tibble)



# Define colors
color_palette2 <- c(
  "#4575B4",  # deep blue
  "#85B6D6",  # slightly lighter/more even blue
  "#E2EFF2",  # less stark pastel blue
  "#FFE3B0",  # warmer, slightly less saturated yellow
  "#EF9651",  # softer orange
  "#D83629"   # red
)
color_palette2_custom <- colorRampPalette(color_palette2)(100)
color_palette2_custom_rev <- colorRampPalette(rev(color_palette2))(100)


# Set path where images are stored and navigate to directory.
output_path <- "/Volumes/broad_chen_imaging/Dawn/241110_pSLED_Ef1a_all_D3/"
setwd(output_path)

plate1_sample_sheet_path <- "/Volumes/broad_chen_imaging/Dawn/241110_pSLED_Ef1a_all_D3/plate1_sample_sheet.csv"
plate2_sample_sheet_path <- "/Volumes/broad_chen_imaging/Dawn/241110_pSLED_Ef1a_all_D3/plate2_sample_sheet.csv"
plate3_sample_sheet_path <- "/Volumes/broad_chen_imaging/Dawn/241110_pSLED_Ef1a_all_D3/plate3_sample_sheet.csv"
plate1_sample_sheet <- read.csv(plate1_sample_sheet_path) %>% mutate(well = str_extract(filename, "(?<=\\d{2}_)[A-Z]\\d{2}")) %>% select(-filename) %>% mutate(well = paste0("plate1", "_", well))
plate2_sample_sheet <- read.csv(plate2_sample_sheet_path) %>% mutate(well = str_extract(filename, "(?<=\\d{2}_)[A-Z]\\d{2}")) %>% select(-filename) %>% mutate(well = paste0("plate2", "_", well))
plate3_sample_sheet <- read.csv(plate3_sample_sheet_path) %>% mutate(well = str_extract(filename, "(?<=\\d{2}_)[A-Z]\\d{2}")) %>% select(-filename) %>% mutate(well = paste0("plate3", "_", well))
sample_sheet <- rbind(plate1_sample_sheet, plate2_sample_sheet, plate3_sample_sheet)

# Constants
min_mcherry <- 130
min_gfp <- 117
min_bfp <- 120

######### Functions ##########
# Background subtraction function
subtract_background_from_plate <- function(plate) {
  plate %>%
    mutate(
      mc_no_background = pmax(mcherry - min_mcherry, 0),
      gfp_no_background = pmax(gfp - min_gfp, 0),
      bfp_no_background = pmax(bfp - min_bfp, 0)
    )
}

# Aggregation function
aggregate_values <- function(plate_bg_subtracted) {
  plate_bg_subtracted %>%
    group_by(well) %>%
    summarise(
      mcherry_sum = sum(mc_no_background * area),
      gfp_sum = sum(gfp_no_background * area),
      bfp_sum = sum(bfp_no_background * area),
      mean_mcherry = mean(mc_no_background * area / sum(area)),
      mean_gfp = mean(gfp_no_background * area / sum(area)),
      mean_bfp = mean(bfp_no_background * area / sum(area)),
      area_sum = sum(area),
      mcherry_gfp_ratio = mean_mcherry / mean_gfp,
      gfp_mcherry_ratio = mean_gfp / mean_mcherry,
      mcherry_bfp_ratio = mean_mcherry / mean_bfp,
      gfp_bfp_ratio = mean_gfp / mean_bfp
    )
}

######### Process images and make overview PSI figures ########
# Cell lines and plate numbers
cell_lines <- c("HEK", "TOV21G", "JHH6", "MCF7", "KELLY", "HCC38", "T47D", "GI1")
plates <- 1:3

# Generate one aggregated variable per cell line
aggregated_data <- list()
for (cell in cell_lines) {
  cell_data <- bind_rows(lapply(plates, function(plate_num) {
    input_parent_path <- file.path(output_path, paste0(cell, "_plate", plate_num))
    multi_plate_image_filenames <- list.files(input_parent_path, pattern = "\\.csv$")
    merged_csvs <- vroom(file.path(input_parent_path, multi_plate_image_filenames), id = "filename") %>%
      mutate(
        filename = basename(filename),
        well = str_extract(filename, "(?<=\\d{2}_)[A-Z]\\d{2}"),
        well = paste0("plate", plate_num, "_", well)
      ) %>%
      select(-filename)
    plate_bg_subtracted <- subtract_background_from_plate(merged_csvs)
    aggregate_values(plate_bg_subtracted)
  }))
  assign(paste0(cell, "_aggregated"), cell_data %>% mutate(celltype = cell))
  aggregated_data[[cell]] <- cell_data
}

all_aggregated <- rbind(
  HEK_aggregated, T47D_aggregated, JHH6_aggregated, MCF7_aggregated, 
  HCC38_aggregated, KELLY_aggregated, GI1_aggregated, TOV21G_aggregated
) %>% 
  filter(area_sum > 10000)

# Merge with sample sheet and adjust cell type factor levels
all_aggregated <- merge(all_aggregated, sample_sheet, by = "well")
all_aggregated$celltype <- factor(all_aggregated$celltype, levels = c("HEK", "T47D", "JHH6", "MCF7", "HCC38", "KELLY", "GI1", "TOV21G"))

# filter out cells that have low GFP/BFP and low mCherry/BFP raito
all_aggregated <- all_aggregated %>%
  filter(mcherry_bfp_ratio >= 0.1 | gfp_bfp_ratio >=0.1)

# Compute summary statistics for each cell type and condition
all_aggregated_avg <- all_aggregated %>% 
  group_by(celltype, condition) %>% 
  summarise(
    mean_mcherry = mean(mean_mcherry), 
    mean_gfp = mean(mean_gfp), 
    mean_bfp = mean(mean_bfp), 
    mcherry_gfp_ratio = mean(mcherry_gfp_ratio), 
    gfp_mcherry_ratio = mean(gfp_mcherry_ratio),
    sd_mcherry = sd(mean_mcherry), 
    sd_gfp = sd(mean_gfp), 
    sd_bfp = sd(mean_bfp), 
    sd_mcherry_gfp_ratio = sd(mcherry_gfp_ratio), 
    sd_gfp_mcherry_ratio = sd(gfp_mcherry_ratio)
  )

# Convert to wide format for SLED data
all_aggregated_wide <- all_aggregated_avg %>% 
  mutate(celltype = paste0(celltype, "_SLED")) %>%
  select(celltype, condition, gfp_mcherry_ratio) %>% 
  pivot_wider(names_from = celltype, values_from = gfp_mcherry_ratio, values_fill = NA)

# Read in the validation labels
cloning_sheet <- read_csv("/Volumes/broad_dawnccle/melange/data/twist_seq_shortlist.csv") %>%
  mutate(condition = str_extract(DCPlasmidNum2, "DC\\d+")) %>% 
  select(condition, ExonID, Mod3)

# Read in PSI values and process
all_sample_reps <- fread("/Volumes/broad_dawnccle/melange/data/250527_FINAL_PSI_TABLE_by_condition.csv")
all_sample_reps <- all_sample_reps %>%
  separate(index_offset, into = c("index", "offset"), sep = "__", remove = FALSE) %>%
  filter(offset == "0:0:0")

  
all_samples_wide <- all_sample_reps %>%
  pivot_wider(names_from = condition, values_from = PSI, values_fill = NA) %>%
  select(-index_offset, -offset)

# Filter PSI for specific columns and merge with cloning sheet
measured_PSI <- all_samples_wide %>%
  select(index, HEK, T47D, JHH6, MCF7, HCC38, KELLY, GI1, TOV21G) %>%
  filter(index %in% cloning_sheet$ExonID) %>%
  merge(cloning_sheet, by.x = "index", by.y = "ExonID") %>%
  arrange(condition)

# Prepare PSI matrix
PSI_mat <- measured_PSI %>%
  select(HEK, T47D, JHH6, MCF7, HCC38, KELLY, GI1, TOV21G) %>%
  as.matrix()
rownames(PSI_mat) <- measured_PSI$index

# Prepare SLED data and align it with PSI data
measured_PSI_SLED <- merge(measured_PSI, all_aggregated_wide, by = "condition") %>%
  mutate(across(ends_with("_SLED"), ~ ifelse(Mod3 == 1, 1 / ., .)))

# Remove "_SLED" suffix from SLED columns for alignment
sled_columns <- names(measured_PSI_SLED)[grepl("SLED$", names(measured_PSI_SLED))]
sled_mat <- measured_PSI_SLED %>%
  select(all_of(sled_columns)) %>%
  as.matrix()
rownames(sled_mat) <- measured_PSI_SLED$index
colnames(sled_mat) <- gsub("_SLED$", "", colnames(sled_mat))

# Ensure both matrices have the same rows and columns
common_rows <- intersect(rownames(PSI_mat), rownames(sled_mat))
common_columns <- union(colnames(PSI_mat), gsub("_SLED$", "", colnames(sled_mat)))

# Reorder rows and columns to ensure alignment
PSI_mat_aligned <- PSI_mat[match(common_rows, rownames(PSI_mat)), match(common_columns, colnames(PSI_mat), nomatch = NA)]
colnames(PSI_mat_aligned) <- common_columns
rownames(PSI_mat_aligned) <- common_rows

# Align rows and columns for SLED data
sled_mat_aligned <- sled_mat[match(common_rows, rownames(sled_mat)), ]
colnames(sled_mat_aligned) <- gsub("_SLED$", "", colnames(sled_mat_aligned))
rownames(sled_mat_aligned) <- common_rows


######### Normalize frameshift results by known controls ####


psi_1_controls <- c("ENSG00000173175.15;ADCY5;chr3-123314234-123314322-123303054-123303219-123319673-123319818",
                    "ENSG00000204308.8;RNF5;chr6-32180008-32180126-32179876-32179931-32180228-32180793",
                    "ENSG00000009830.13;POMT2;chr14-77304691-77304800-77291313-77291380-77320433-77320860")


psi_0_controls <-c("ENSG00000177830.18;CHID1;chr11-908214-908296-904705-904859-908541-908645",
                   "ENSG00000174669.12;SLC29A2;chr11-66366125-66366231-66365935-66366021-66366430-66366564",
                   "ENSG00000236532.6;LINC01695;chr21-28167225-28167319-28163601-28163717-28170077-28170143")

# Create a new dataframe to store normalized values
sled_mat_aligned_df <- as.data.frame(sled_mat_aligned) 
normalized_sled_df <- sled_mat_aligned_df

# Get row indices for PSI=1 controls
psi_1_rows <- which(rownames(sled_mat_aligned_df) %in% psi_1_controls)

# Get row indices for PSI=0 controls
psi_0_rows <- which(rownames(sled_mat_aligned_df) %in% psi_0_controls)


# Small constant to avoid log(0)
epsilon <- 1e-6

for (col_name in colnames(sled_mat_aligned_df)) {
  sled_col <- sled_mat_aligned_df[[col_name]]
  # Add epsilon to avoid log of zero or negative numbers
  sled_col_adjusted <- sled_col + epsilon
  # Take the natural logarithm of the column
  log_sled_col <- log2(sled_col_adjusted)
  # Calculate mean values for controls on the log-transformed data
  mean_log_psi_1 <- mean(log_sled_col[psi_1_rows], na.rm = TRUE)
  mean_log_psi_0 <- mean(log_sled_col[psi_0_rows], na.rm = TRUE)
  # Perform normalization
  normalized_sled_col <- (log_sled_col - mean_log_psi_0) / (mean_log_psi_1 - mean_log_psi_0)

  # Assign back to the dataframe
  normalized_sled_df[[col_name]] <- normalized_sled_col
}

# re-scale between 0 and 1
normalized_sled_df <- apply(normalized_sled_df, 2, function(x) {
  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x, na.rm = TRUE)
  (x - x_min) / (x_max - x_min)
})

# Convert back to data frame if you started with a data frame
normalized_sled_df <- as.data.frame(normalized_sled_df)
normalized_sled_mat <- as.matrix(normalized_sled_df)

colnames(normalized_sled_mat) <- gsub("_SLED$", "", colnames(sled_mat_aligned))
rownames(normalized_sled_mat) <- common_rows



######### Filter to smaller to display #######

# Identify the indices of the selected rows in cloning_sheet
selected_row_labels <- c("ENSG00000262560.1;RP11-296A16.1;chr15-43774475-43774518-43774181-43774353-43774595-43774773", 
                         "ENSG00000134490.14;TMEM241;chr18-23399565-23399647-23377645-23377671-23421366-23421444",
                         "ENSG00000155966.14;AFF2;chrX-148842965-148843002-148809875-148809920-148843381-148843415",
                         "ENSG00000132388.13;UBE2G1;chr17-4331773-4331888-4307020-4307123-4366270-4366628",
                         "ENSG00000028203.18;VEZT;chr12-95258244-95258282-95257149-95257239-95262905-95263081",
                         "ENSG00000126777.18;KTN1;chr14-55646972-55647007-55641691-55641760-55648024-55648088",
                         "ENSG00000171503.13;ETFDH;chr4-158681570-158681614-158672306-158672490-158682194-158682424",
                         "ENSG00000105671.12;DDX49;chr19-18921662-18921748-18920579-18920703-18921842-18921964",
                         "ENSG00000138821.14;SLC39A8;chr4-102285948-102285973-102267871-102268079-102304316-102304481",
                         "ENSG00000112818.11;MEP1A;chr6-46793551-46793576-46793388-46793467-46793665-46793716",
                         "ENSG00000135365.16;PHF21A;chr11-45946075-45946098-45938156-45938312-45948885-45948946",
                         "ENSG00000177483.12;RBM44;chr2-237821345-237821368-237820174-237821254-237821742-237821827")



# Filter and reorder PSI_final
PSI_final <- as.data.frame(PSI_mat_aligned) %>%
  rownames_to_column(var = "ExonID") %>% # Convert row names to a column named 'ExonID'
  filter(ExonID %in% selected_row_labels) %>%
  select(-GI1, -TOV21G, -T47D) 


# Filter and reorder SLED_final
SLED_final <- as.data.frame(normalized_sled_mat) %>%
  rownames_to_column(var = "ExonID") %>% # Convert row names to a column named 'ExonID'
  filter(ExonID %in% selected_row_labels) %>%
  select(-GI1, -TOV21G, -T47D) 


# Convert to matrix and handle missing values
PSI_matrix <- PSI_final %>%
  as.data.frame() %>%                 
  mutate(ExonID = as.character(ExonID)) %>%  
  column_to_rownames("ExonID") %>%    
  as.matrix()                         

SLED_matrix <- SLED_final %>%
  as.data.frame() %>%                 
  mutate(ExonID = as.character(ExonID)) %>%  
  column_to_rownames("ExonID") %>%    
  as.matrix()                         

# Re-order rows
desired_order <-  c(10, 11, 12, 4, 1, 5, 6, 8, 7, 2, 3, 9) 
PSI_matrix <- PSI_matrix[desired_order, , drop = FALSE]
SLED_matrix <- SLED_matrix[desired_order, , drop = FALSE]


piyg_colors <- colorRampPalette(rev(brewer.pal(11, "PiYG")))(100)

# Generate the heatmaps
PSI_heatmap <- pheatmap(
  PSI_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "Sequencing PSI values",
  show_rownames = FALSE,
  silent = TRUE,
  color = color_palette2_custom
)

SLED_normalized_heatmap <- pheatmap(
  SLED_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "Normalized SLED values",
  show_rownames = TRUE,
  silent = TRUE,
  color = piyg_colors 
)

# Save both to a PDF
pdf("/Volumes/broad_dawnccle/test/fig02_frameshift_heatmap.pdf", width = 12, height = 6)  
grid.arrange(PSI_heatmap$gtable, SLED_normalized_heatmap$gtable, ncol = 2, widths = c(1, 3))
dev.off()

######### Calculate the correaltion values between the PSI and sled values #########
# Flatten matrices into vectors for comparison
PSI_values <- as.vector(PSI_matrix)
SLED_values <- as.vector(SLED_matrix)

# Combine into a dataframe
comparison_df <- data.frame(
  PSI = PSI_values,
  SLED = SLED_values
)

# Calculate correlations
pearson_cor <- cor(comparison_df$PSI, comparison_df$SLED, method = "pearson")
spearman_cor <- cor(comparison_df$PSI, comparison_df$SLED, method = "spearman")

# Print correlation results
cat("Pearson correlation:", round(pearson_cor, 3), "\n")
cat("Spearman correlation:", round(spearman_cor, 3), "\n")



