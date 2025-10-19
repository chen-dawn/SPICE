library(data.table)
library(tidyverse)

dir <- "U:/processed_data/celltype47_230524/umi_dedup/"
MEK1_pileup_filenames <- list.files(file.path(dir), pattern="umi_dedup.csv", full.names = T)

MEK1_pileup <- read_csv(MEK1_pileup_filenames[1]) %>% 
  mutate(sample = str_extract(basename(MEK1_pileup_filenames[1]), ".+(S\\d+)")) %>% 
  mutate(condition = str_extract(sample, "^.+(?=-[:upper:]\\d+_S\\d+)"))

for (i in 2:length(MEK1_pileup_filenames)){
  print(MEK1_pileup_filenames[i])
  tmp <- read_csv(MEK1_pileup_filenames[i]) %>% 
    mutate(sample = str_extract(basename(MEK1_pileup_filenames[i]), ".+(S\\d+)")) %>% 
    mutate(condition = str_extract(sample, "^.+(?=-[:upper:]\\d+_S\\d+)"))
  MEK1_pileup <- rbind(MEK1_pileup, tmp)
}

all_umi_dedup <- MEK1_pileup
write.csv(all_umi_dedup, "U:/processed_data/celltype47_230524/umi_dedup/umi_dedup_all.csv")

good_samples <- read_csv("U:/processed_data/celltype47_230524/good_samples.csv")

good_umi_dedup <- all_umi_dedup %>% filter(sample %in% good_samples$sample)
bad_umi_dedup <- all_umi_dedup %>% filter(!(sample %in% good_samples$sample))


test <- good_umi_dedup %>% 
  filter((included + skipped) >20) %>%
  group_by(sample, condition) %>% summarise(n = n())#  %>% filter(grepl("K562", sample))
ggplot(test, aes(sample, n, fill = condition)) + geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + guides(fill = "none") + 
  ylab("Num Good Cell Barcodes Detected (at least 20 UMIs)")

# Change the names of the samples.
good_umi_dedup_clean <- good_umi_dedup %>% filter(condition != "K562_K700E")
good_umi_dedup_clean$condition[which(grepl("K562", good_umi_dedup_clean$condition))] <- "K562"
good_umi_dedup_clean$condition[which(grepl("A375", good_umi_dedup_clean$condition))] <- "A375"

names(good_umi_dedup_clean)[1] <- "cb"
write.csv(good_umi_dedup_clean, "U:/processed_data/celltype47_230524/umi_dedup/umi_dedup_good_samples_cleaned.csv")




HEK <- good_umi_dedup %>% filter(sample %in% c("HEK-A02_S194", "HEK-A03_S195", "K562-C10_S226", "K562_WT-H01_S259")) %>% 
  filter((included + skipped) >30) %>%
  mutate(PSI = included/(included + skipped)) 

names(HEK)[1] <- "cb"
HEK<- HEK%>% select(cb, PSI, sample) %>% pivot_wider(values_from = PSI, names_from = sample)

ggplot(HEK, aes(`HEK-A03_S195`, `K562-C10_S226`)) + geom_point(alpha = 0.3)

