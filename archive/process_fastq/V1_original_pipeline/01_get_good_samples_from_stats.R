library(data.table)
library(tidyverse)

dir <- "U:/processed_data/celltype47_230524"
MEK1_pileup_filenames <- list.files(file.path(dir), pattern="stats_log.txt", full.names = T)

MEK1_pileup <- read_csv(MEK1_pileup_filenames[1], col_names = F) %>% 
  mutate(sample = str_extract(basename(MEK1_pileup_filenames[1]), ".+(S\\d+)")) %>% 
  mutate(condition = str_extract(sample, "^.+(?=-[:upper:]\\d+_S\\d+)"))

for (i in 2:length(MEK1_pileup_filenames)){
  print(MEK1_pileup_filenames[i])
  tmp <- read_csv(MEK1_pileup_filenames[i], col_names = F) %>% 
    mutate(sample = str_extract(basename(MEK1_pileup_filenames[i]), ".+(S\\d+)")) %>% 
    mutate(condition = str_extract(sample, "^.+(?=-[:upper:]\\d+_S\\d+)"))
  MEK1_pileup <- rbind(MEK1_pileup, tmp)
}

all_stats_log <- MEK1_pileup
names(all_stats_log) <- c("metric", "value", "sample", "condition")
write.csv(all_stats_log, "U:/processed_data/celltype47_230524/all_stats_log.txt")

stats_log_average <- all_stats_log %>% group_by(condition, metric) %>% 
  summarise(mean = mean(value), sd = sd(value))

ggplot(stats_log_average, aes(condition, mean)) + geom_bar(stat = "identity") + 
  facet_wrap(~metric, scale = "free_y") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(all_stats_log %>% filter(metric == "perc_chimera_reads"), aes(sample, value, fill = condition)) +
  geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  guides(fill = "none")

all_stats_pivot <- all_stats_log %>% pivot_wider(values_from = value, names_from = metric)
ggplot(all_stats_pivot%>% filter(total_reads > 10000), aes(total_reads, perc_chimera_reads)) + geom_point() + scale_x_log10()

good_samples <- all_stats_pivot %>% filter(total_reads > 500000 & perc_chimera_reads < 0.15)
write.csv(good_samples, "U:/processed_data/celltype47_230524/good_samples.csv")

