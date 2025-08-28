
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(scales)

#
## All conf ####
#
df_allconf <- read.table(
  "Benchmarking_results/results_offtarget_byConf_sub_count_rel_abund_0.0004_summarized.txt",  # <-- replace with your actual path
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

df_allconf <- df_allconf %>%
  mutate(
    Conf = sub(".*(conf[^_/]*)$", "\\1", sub("\\.k_\\d+$", "", file))
  )

str(df_allconf)
unique(df_allconf$status)

# Summarize by read type
df_allconf %>%
  group_by(read_type, numPSVs) %>%
  summarise(
    num_aligned_mean = mean(num_aligned, na.rm = TRUE),
    num_aligned_sd   = sd(num_aligned, na.rm = TRUE),
    .groups = "drop"
  )

df_allconf

# make final object adjustments
df_allconf_merged <- df_allconf %>%
  filter(read_type == "combined" ) %>%
  filter(k_col == "k_9" ) %>%
  mutate(
    Conf = factor(
      Conf,
      levels = c("conf0","conf0.1","conf0.2","conf0.3","conf0.4","conf0.5")
    )
  )


df_allconf_merged %>%
  group_by(Conf, numPSVs) %>%
  summarise(
    num_aligned_mean = mean(num_aligned, na.rm = TRUE),
    num_aligned_sd   = sd(num_aligned, na.rm = TRUE),
    num_samples = n(),
    .groups = "drop"
  )

df_allconf_merged %>% 
  group_by(Conf) %>% 
  summarise(
    num_aligned_mean = mean(num_aligned, na.rm = TRUE),
    num_aligned_sd   = sd(num_aligned,  na.rm = TRUE),
    num_samples      = n(),
    num_passed       = sum(status == "passed", na.rm = TRUE),
    num_filtered = sum(status == "filtered", na.rm = TRUE),
    percent_FP   = (num_passed / num_samples) * 100,
    .groups = "drop"
  )

