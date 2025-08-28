
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(scales)

#
## All conf ####
#
df_allconf <- read.table(
  "Benchmarking_results/results_offtarget_byConf_rel_by_GSV_c0_99q_summarized.txt",  # <-- replace with your actual path
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
    mean_input_reads = mean(total_input_reads),
    .groups = "drop"
  )

names(df_allconf)

# make final object adjustments
df_allconf_merged <- df_allconf %>%
  filter(read_type == "combined") %>%
  #filter(k_col == "k_9" ) %>%
  mutate(
    Conf = factor(
      Conf,
      levels = c("conf0","conf0.1","conf0.2","conf0.3","conf0.4","conf0.5")
    )
  )

# This adds up the psv counts for whatever made it past the filter
df_allconf_merged <- df_allconf_merged %>% 
  mutate(
    filtered_num_aligned = str_extract_all(          # pull every “… : …” pair
      filtered_psv_groups,
      "\\d+\\s*:\\s*([0-9.]+(?:[eE][+-]?[0-9]+)?)"
    ) %>%                                            # list-column of strings
      map_dbl(~ {                                    # → numeric vector per row
        if (length(.x) == 0) return(0)               # keep 0 when no pairs
        str_extract(.x, "(?<=:)\\s*[0-9.]+(?:[eE][+-]?[0-9]+)?") %>% 
          as.numeric() %>% 
          sum(na.rm = TRUE)
      })
  )

## Start of summary stats ####

df_allconf_merged %>% 
  filter(read_type == "combined") %>% 
  group_by(Conf) %>% 
  filter(k_col == "k_9" ) %>% 
  summarise(
    num_aligned_mean = mean(num_aligned, na.rm = TRUE),
    #num_aligned_sd   = sd(num_aligned,  na.rm = TRUE),
    filtered_num_aligned_mean = mean(filtered_num_aligned, na.rm = TRUE),
    num_samples      = n(),
    num_passed       = sum(status == "passed", na.rm = TRUE),
    num_filtered = sum(status == "filtered", na.rm = TRUE),
    percent_FP   = (num_passed / 500) * 100,
    .groups = "drop"
  )

df_allconf_merged %>% 
  filter(read_type == "combined") %>% 
  filter(status == "passed")

df_allconf_merged %>% 
  filter(read_type == "combined") %>% 
  filter(k_col == "k_9" ) %>%
  group_by(Conf) %>% 
  summarise(
    #num_aligned_mean = mean(num_aligned, na.rm = TRUE),
    #proportion_aligned = num_aligned/total_input_reads  ,
    mean_proportion_aligned = mean(num_aligned/total_input_reads, na.rm = TRUE),
    num_samples      = n(),
    num_passed       = sum(status == "passed", na.rm = TRUE),
    num_filtered = sum(status == "filtered", na.rm = TRUE),
    percent_FP   = (num_passed / num_samples) * 100,
    .groups = "drop"
  )


## 1 · build a per-sample proportion column  ---------------------------
df_plot <- df_allconf_merged %>% 
  filter(read_type == "combined") %>%           # one row per sample
  mutate(prop_aligned = num_aligned / total_input_reads)

## 2 · box-and-whisker plot by Conf  -----------------------------------
ggplot(df_plot, aes(x = Conf, y = prop_aligned)) +
  geom_boxplot(outlier.shape = 21, fill = "grey80") +
  scale_y_continuous(
    #labels = scales::percent_format(accuracy = 0.01)   # show to 0.01 %
    # or, if you prefer raw proportions:
    labels = scales::number_format(accuracy = 0.001)
  ) +
  labs(x = "Confidence threshold",
       y = "Proportion of reads aligned",
       title = "FP read proportion by kraken confidence level") +
  theme_classic()




### Different figure #####
## 0 · data with per–sample proportion
df_plot <- df_allconf_merged %>% 
  filter(read_type == "combined") %>% 
  mutate(prop_aligned = num_aligned / total_input_reads)

## 1 · mean per confidence level  (one row per Conf)
df_mean <- df_plot %>% 
  group_by(Conf) %>% 
  summarise(y = mean(prop_aligned, na.rm = TRUE), .groups = "drop")

## 2 · plot
ggplot(df_plot, aes(x = Conf, y = prop_aligned, fill = Conf)) +
  geom_boxplot(outlier.shape = 21, alpha = .65, width = .7, colour = "grey40") +
  
  ## red horizontal bar at the mean
  geom_segment(data = df_mean,
               aes(x = as.numeric(Conf) - .35,
                   xend = as.numeric(Conf) + .35,
                   y  = y, yend = y),
               colour = "red", linewidth = 1.1) +
  
  ## numeric label just above the bar
  geom_text(data = df_mean,
            aes(x = Conf, y = y, label = number(y, accuracy = 0.0001)),
            colour = "red", vjust = -0.6, size = 3) +
  
  scale_y_continuous(labels = number_format(accuracy = 0.001)) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(x = "Confidence threshold",
       y = "Proportion of reads aligned",
       title = "FP read proportion by Kraken confidence level") +
  theme_classic()


## per–sample proportion (as before)
df_plot <- df_allconf_merged %>% 
  filter(read_type == "combined") %>% 
  mutate(prop_aligned = num_aligned / total_input_reads)

## summary: mean, total n, and n above the mean
df_above_mean <- df_plot %>% 
  group_by(Conf) %>% 
  mutate(mean_conf = mean(prop_aligned, na.rm = TRUE)) %>%   # same mean for every row in group
  summarise(
    mean_prop   = mean_conf[1],                    # the group mean itself
    n_total     = n(),                             # number of samples in this Conf
    n_above_mean = sum(prop_aligned > mean_conf),  # samples whose value > mean
    .groups = "drop"
  )

df_above_mean

505/2320 * 100 # Percent samples with FP greater than average
