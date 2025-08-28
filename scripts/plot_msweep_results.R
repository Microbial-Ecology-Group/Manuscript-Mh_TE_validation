# Plot msweep results

#!/usr/bin/env Rscript

# ----------------------------------------------------
# Extended R script to produce:
#   1) A plot with Precision & Recall by k_col + numPSVs
#      (multiple lines: one for each (Precision/Recall, numPSVs) combo).
#   2) A plot with Precision & Recall by k_col only
#      (aggregated across all numPSVs => single line for Precision, single line for Recall).
#
# Final ggplot objects => 'p1' and 'p2'
# ----------------------------------------------------

# (Optional) install required libraries:
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("tidyr")
# install.packages("stringr")
# install.packages("scales") # if needed for colorRampPalette

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(scales)

# --------------------------------------------------
# 1) Read the summarized file
# --------------------------------------------------
df <- read.table(
  "Benchmarking_results/results_100i_0.25conf_summarized.txt",  # <-- replace with your actual path
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

#results_900iters_ntcore_0.5conf_summarized.txt
#results_900iters_0.001abund_ntcore_0.5conf_summarized.txt
#results_900iters_10count_ntcore_0.5conf_summarized.txt
#results_900iters_100count_ntcore_0.5conf_summarized.tx



# --------------------------------------------------
# 2) (Optional) Filter out rows where k_col == "off_target"
# --------------------------------------------------
df_filt <- df %>%
  filter(read_type == "merged")

df_filt <- df_filt %>%
  filter(k_col != "off_target")

# Also remove cases where numPSVs is greater than k_col
nrow(df_filt) # 55345

df_filt <- df_filt %>%
  mutate(
    nPSVs = as.numeric(str_extract(numPSVs, "\\d+")),
    nK    = as.numeric(str_extract(k_col, "\\d+"))
  ) %>%
  filter(nK >= nPSVs)  

nrow(df_filt) #49609

# --------------------------------------------------
# 3) FIGURE 1:
#    Precision & Recall by (k_col, numPSVs).
#    Each combination => separate line.
# --------------------------------------------------

df_filt %>%
  group_by(numPSVs) %>%
  summarise(
    num_aligned_mean = mean(num_aligned, na.rm = TRUE),
    num_aligned_sd   = sd(num_aligned, na.rm = TRUE),
    .groups = "drop"
  )


# 3a) Group by k_col, numPSVs => compute mean & SD for Precision, Recall
df_summary <- df_filt %>%
  group_by(k_col, numPSVs) %>%
  summarise(
    Precision_mean = mean(Precision, na.rm = TRUE),
    Precision_sd   = sd(Precision, na.rm = TRUE),
    Recall_mean    = mean(Recall, na.rm = TRUE),
    Recall_sd      = sd(Recall, na.rm = TRUE),
    .groups = "drop"
  )

# 3b) Convert 'k_col' to numeric so we can plot on x-axis in ascending order
df_summary <- df_summary %>%
  mutate(k_col_num = as.numeric(str_extract(k_col, "\\d+"))) %>%
  arrange(k_col_num, numPSVs)

# 3c) Pivot to a "long" format for easy plotting of multiple metrics (Precision, Recall)
df_long <- df_summary %>%
  pivot_longer(
    cols         = c(Precision_mean, Precision_sd, Recall_mean, Recall_sd),
    names_to     = c("metric", "stat"),   # e.g. metric="Precision", stat="mean"
    names_pattern= "^(.*)_(mean|sd)$",    # parse "Precision_mean" => metric=Precision, stat=mean
    values_to    = "value"
  ) %>%
  pivot_wider(
    names_from  = stat,    # "mean" vs "sd"
    values_from = value
  )
# df_long has columns: k_col, numPSVs, k_col_num, metric, mean, sd

# 3d) Build a color map so each (metric, numPSVs) combo has a distinct color
unique_numpsvs <- sort(unique(df_long$numPSVs))

# Two color gradients (optional):
purple_palette <- colorRampPalette(c("#E6CCFF", "#6B00B6"))(length(unique_numpsvs))
yellow_palette <- colorRampPalette(c("#FFF7CC", "#B3A300"))(length(unique_numpsvs))

make_key <- function(metric, psv) paste(metric, psv, sep=":")

color_map <- c()
for (i in seq_along(unique_numpsvs)) {
  psv_val <- unique_numpsvs[i]
  color_map[ make_key("Precision", psv_val) ] <- purple_palette[i]
  color_map[ make_key("Recall",    psv_val) ] <- yellow_palette[i]
}

df_long <- df_long %>%
  mutate(combo = paste(metric, numPSVs, sep=":"))

# 3e) Plot => multiple lines
p1 <- ggplot(df_long, aes(x = k_col_num, y = mean, group = combo, color = combo)) +
  geom_line() +
  geom_point() +
  #geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
  
  # Use numeric x for spacing, but label with original k_col
  scale_x_continuous(
    breaks = df_long$k_col_num,
    labels = df_long$k_col
  ) +
  scale_color_manual(values = color_map) +
  
  labs(
    x = "k_col",
    y = "Value",
    title = "Figure 1: Precision & Recall by k_col + numPSVs (100 iterations)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))

print(p1)

# --------------------------------------------------
# 4) FIGURE 2:
#    Single line for Precision & Recall across all numPSVs
#    i.e., aggregated by k_col only.
# --------------------------------------------------

# 4a) Group by k_col => compute mean & SD for Precision, Recall (pooling all numPSVs)
df_agg <- df_filt %>%
  group_by(k_col) %>%
  summarise(
    Precision_mean = mean(Precision, na.rm = TRUE),
    Precision_sd   = sd(Precision,   na.rm = TRUE),
    Recall_mean    = mean(Recall,    na.rm = TRUE),
    Recall_sd      = sd(Recall,      na.rm = TRUE),
    .groups = "drop"
  )

# 4b) Convert 'k_col' to numeric
df_agg <- df_agg %>%
  mutate(k_col_num = as.numeric(str_extract(k_col, "\\d+"))) %>%
  arrange(k_col_num)

# 4c) Plot: one line for Precision, one for Recall, each with SD error bars
# We'll do two separate geom_line() calls with fixed colors.
p2_filt <- ggplot(df_agg, aes(x = k_col_num)) +
  # Precision line + points + error bars
  geom_line(aes(y = Precision_mean, color = "Precision")) +
  geom_point(aes(y = Precision_mean, color = "Precision")) +
  #geom_errorbar(aes(
  #  ymin = Precision_mean - Precision_sd,
  #  ymax = Precision_mean + Precision_sd,
  #  color = "Precision"
  #), width = 0.2) +
  
  # Recall line + points + error bars
  geom_line(aes(y = Recall_mean, color = "Recall")) +
  geom_point(aes(y = Recall_mean, color = "Recall")) +
  #geom_errorbar(aes(
  #  ymin = Recall_mean - Recall_sd,
  #  ymax = Recall_mean + Recall_sd,
  #  color = "Recall"
  #), width = 0.2) +
  
  # Label the x-axis with the original k_col
  scale_x_continuous(
    breaks = df_agg$k_col_num,
    labels = df_agg$k_col
  ) +
  
  # Manual color assignments for the 2 lines
  scale_color_manual(
    name   = "Metric",
    values = c("Precision" = "#6B00B6", "Recall" = "#B3A300")
  ) +
  
  labs(
    x = "k_col",
    y = "Value",
    title = "Figure 2: Precision & Recall by k_col (Aggregated over numPSVs, 100 iterations)"
  ) +
  
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))

# --------------------------------------------------
# 5) Print or save the plots
# --------------------------------------------------
print(p1)
print(p2_filt)



## 

df_agg

df_agg_conf0 <- df_filt
df_agg_conf0$Conf <- "0 conf"

df_agg_conf0.25 <- df_filt
df_agg_conf0.25$Conf <- "0.25 conf"

df_agg_conf0.5 <- df_filt
df_agg_conf0.5$Conf <- "0.5 conf"


combined_filtering_results <- rbind(df_agg_conf0,df_agg_conf0.25,df_agg_conf0.5)

write.csv(combined_filtering_results, "combined_confidence_results.txt")
# Optionally:
# ggsave("figure1_by_numPSVs.png", p1, width = 8, height = 6)
# ggsave("figure2_aggregated.png", p2, width = 8, height = 6)

# Now you have two ggplot objects:
# - p1 => multi-line (split by numPSVs)
# - p2 => single line for Precision, single for Recall

# Reshape the data: separate out "Precision" and "Recall" metrics with their mean and sd
df_long <- combined_filtering_results %>%
  pivot_longer(
    cols = c(Precision_mean, Recall_mean, Precision_sd, Recall_sd),
    names_to = c("metric", ".value"),
    names_pattern = "(.*)_(.*)"
  )

# Now df_long will have columns like: k_col, k_col_num, Params, metric, mean, and sd

# Plot with ggplot2: x-axis is k_col_num, and we facet by Params.
ggplot(df_long, aes(x = k_col_num, y = mean, color = metric)) +
  geom_line() +
  geom_point() +
  scale_color_manual(
    name   = "metric",
    values = c("Precision" = "#6B00B6", "Recall" = "#B3A300")
  ) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
  facet_wrap(~ Params) +
  labs(
    x = "k (as numeric)",
    y = "Mean Value",
    color = "Metric",
    title = "Precision & Recall across different k values - comparing filtering thresholds"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))