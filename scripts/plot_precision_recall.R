library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(scales)

# --------------------------------------------------
# 1) Read the summarized file
# --------------------------------------------------
df <- read.table(
  "Benchmarking_results/results_100i_0.5conf_summarized.txt",  # <-- replace with your actual path
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)


# --------------------------------------------------
# 2) (Optional) Filter out rows where k_col == "off_target"
# --------------------------------------------------
df_filt <- df %>%
  filter(read_type == "merged")

df_filt <- df_filt %>%
  filter(k_col != "off_target")

# Also remove cases where numPSVs is greater than k_col
nrow(df_filt) # 14442

df_filt <- df_filt %>%
  mutate(
    nPSVs = as.numeric(str_extract(numPSVs, "\\d+")),
    nK    = as.numeric(str_extract(k_col, "\\d+"))
  ) %>%
  filter(nK >= nPSVs)  

nrow(df_filt) #12842

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

str(df_filt)




## 1) Group by k_col: compute the average Precision & Recall
df_pr <- df_filt %>%
  group_by(k_col, nK) %>%  # assuming `nK` is directly available in df_filt
  summarise(
    mean_precision = mean(Precision, na.rm = TRUE),
    mean_recall    = mean(Recall,    na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Compute F1 if desired
  mutate(mean_f1 = 2 * mean_precision * mean_recall / (mean_precision + mean_recall)) %>%
  # Sort in ascending order of Recall for a left-to-right line
  arrange(mean_recall)

## 2) Plot: 
##    - Precision vs. Recall
##    - Color points on a gradient from light to dark blue by numeric nK
##    - Label each point with k_col
ggplot(df_pr, aes(x = mean_recall, y = mean_precision, color = nK)) +
  # optional line to connect thresholds in ascending recall
  geom_line(alpha = 0.3) +
  geom_point(size = 3) +
  geom_text(
    aes(label = k_col),
    vjust = 2, alpha = 0.7, check_overlap = TRUE
  ) +
  scale_color_gradient(low = "blue", high = "lightblue") +
  labs(
    title = "Precision–Recall by k_col Threshold",
    x = "Mean Recall",
    y = "Mean Precision",
    color = "nK"
  ) +
  theme_minimal()

View(df_pr)
