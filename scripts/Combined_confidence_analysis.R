library(dplyr)
library(rstatix) 

combined_filtering_results <- read.csv("combined_confidence_results.txt")

# 1) Filter to keep k_col == "k_9"
df_k9 <- combined_filtering_results %>%
  filter(k_col == "k_9")

# 2) Quick descriptive stats: mean, sd, etc., grouped by Conf
df_k9_summary <- df_k9 %>%
  group_by(Conf) %>%
  summarise(
    n = n(),
    mean_precision = mean(Precision, na.rm = TRUE),
    sd_precision   = sd(Precision,   na.rm = TRUE),
    mean_Recall = mean(Recall, na.rm = TRUE),
    sd_Recall   = sd(Recall,   na.rm = TRUE),
    .groups        = "drop"
  )

df_k9_summary

0.144
0.072 halfway
0.0288 per conf

# 3) Statistical comparison across Conf groups.
#    Example: a nonparametric test, comparing each pair of Conf levels.
res_wilcox <- df_k9 %>%
  wilcox_test(Precision ~ Conf) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

res_wilcox

# 4) Optional: a boxplot to visualize Precision by Conf
library(ggplot2)
ggplot(df_k9, aes(x = Conf, y = Precision)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Precision by Conf (k_9 only)")
