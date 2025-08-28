# Comparing confidence values for off target confidence values


library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(scales)

#
## 0conf ####
#
df_0conf <- read.table(
  "Benchmarking_results/results_offtarget_0conf_summarized.txt",  # <-- replace with your actual path
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

df_0conf <- df_0conf %>%
  filter(read_type == "merged")

df_0conf %>%
  group_by(read_type, numPSVs) %>%
  summarise(
    num_aligned_mean = mean(num_aligned, na.rm = TRUE),
    num_aligned_sd   = sd(num_aligned, na.rm = TRUE),
    .groups = "drop"
  )

df_0conf$Conf <- "0conf"

#
## 0.05conf ####
#
df_0.05conf <- read.table(
  "Benchmarking_results/results_offtarget_0.05conf_summarized.txt",  # <-- replace with your actual path
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

df_0.05conf <- df_0.05conf %>%
  filter(read_type == "merged")

df_0.05conf %>%
  group_by(numPSVs) %>%
  summarise(
    num_aligned_mean = mean(num_aligned, na.rm = TRUE),
    num_aligned_sd   = sd(num_aligned, na.rm = TRUE),
    .groups = "drop"
  )

df_0.05conf$Conf <- "0.05conf"


#
## 0.1conf ####
#
df_0.1conf <- read.table(
  "Benchmarking_results/results_offtarget_0.1conf_summarized.txt",  # <-- replace with your actual path
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

df_0.1conf <- df_0.1conf %>%
  filter(read_type == "merged")

df_0.1conf %>%
  group_by(numPSVs) %>%
  summarise(
    num_aligned_mean = mean(num_aligned, na.rm = TRUE),
    num_aligned_sd   = sd(num_aligned, na.rm = TRUE),
    .groups = "drop"
  )

df_0.1conf$Conf <- "0.1conf"

#
## 0.15conf ####
#
df_0.15conf <- read.table(
  "Benchmarking_results/results_offtarget_0.15conf_summarized.txt",  # <-- replace with your actual path
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

df_0.15conf <- df_0.15conf %>%
  filter(read_type == "merged")

df_0.15conf %>%
  group_by(numPSVs) %>%
  summarise(
    num_aligned_mean = mean(num_aligned, na.rm = TRUE),
    num_aligned_sd   = sd(num_aligned, na.rm = TRUE),
    .groups = "drop"
  )

df_0.15conf$Conf <- "0.15conf"


#
## 0.2conf ####
#
df_0.2conf <- read.table(
  "Benchmarking_results/results_offtarget_0.2conf_summarized.txt",  # <-- replace with your actual path
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

df_0.2conf <- df_0.2conf %>%
  filter(read_type == "merged")

df_0.2conf %>%
  group_by(numPSVs) %>%
  summarise(
    num_aligned_mean = mean(num_aligned, na.rm = TRUE),
    num_aligned_sd   = sd(num_aligned, na.rm = TRUE),
    .groups = "drop"
  )

df_0.2conf$Conf <- "0.2conf"


#
## 0.3conf ####
#
df_0.3conf <- read.table(
  "Benchmarking_results/results_offtarget_0.3conf_summarized.txt",  # <-- replace with your actual path
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

df_0.3conf <- df_0.3conf %>%
  filter(read_type == "merged")

df_0.3conf %>%
  group_by(numPSVs) %>%
  summarise(
    num_aligned_mean = mean(num_aligned, na.rm = TRUE),
    num_aligned_sd   = sd(num_aligned, na.rm = TRUE),
    .groups = "drop"
  )

df_0.3conf$Conf <- "0.3conf"


#
## 0.4conf ####
#
df_0.4conf <- read.table(
  "Benchmarking_results/results_offtarget_0.4conf_summarized.txt",  # <-- replace with your actual path
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

df_0.4conf <- df_0.4conf %>%
  filter(read_type == "merged")

df_0.4conf %>%
  group_by(numPSVs) %>%
  summarise(
    num_aligned_mean = mean(num_aligned, na.rm = TRUE),
    num_aligned_sd   = sd(num_aligned, na.rm = TRUE),
    .groups = "drop"
  )

df_0.4conf$Conf <- "0.4conf"


#
## 0.5conf ####
#
df_0.5conf <- read.table(
  "Benchmarking_results/results_offtarget_0.5conf_summarized.txt",  # <-- replace with your actual path
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

df_0.5conf <- df_0.5conf %>%
  filter(read_type == "merged")

df_0.5conf %>%
  group_by(numPSVs) %>%
  summarise(
    num_aligned_mean = mean(num_aligned, na.rm = TRUE),
    num_aligned_sd   = sd(num_aligned, na.rm = TRUE),
    .groups = "drop"
  )

df_0.5conf$Conf <- "0.5conf"


combined_offtarget_filtering_results <- rbind(df_0conf,df_0.05conf,df_0.1conf,df_0.15conf,df_0.2conf,df_0.3conf,df_0.4conf,df_0.5conf)



combined_offtarget_filtering_results %>%
  group_by(Conf, numPSVs) %>%
  summarise(
    mean_aligned_reads = mean(num_aligned, na.rm = TRUE),
    median_aligned_reads = median(num_aligned),
    sd_aligned_reads   = sd(num_aligned, na.rm = TRUE),
    min_aligned_reads = min(num_aligned),
    max_aligned_reads = max(num_aligned),
    number_of_samples = n(),
    .groups = "drop"
  )


#
##
# Make boxplot ####
##
#

# 1) Make Conf an ordered factor
combined_offtarget_filtering_results <- combined_offtarget_filtering_results %>%
  mutate(
    Conf = factor(
      Conf,
      levels = c("0conf","0.05conf","0.1conf","0.15conf","0.2conf","0.3conf","0.4conf","0.5conf")
    )
  )

# 2) Calculate how many samples per Conf for labeling
counts_df <- combined_offtarget_filtering_results %>%
  group_by(Conf) %>%
  summarise(
    n_samps = n(),
    # We also grab the max(num_aligned) so we can place text just above each box
    max_aligned = max(num_aligned, na.rm = TRUE)
  ) %>%
  # For display, clamp the label's y-position to 6500
  mutate(
    label_ypos = if_else(max_aligned > 1000, 900, max_aligned)
  )

# 3) Boxplot with label text
ggplot(combined_offtarget_filtering_results, aes(x = Conf, y = num_aligned)) +
  geom_boxplot() +
  # Add text for sample size, placed slightly above each box’s max
  geom_text(
    data = counts_df,
    aes(x = Conf, y = max_aligned, label =  paste0("n=", n_samps)),
    vjust = -0.5  # move text a bit above the max
  ) +
  labs(
    x = "Conf",
    y = "num_aligned",
    title = "Boxplot of num_aligned by Conf - 100 iterations"
  ) +
  theme_minimal()


 ggplot(combined_offtarget_filtering_results, aes(x = Conf, y = num_aligned)) +
  geom_boxplot() +
  geom_text(
    data = counts_df,
    aes(x = Conf, y = label_ypos, label = paste0("n=", n_samps)),
    vjust = -0.5
  ) +
  labs(
    x = "Conf",
    y = "num_aligned",
    title = "Boxplot of num_aligned by Conf (zoomed 0–1000)"
  ) +
  theme_minimal() +
  # "Zoom" the Y axis from 0 to 6500
  coord_cartesian(ylim = c(0, 1000))

 