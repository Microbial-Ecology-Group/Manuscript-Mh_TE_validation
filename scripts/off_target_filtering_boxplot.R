library(dplyr)
library(stringr)
library(purrr)  # for map() and related list operations

df_allconf_merged_parsed <- df_allconf_merged %>%
  # 1) Remove the curly braces { } and split by comma
  mutate(
    psv_groups_clean = str_remove_all(filtered_psv_groups, "[{}]"),
    psv_list = str_split(psv_groups_clean, ",")
  ) %>%
  # 2) For each row, parse the read counts, then compute mean, sd, min, max
  rowwise() %>%
  mutate(
    # Convert each "PSV: readcount" string into numeric read counts
    read_values = list(
      if (length(psv_list) == 0) {
        # If there's no data, return an empty numeric vector
        numeric(0)
      } else {
        sapply(psv_list, function(item) {
          # Trim whitespace
          item <- str_trim(item)
          # Split on ":" => c("PSV_ID", "reads")
          parts <- str_split(item, ":", n = 2)[[1]]
          # Convert the read count to numeric
          as.numeric(str_trim(parts[2]))
        })
      }
    ),
    # Now calculate your statistics from read_values
    mean_reads_per_PSV = mean(read_values, na.rm = TRUE),
    sd_reads_per_PSV   = sd(read_values, na.rm = TRUE),
    min_reads_per_PSV  = if (length(read_values) > 0) min(read_values) else NA_real_,
    max_reads_per_PSV  = if (length(read_values) > 0) max(read_values) else NA_real_
  ) %>%
  ungroup()

# Inspect the new columns
df_allconf_merged_parsed %>%
  select(filtered_psv_groups, mean_reads_per_PSV, sd_reads_per_PSV, min_reads_per_PSV, max_reads_per_PSV) %>%
  head()

##
### Summary view
##

df_allconf_merged_parsed %>%
  group_by(Conf) %>%
  summarise(
    mean_aligned_reads = mean(num_aligned, na.rm = TRUE),
    #median_aligned_reads = median(num_aligned),
    sd_aligned_reads   = sd(num_aligned, na.rm = TRUE),
    min_aligned_reads = min(num_aligned),
    max_aligned_reads = max(num_aligned),
    number_of_samples = n(),
    .groups = "drop"
  )


df_allconf_merged_parsed %>%
  group_by(Conf) %>%
  summarise(
    mean_reads_per_PSV = mean(mean_reads_per_PSV),
    mean_aligned_reads = mean(num_aligned, na.rm = TRUE),
    sd_reads_per_PSV = median(sd_reads_per_PSV, na.rm = TRUE),
    min_reads_per_PSV = min(min_reads_per_PSV, na.rm = TRUE),
    max_reads_per_PSV = max(max_reads_per_PSV, na.rm = TRUE),
    number_of_samples = n(),
    .groups = "drop"
  )

2445/2500 * 100
1490/2500 * 100
855/2500 * 100
455/2500 * 100
340/2500 * 100
280/2500 * 100




# With sub count rel abund
2325/2500 * 100
1770/2500 * 100
1175/2500 * 100
785/2500 * 100
625/2500 * 100
445/2500 * 100
350/2500 * 100
280/2500 * 100

10/2500 * 100


1000 / 250000 * 100

results_filtered_summary <- df_allconf_merged_parsed %>%
  # Subtract 1000 from num_aligned
  mutate(num_aligned_adj = num_aligned - 1000) %>%
  # Keep only rows where adjusted num_aligned > 0
  filter(num_aligned_adj > 0) %>%
  # Now group and summarize using the adjusted alignment counts
  group_by(Conf, numPSVs) %>%
  summarise(
    mean_reads_per_PSV  = mean(mean_reads_per_PSV, na.rm = TRUE),
    mean_aligned_reads  = mean(num_aligned_adj, na.rm = TRUE),
    #sd_reads_per_PSV    = median(sd_reads_per_PSV, na.rm = TRUE),
    min_reads_per_PSV   = min(min_reads_per_PSV, na.rm = TRUE),
    max_reads_per_PSV   = max(max_reads_per_PSV, na.rm = TRUE),
    number_of_samples   = n(),
    .groups = "drop"
  )

results_filtered_summary


S##
# Make boxplot ####
##
#

# 1) Make Conf an ordered factor
df_allconf_merged <- df_allconf_merged_parsed %>%
  mutate(
    Conf = factor(
      Conf,
      levels = c("0conf","0.05conf","0.1conf","0.15conf","0.2conf","0.3conf","0.4conf","0.5conf")
    )
  )

# 2) Calculate how many samples per Conf for labeling
counts_df <- df_allconf_merged %>%
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
ggplot(df_allconf_merged, aes(x = Conf, y = num_aligned)) +
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


# 3) Boxplot with label text, mean_reads_per_PSV
ggplot(df_allconf_merged, aes(x = Conf, y = mean_reads_per_PSV)) +
  geom_boxplot() +
  # Add text for sample size, placed slightly above each box’s max
  geom_text(
    data = counts_df,
    aes(x = Conf, y = max_aligned, label =  paste0("n=", n_samps)),
    vjust = -0.5  # move text a bit above the max
  ) +
  labs(
    x = "Conf",
    y = "mean_reads_per_PSV",
    title = "Boxplot of mean reads per GSV by Conf - 100 iterations"
  ) +
  theme_minimal()



ggplot(df_allconf_merged, aes(x = Conf, y = num_aligned)) +
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

