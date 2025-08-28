library(dplyr)

df_filt %>%
  group_by(k_col) %>%
  summarise(
    n_observations = n(),
    n_unique_values = n_distinct(Precision),
    .groups = "drop"
  )

# df_filt_clean <- df_filt %>%
#   group_by(k_col) %>%
#   filter(n() >= 2) %>%          # keep only groups with 2+ obs
#   filter(n_distinct(Precision) > 1)  # optionally also ensure >1 unique value
# ungroup()

pwc <- pairwise_wilcox_test(
  data = df_filt_clean,
  formula = Precision ~ k_col,
  p.adjust.method = "BH",
  paired = FALSE
)
pwc



# Unique factor levels
levels_k <- unique(df_filt$k_col)

# All pairwise combinations
pairs <- combn(levels_k, 2, simplify = FALSE)

results_list <- lapply(pairs, function(pair) {
  sub_df <-  %>% 
    filter(k_col %in% pair) %>%
    mutate(k_col = factor(k_col, levels = pair))
  
  # Attempt an unpaired Wilcoxon test
  wilcox_res <- try(wilcox.test(Precision ~ k_col, data = sub_df))
  
  # If there's an error, store it
  if (inherits(wilcox_res, "try-error")) {
    list(pair = pair, error = wilcox_res)
  } else {
    list(pair = pair, p_value = wilcox_res$p.value)
  }
})

results_list

results_list_k2 <- results_list[
  sapply(results_list, function(x) {
    "pair" %in% names(x) && "k_2" %in% x$pair
  })
]

results_list_k2


##
## Compare between confidence scores ####
##

#df_filt_0conf <- df_filt
#df_filt_0.5conf <- df_filt
#df_filt_0.25conf <- df_filt

df_filt_0conf$Conf <- "0conf"
df_filt_0.25conf$Conf <- "0.25conf"
df_filt_0.5conf$Conf <- "0.5conf"

combined_filtering_results <- rbind(df_filt_0conf,df_filt_0.25conf,df_filt_0.5conf)


combined_filtering_results <- combined_filtering_results %>%
  filter(k_col == "k_9")

wilcox.test(Precision ~ Conf, data = combined_filtering_results,)
