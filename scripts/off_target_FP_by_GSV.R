# Goal is to view raw results from benchmarking and determine trends in FP rates

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(scales)

#
## All conf ####
#
df_allconf <- read.table(
  "Benchmarking_results/results_offtarget_byConf_raw_summarized.txt",  # <-- replace with your actual path
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

df_allconf <- df_allconf %>%
  mutate(
    Conf = sub(".*(conf[^_/]*)$", "\\1", sub("\\.k_\\d+$", "", file)),
    ANI_k = str_extract(file, "(?<=\\.)[^.]+$")
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
  filter(read_type == "combined") %>%
  filter(ANI_k == "k_9" ) %>%
  mutate(
    Conf = factor(
      Conf,
      levels = c("conf0","conf0.1","conf0.2","conf0.3","conf0.4","conf0.5")
    )
  )


## ---------------------------------------------------------------
## 1 · tidy the PSV dictionary column
## ---------------------------------------------------------------
df_long <- df_allconf_merged %>% 
  select(file, Conf, filtered_psv_groups,
         num_aligned, total_input_reads) %>% 
  
  ## pull every  "<digits> : <number>"  pair
  mutate(pair_vec = str_extract_all(
    filtered_psv_groups,
    "\\d+\\s*:\\s*[0-9.]+(?:[eE][+-]?[0-9]+)?"
  )) %>% 
  unnest(pair_vec) %>% 
  separate(pair_vec, into = c("psv","abs_count"), sep="\\s*:\\s*") %>% 
  mutate(
    psv        = as.integer(psv),
    abs_count  = as.numeric(abs_count),
    FP_by_psv  = abs_count / total_input_reads      # ⟵ no extra multiply
  )

# df_long <- df_allconf_merged %>%                       # only combined rows
#   select(file, Conf, filtered_psv_groups,
#          num_aligned, total_input_reads) %>% 
#   
#   mutate(pair_list = str_match_all(
#     filtered_psv_groups,
#     "(\\d+)\\s*:\\s*([0-9.]+(?:[eE][+-]?[0-9]+)?)")) %>% 
#   unnest_wider(pair_list, names_sep = "_") %>%         # 1 col per match piece
#   pivot_longer(starts_with("pair_list_"),               # long format
#                names_to = ".discard",
#                values_to = "kv_pair") %>% 
#   filter(!is.na(kv_pair)) %>% 
#   separate(kv_pair, into = c("psv","rel_abund"), sep=":", convert = TRUE) %>% 
#   mutate(
#     psv_counts = rel_abund * num_aligned,
#     FP_by_psv  = psv_counts / total_input_reads
#   )
## ---------------------------------------------------------------
## 2 · box-plot   FP_by_psv  vs Conf
## ---------------------------------------------------------------

quantile(df_long$FP_by_psv)


## ──────────────────────────────────────────────
## 1 · summary table (mean FP per Conf × PSV)
## ──────────────────────────────────────────────
conf_psv_summary <- df_long %>%             # df_long already contains one row per PSV
  group_by(Conf, psv) %>%                   # → group by BOTH factors
  summarise(
    mean_FP = mean(FP_by_psv, na.rm = TRUE),
    max_FP = max(FP_by_psv, na.rm = TRUE),
    median_FP = median(FP_by_psv, na.rm = TRUE),
    .groups = "drop"
  )

conf_psv_summary <- df_long %>% 
  group_by(Conf, psv) %>%                               # <- keep both factors
  summarise(                                            # one row per Conf × PSV
    mean_FP    = mean(FP_by_psv, na.rm = TRUE),
    max_FP     = max(FP_by_psv,  na.rm = TRUE),
    median_FP  = median(FP_by_psv, na.rm = TRUE),
    
    ## --- quantiles ---------------------------------------------------
    q50_FP     = quantile(FP_by_psv, 0.50, na.rm = TRUE),   # == median_FP
    q75_FP     = quantile(FP_by_psv, 0.75, na.rm = TRUE),
    q95_FP     = quantile(FP_by_psv, 0.95, na.rm = TRUE),
    q99_FP     = quantile(FP_by_psv, 0.99, na.rm = TRUE),
    q99.5_FP     = quantile(FP_by_psv, 0.99, na.rm = TRUE),
    q99.9_FP     = quantile(FP_by_psv, 0.99, na.rm = TRUE),
    .groups    = "drop"
  )

conf_psv_summary %>% 
  filter(Conf == "conf0.1" )

## ──────────────────────────────────────────────
## 2 · plotting
## ──────────────────────────────────────────────



ggplot(df_long, aes(x = Conf, y = FP_by_psv)) +
  ## neutral background distribution
  geom_boxplot(width = .65, fill = "grey88", colour = "grey40",
               outlier.shape = NA) +
  
  ## every observation, PSV–coloured
  geom_jitter(aes(fill = factor(psv), colour = factor(psv)),
              width = .25, height = 0,
              shape = 21, size = 1.8, alpha = .6, stroke = .25) +
  
  ## mean points per Conf×PSV
  geom_point(data = conf_psv_summary,
             aes(x = Conf, y = mean_FP, fill = factor(psv)),
             shape = 23, colour = "black",
             size  = 3, stroke = .3,
             position = dodge) +
  
  ## mean labels
  geom_text(data = conf_psv_summary,
            aes(x = Conf, y = mean_FP,
                label = format(mean_FP, digits = 2, scientific = TRUE),
                colour = factor(psv)),
            position = dodge, vjust = -0.8, size = 2.6,
            show.legend = FALSE) +
  
  scale_fill_brewer(palette = "Set2", name = "PSV") +
  scale_colour_brewer(palette = "Set2", name = "PSV") +
  scale_y_continuous(labels = scales::number_format(accuracy = .0001)) +
  labs(x = "Confidence threshold",
       y = "FP read proportion",
       title = "GSV FP read proportion by Kraken confidence") +
  theme_classic() +
  theme(legend.position = "right")



## ───────────────────────────────────────────────────────────────
## 1 · mean FP for every (Conf × PSV)  ───────────────────────────
## ───────────────────────────────────────────────────────────────
mean_df <- df_long %>%
  group_by(Conf, psv) %>%
  summarise(mean_FP = mean(FP_by_psv, na.rm = TRUE), .groups = "drop")

## ───────────────────────────────────────────────────────────────
## 2 · plot: box-plots + jittered points + mean markers / labels ─
## ───────────────────────────────────────────────────────────────
dodge <- position_dodge(width = .3)          # reuse for text & triangles

ggplot(df_long,
       aes(x = factor(psv),
           y = FP_by_psv,
           fill = factor(psv))) +
  
  ## base boxes
  geom_boxplot(outlier.shape = NA,
               alpha = .7) +
  
  ## raw points (optional – comment out if cluttered)
  geom_jitter(aes(colour = factor(psv)),
              width = .18, height = 0,
              size = 1.0, alpha = .35,
              show.legend = FALSE) +
  
  ## mean FP – red triangles
  geom_point(data = mean_df,
             aes(x = factor(psv), y = mean_FP),
             colour = "red",
             shape  = 17, size = 1.3,
             position = dodge) +
  
  ## mean labels just above the triangles
  geom_text(
    data      = mean_df,
    aes(x = factor(psv), y = mean_FP,
        label = sprintf("%.5f", mean_FP)),
    colour    = "red",
    angle     = 45,        # ← tilt 45°
    hjust     = 0,         #   align left edge of the text
    vjust     = -0.5,      #   keep it above the box-whisker
    position  = position_dodge(width = .65),
    size      = 2.5,
    show.legend = FALSE
  ) +
  
  ## facet one strip per confidence
  facet_wrap(~ Conf, nrow = 1, strip.position = "bottom") +
  
  ## cosmetics
  scale_fill_brewer(palette = "Set2", name = "PSV") +
  scale_colour_brewer(palette = "Set2", guide = "none") +
  scale_y_continuous(labels = number_format(accuracy = .0001)) +
  
  labs(x = "PSV ID",
       y = "FP read proportion",
       title = "FP proportion by PSV (faceted by Kraken confidence)") +
  
  theme_classic(base_size = 11) +
  theme(
    strip.placement  = "outside",
    strip.background = element_blank(),
    legend.position  = "right",
    
    ## thin grey borders and some space between facets
    panel.border     = element_rect(colour = "grey60", fill = NA, linewidth = 0.3),
    panel.spacing.x  = unit(1, "lines")
  )



## Facet by PSV instead

dodge <- position_dodge(width = .55)    # reuse

ggplot(df_long,
       aes(x = factor(Conf),                 # ➊ x is Conf now
           y = FP_by_psv,
           fill = factor(Conf))) +           # keep fill = Conf (optional)
  
  ## ── base boxes ──────────────────────────────────────────
  geom_boxplot(outlier.shape = NA,
               alpha = .7) +
  
  ## ── raw points ─────────────────────────────────────────
  geom_jitter(aes(colour = factor(Conf)),
              width = .18, height = 0,
              size  = 1.0, alpha = .35,
              show.legend = FALSE) +
  
  ## ── red triangles at the mean ─────────────────────────
  geom_point(data = mean_df,
             aes(x = factor(Conf),            # ➋ update here
                 y = mean_FP),
             colour = "red",
             shape  = 17, size = 1.3,
             position = dodge) +
  
  ## ── numeric labels above those triangles ──────────────
  geom_text(
    data   = mean_df,
    aes(x = factor(Conf), y = mean_FP,        # ➌ update here
        label = sprintf("%.5f", mean_FP)),
    colour = "red",
    angle  = 45,  hjust = 0,  vjust = -0.5,
    position = dodge,
    size   = 2.5,
    show.legend = FALSE
  ) +
  
  ## ── facet by PSV instead of Conf ──────────────────────
  facet_wrap(~ psv, nrow = 1, strip.position = "bottom") +   # ➍
  
  ## ── cosmetics (unchanged) ─────────────────────────────
  scale_fill_brewer(palette = "Set2", name = "Conf") +
  scale_colour_brewer(palette = "Set2", guide = "none") +
  scale_y_continuous(labels = scales::number_format(accuracy = .0001)) +
  
  labs(x = "Kraken confidence",             # axis labels flipped
       y = "FP read proportion",
       title = "FP proportion by confidence (faceted by PSV)",
       fill = "Conf") +
  
  theme_classic(base_size = 11) +
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        legend.position  = "right",
        panel.border     = element_rect(colour = "grey60", fill = NA, linewidth = 0.3),
        panel.spacing.x  = unit(1, "lines"),
        axis.text.x = element_text(angle=20)
        )



## ── 3.  how many PSVs / samples are > mean?  ──────────────────────
above_mean_tbl <- df_long %>% 
  left_join(conf_summary, by = "Conf") %>%          # bring in the mean per Conf
  filter(FP_by_psv > mean_FP) %>%                   # keep only rows above the mean
  group_by(Conf) %>% 
  summarise(
    psvs_above_mean     = n(),                      # number of PSV entries
    samples_with_above  = n_distinct(file),         # how many samples they appear in
    .groups             = "drop"
  )

print(above_mean_tbl)
## ---------------------------------------------------------------
## 3 · which PSV has the highest total count?
## ---------------------------------------------------------------
psv_summary <- df_long %>% 
  group_by(psv, Conf) %>% 
  summarise(
    total_psv_counts = sum(abs_count, na.rm = TRUE),
    mean_psv_counts  = mean(abs_count, na.rm = TRUE),   # ← new
    n_samples        = sum(abs_count > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  arrange(desc(total_psv_counts))

psv_summary

## PSV mean counts by confidence ####
ggplot(psv_summary,
       aes(x = Conf,
           y = mean_psv_counts,
           colour = factor(psv),
           group  = psv)) +                 #   ← tell ggplot which points to connect
  geom_line(position = position_dodge(width = .45),  # connect them
            linewidth = .6, alpha = .8) +
  geom_point(position = position_dodge(width = .45), # the dots
             size = 2.3, alpha = .8) +
  scale_y_continuous(labels = comma,
                     trans  = "log10") +
  labs(x      = "Confidence threshold",
       y      = "Mean PSV counts / sample",
       colour = "PSV",
       title  = "Mean PSV counts by Kraken confidence") +
  theme_classic()


psv_conf_matrix <- psv_summary %>% 
  select(Conf, psv, mean_psv_counts) %>%               # keep what you want
  pivot_wider(                                         # wide format
    names_from  = psv,
    values_from = mean_psv_counts,
    values_fill = 0                                    # 0 where missing
  )

psv_conf_matrix

## Only top IDs #####


top_ids <- psv_summary %>% 
  group_by(psv) %>% 
  summarise(total = sum(total_psv_counts)) %>% 
  top_n(10, total) %>%              # keep top 10
  pull(psv)

psv_summary %>% 
  filter(psv %in% top_ids) %>% 
  ggplot( aes(x = Conf,
             y = mean_psv_counts,
             colour = factor(psv))) +          # colour by PSV ID
  geom_point(size = 2.3, alpha = .8,
             position = position_dodge(width = .45)) +
  scale_y_continuous(labels = comma,
                     trans  = "log10") +     # ← drop if you prefer linear scale
  labs(x        = "Confidence threshold",
       y        = "Mean PSV counts / sample",
       colour   = "PSV",
       title    = "Mean PSV counts by Kraken confidence") +
  theme_classic()

psv_conf_matrix <- psv_summary %>% 
  select(Conf, psv, mean_psv_counts) %>%               # keep what you want
  pivot_wider(                                         # wide format
    names_from  = psv,
    values_from = mean_psv_counts,
    values_fill = 0                                    # 0 where missing
  )
