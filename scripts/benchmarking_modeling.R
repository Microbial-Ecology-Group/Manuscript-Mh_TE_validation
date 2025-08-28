
### Linear Mixed-Effects Model (Parametric, flexible) ####

## Install if not present: install.packages("lme4"); install.packages("lmerTest")
library(lme4)
library(lmerTest)  # adds p-values to lme4 models
library(emmeans)   # for post-hoc comparisons

# 1) Fit a mixed model, treating k_col as a fixed effect, iter_num as a random effect
mod <- lmer(Precision ~ k_col + (1 | iter_num), data = df_filt)

# 2) ANOVA table (tests if there's an overall k_col effect)
anova(mod)

# 3) Post-hoc pairwise comparisons among k_col levels
#    'emmeans' gives you estimated marginal means for each k_col, with pairwise contrasts
pairwise <- emmeans(mod, pairwise ~ k_col, adjust = "tukey")
pairwise

### Friedman Test (Non-Parametric Repeated-Measures) ####

# Suppose each iter_num has exactly one Precision measure for each k_col
# 1) Reshape: wide format => one row per iter_num, columns = k_col
df_wide <- df_filt %>%
  # 1) Group by iter_num, k_col (ignoring numPSVs)
  group_by(iter_num, k_col) %>%
  # 2) Average across duplicates if they differ by numPSVs
  summarise(
    Precision = mean(Precision, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  # 3) Pivot: one row per iter_num, one column per k_col
  pivot_wider(
    names_from  = k_col,
    values_from = Precision
  )
# Remove iter_num column from test input
df_mat <- as.matrix(df_wide[ , -1])

# 2) Friedman test
friedman.test(df_mat)

library(rstatix)

df_long_again <- df_wide %>%
  pivot_longer(
    cols = -iter_num,
    names_to = "k_col",
    values_to = "precision"
  )

# Friedman test with rstatix
#friedman_test(precision ~ k_col | iter_num, data = df_long_again)

# Then a pairwise Wilcoxon post-hoc with BH adjustment:
pwc <- wilcox_test(precision ~ k_col, paired = FALSE, data = df_long_again) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
pwc

# 1) Kruskal–Wallis test (non-parametric, for independent groups)
kruskal_res <- kruskal_test(Precision ~ k_col, data = df_filt)
kruskal_res



df_filt %>%
  group_by(k_col) %>%
  summarise(
    n_total  = n(),
    n_non_na = sum(!is.na(Precision)),
    unique_values = n_distinct(Precision, na.rm = TRUE)
  ) %>%
  arrange(n_non_na)

df_filt_clean <- df_filt %>%
  filter(!is.na(Precision), !is.na(k_col)) %>%
  group_by(k_col) %>%
  filter(n() >= 2) %>%          # keep only groups with at least 2 observations
  ungroup()

pwc <- pairwise_wilcox_test(
  data            = df_filt_clean,
  formula         = Precision ~ k_col,
  p.adjust.method = "BH",
  paired          = FALSE,   # unpaired if your data are genuinely independent
  exact = FALSE
)
pwc

df_filt_clean %>%
  summarise(
    any_na_precision = any(is.na(Precision)),
    any_na_k_col     = any(is.na(k_col))
  )
