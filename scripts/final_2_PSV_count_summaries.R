mh_TE_PSV_transformed.ps

mean(mh_TE_map$raw_paired_reads)
min(mh_TE_map$raw_paired_reads)
max(mh_TE_map$raw_paired_reads)

mean(mh_TE_map$Extracted_reads)
min(mh_TE_map$Extracted_reads)
max(mh_TE_map$Extracted_reads)


# Sequencing results by Mh_level ####
result_summary <- mh_TE_map %>%
  dplyr::group_by(as.factor(Mh_level)) %>%
  dplyr::summarise(
    num_samples = n(),
    mean_Raw_reads_paired_added   = mean(Raw_reads_paired_added   , na.rm = TRUE),
    median_Raw_reads_paired_added   = median(Raw_reads_paired_added   , na.rm = TRUE),
    min_Raw_reads_paired_added   = min(Raw_reads_paired_added   , na.rm = TRUE),
    max_Raw_reads_paired_added   = max(Raw_reads_paired_added   , na.rm = TRUE),
    
    mean_Nonhost_reads = mean( Nonhost_reads , na.rm = TRUE),
    median_Nonhost_reads = median( Nonhost_reads , na.rm = TRUE),
    min_Nonhost_reads = min( Nonhost_reads , na.rm = TRUE),
    max_Nonhost_reads = max( Nonhost_reads , na.rm = TRUE),
    
    mean_Extracted_reads = mean( Extracted_reads , na.rm = TRUE),
    median_Extracted_reads = median( Extracted_reads , na.rm = TRUE),
    min_Extracted_reads = min( Extracted_reads , na.rm = TRUE),
    max_Extracted_reads = max( Extracted_reads , na.rm = TRUE)
  )
result_summary


## Raw reads ####
Raw_reads_paired_added_mhpp_Mh_level <- pairwise.wilcox.test(mh_TE_map$Raw_reads_paired_added  , mh_TE_map$Mh_level, p.adjust.method = "BH")
Raw_reads_paired_added_mhpp_Mh_level$p.value # p-values for raw read comparisons


boxplot(mh_TE_map$Raw_reads_paired_added   ~ mh_TE_map$Mh_level)
ggplot(mh_TE_map, aes(x = as.factor(Mh_level), y = Raw_reads_paired_added  )) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.5)


## Non-host reads ####
Nonhost_reads_mhpp_Mh_level <- pairwise.wilcox.test(mh_TE_map$Nonhost_reads, mh_TE_map$Mh_level, p.adjust.method = "BH")
Nonhost_reads_mhpp_Mh_level$p.value # p-values for Observed PSVs

boxplot(mh_TE_map$Nonhost_reads ~ mh_TE_map$Mh_level)
ggplot(mh_TE_map, aes(x = as.factor(Mh_level), y = Nonhost_reads)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.5)


## Screening reads, including PSV ####
Kraken_nt_MH_counts_mhpp_Mh_level <- pairwise.wilcox.test(mh_TE_map$Kraken_nt_MH_counts, mh_TE_map$Mh_level, p.adjust.method = "BH")
Kraken_nt_MH_counts_mhpp_Mh_level


boxplot(mh_TE_map$Kraken_nt_MH_counts ~ mh_TE_map$Mh_level)
ggplot(mh_TE_map, aes(x = as.factor(Mh_level), y = Kraken_nt_MH_counts)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.5)


#write.csv(sample_sums(mhpp.ps),"BRD_confirmatory_hits.csv")

## PSV counts ####

sum(mh_TE_map$Extracted_reads)
mean(mh_TE_map$Extracted_reads)
mean(mh_TE_map$Percent_Mh_nonhost)



min(mh_TE_map$Extracted_reads)
max(mh_TE_map$Extracted_reads)

Extracted_reads_mhpp_Mh_level <- pairwise.wilcox.test(mh_TE_map$Extracted_reads, mh_TE_map$Mh_level, p.adjust.method = "BH")
Extracted_reads_mhpp_Mh_level$p.value # p-values for Observed PSVs

boxplot(mh_TE_map$Extracted_reads ~ mh_TE_map$Mh_level)
ggplot(mh_TE_map, aes(x = as.factor(Mh_level), y = Extracted_reads)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.5)


boxplot(mh_TE_map$Extracted_reads ~ mh_TE_map$BRD)
ggplot(mh_TE_map, aes(x = as.factor(BRD), y = Extracted_reads)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.5)

# PSV results #### 

PSV_rel_abund <- transform_sample_counts(mh_TE_PSV_transformed.ps, function(x) {x/sum(x)}*100)

PSV_rel_abund_melt <- psmelt(PSV_rel_abund)

# calculate percentage across all counts
all_sample_PSV_by_abund <- PSV_rel_abund_melt %>%
  dplyr::group_by(OTU) %>%
  dplyr::summarize(median_PSV = median(Abundance)) %>%
  arrange(-median_PSV)

all_sample_PSV_by_abund

# 
# # Example for doing at other taxa level, this doesn't work for PSVs
# ra_class <- tax_glom(PSV_rel_abund, taxrank = "Class")
# ra_class 
# 
# ra_class_melt <- psmelt(ra_class)
# 
# # calculate percentage across all counts
# all_sample_PSV_class_by_abund <- ra_class_melt %>%
#   dplyr::group_by(class) %>%
#   dplyr::summarize(median_PSV = median(Abundance)) %>%
#   arrange(-median_PSV)
# 
# all_sample_PSV_class_by_abund
# Group by Mh_level and apply the Wilcoxon test for each group

library(purrr)
Extracted_reads_mhpp_BRD_by_Mh_level <- mh_TE_map %>%
  group_by(Mh_level) %>%
  group_split() %>%
  map(~ pairwise.wilcox.test(
    .x$Extracted_reads,  # Use Extracted_reads for this group
    .x$BRD,                     # Test by BRD status
    p.adjust.method = "BH"       # Adjust p-values by Benjamini-Hochberg method
  ))

# Optional: Give the results meaningful names (each unique Mh_level)
names(Extracted_reads_mhpp_BRD_by_Mh_level) <- unique(mh_TE_map$Mh_level)

# Access the results for each group:
# Example: results for the first Mh_level
Extracted_reads_mhpp_BRD_by_Mh_level$zero
Extracted_reads_mhpp_BRD_by_Mh_level$low
Extracted_reads_mhpp_BRD_by_Mh_level$medium
Extracted_reads_mhpp_BRD_by_Mh_level$high
Extracted_reads_mhpp_BRD_by_Mh_level$highest

Extracted_reads_mhpp_Concentration_by_Mh_level <- mh_TE_map %>%
  group_by(Mh_level) %>%
  group_split() %>%
  map(~ pairwise.wilcox.test(
    .x$Extracted_reads,  # Use Extracted_reads for this group
    .x$Concentration,                     # Test by Concentration status
    p.adjust.method = "BH"       # Adjust p-values by Benjamini-Hochberg method
  ))

# Optional: Give the results meaningful names (each unique Mh_level)
names(Extracted_reads_mhpp_Concentration_by_Mh_level) <- unique(mh_TE_map$Mh_level)

# Access the results for each group:
# Example: results for the first Mh_level
Extracted_reads_mhpp_Concentration_by_Mh_level$zero
Extracted_reads_mhpp_Concentration_by_Mh_level$low
Extracted_reads_mhpp_Concentration_by_Mh_level$medium
Extracted_reads_mhpp_Concentration_by_Mh_level$high
Extracted_reads_mhpp_Concentration_by_Mh_level$highest


### average classified PSVs ####
# Calculate percent_PSV_classified
mh_TE_map$percent_PSV_classified <- (mh_TE_map$Extracted_reads / mh_TE_map$Raw_reads_paired_added  ) * 100
mean(mh_TE_map$percent_PSV_classified)

