# Sequencing results ####

result_summary <- mh_TE_map %>%
  group_by(as.factor(Orig_Mh_level)) %>%
  summarise(
    num_samples = n(),
    
    mean_raw_reads = mean(raw_reads, na.rm = TRUE),
    median_raw_reads = median(raw_reads, na.rm = TRUE),
    min_raw_reads = min(raw_reads, na.rm = TRUE),
    max_raw_reads = max(raw_reads, na.rm = TRUE),
    
    mean_nonhost_reads = mean(nonhost_reads, na.rm = TRUE),
    median_nonhost_reads = median(nonhost_reads, na.rm = TRUE),
    min_nonhost_reads = min(nonhost_reads, na.rm = TRUE),
    max_nonhost_reads = max(nonhost_reads, na.rm = TRUE),
    
    mean_passed_screening = mean(passed_screening, na.rm = TRUE),
    median_passed_screening = median(passed_screening, na.rm = TRUE),
    min_passed_screening = min(passed_screening, na.rm = TRUE),
    max_passed_screening = max(passed_screening, na.rm = TRUE),
    
    mean_PSV_counts = mean(PSV_counts, na.rm = TRUE),
    median_PSV_counts = median(PSV_counts, na.rm = TRUE),
    min_PSV_counts = min(PSV_counts, na.rm = TRUE),
    max_PSV_counts = max(PSV_counts, na.rm = TRUE)
  )
result_summary


## Raw reads ####
raw_reads_mhpp_Mh_level <- pairwise.wilcox.test(mh_TE_map$raw_reads, mh_TE_map$Mh_level, p.adjust.method = "BH")
raw_reads_mhpp_Mh_level$p.value # p-values for Observed PSVs

boxplot(mh_TE_map$raw_reads ~ mh_TE_map$Mh_level)
ggplot(mh_TE_map, aes(x = as.factor(Mh_level), y = raw_reads)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.5)


## Non-host reads ####
nonhost_reads_mhpp_Mh_level <- pairwise.wilcox.test(mh_TE_map$nonhost_reads, mh_TE_map$Mh_level, p.adjust.method = "BH")
nonhost_reads_mhpp_Mh_level$p.value # p-values for Observed PSVs

boxplot(mh_TE_map$nonhost_reads ~ mh_TE_map$Mh_level)
ggplot(mh_TE_map, aes(x = as.factor(Mh_level), y = nonhost_reads)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.5)


## Screening reads, including PSV ####
passed_screening_mhpp_Mh_level <- pairwise.wilcox.test(mh_TE_map$passed_screening, mh_TE_map$Mh_level, p.adjust.method = "BH")
passed_screening_mhpp_Mh_level$p.value # p-values for Observed PSVs

boxplot(mh_TE_map$passed_screening ~ mh_TE_map$Mh_level)
ggplot(mh_TE_map, aes(x = as.factor(Mh_level), y = passed_screening)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.5)


#write.csv(sample_sums(mhpp.ps),"BRD_confirmatory_hits.csv")

## Screening reads, including PSV ####
passed_screening_mhpp_Concentration <- pairwise.wilcox.test(mh_TE_map$passed_screening, mh_TE_map$Concentration, p.adjust.method = "BH")
passed_screening_mhpp_Concentration$p.value # p-values for Observed PSVs


## PSV counts ####

sum(mh_TE_map$PSV_counts)
mean(mh_TE_map$PSV_counts)
min(mh_TE_map$PSV_counts)
max(mh_TE_map$PSV_counts)

PSV_counts_mhpp_Mh_level <- pairwise.wilcox.test(mh_TE_map$PSV_counts, mh_TE_map$Mh_level, p.adjust.method = "BH")
PSV_counts_mhpp_Mh_level$p.value # p-values for Observed PSVs

boxplot(mh_TE_map$PSV_counts ~ mh_TE_map$Mh_level)
ggplot(mh_TE_map, aes(x = as.factor(Mh_level), y = PSV_counts)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.5)

# PSV results #### 

PSV_rel_abund <- transform_sample_counts(mhpp_PSV.ps, function(x) {x/sum(x)}*100)

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

