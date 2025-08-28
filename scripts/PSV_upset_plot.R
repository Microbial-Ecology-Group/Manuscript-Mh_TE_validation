# PSV upset plot

library(UpSetR)
library("MicrobiotaProcess")


sample_data(mh_TE_PSV.ps)$Mh_level

sample_data(mh_TE_PSV.ps)$Mh_level <- factor(sample_data(mh_TE_PSV.ps)$Mh_level , levels = c("zero","low","medium","high","highest"))
upsetda_edit <- get_upset(mh_TE_PSV.ps, factorNames="Mh_level") ## ASV

filtered_data <- upsetda_edit[upsetda_edit$zero == 0, ]

upset(upsetda_edit, sets=c("zero","low","medium", "high","highest"),
      
      sets.bar.color = c("#ff3333", "#ff0000", "#cc0000", "#990000", "#808080")
      ,text.scale = 2,
      order.by = "freq")



taxa_names <- row.names(filtered_data)

# Now, create a logical vector indicating whether each taxa in the phyloseq object
# is in the taxa_names vector
taxa_to_keep <- taxa_names(mh_TE_PSV.ps) %in% taxa_names

# Use prune_taxa with the logical vector to subset the phyloseq object
mh_TE_PSV_filtered <- prune_taxa(taxa_to_keep, mh_TE_PSV.ps)

boxplot(sample_sums(mh_TE_PSV_filtered) ~ sample_data(mh_TE_PSV_filtered)$Mh_level)


### Splitting by abundance #####
# Calculate the abundance threshold
# Here, we're using the median abundance as an example threshold
abundance_threshold <- median(taxa_sums(mh_TE_PSV.ps))
#abundance_threshold <- 1000


# Identify low and high abundance taxa
low_abundance_taxa <- taxa_sums(mh_TE_PSV.ps) <= abundance_threshold
high_abundance_taxa <- taxa_sums(mh_TE_PSV.ps) > abundance_threshold

# Create two new phyloseq objects
mh_TE_PSV_low_abundance <- prune_taxa(low_abundance_taxa, mh_TE_PSV.ps)
mh_TE_PSV_high_abundance <- prune_taxa(high_abundance_taxa, mh_TE_PSV.ps)

# mh_TE_PSV_low_abundance now contains only the low abundance taxa
# mh_TE_PSV_high_abundance now contains only the high abundance taxa

upsetda_edit_low <- get_upset(mh_TE_PSV_low_abundance, factorNames="Mh_level") ## ASV
upset(upsetda_edit_low, sets=c("zero","low","medium", "high","highest"),
      
      sets.bar.color = c("#ff3333", "#ff0000", "#cc0000", "#990000", "#808080")
      ,text.scale = 2,
      order.by = "freq")


upsetda_edit_high <- get_upset(mh_TE_PSV_high_abundance, factorNames="Mh_level") ## ASV
upset(upsetda_edit_high, sets=c("zero","low","medium", "high","highest"),
      
      sets.bar.color = c("#ff3333", "#ff0000", "#cc0000", "#990000", "#808080")
      ,text.scale = 2,
      order.by = "freq")
