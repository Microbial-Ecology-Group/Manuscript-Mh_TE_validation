library(phyloseq)
library(ggplot2)

# Ensure OTU table is properly formatted
otu_mat <- otu_table(mh_TE_PSV_transformed.ps)

# Check if taxa are rows; transpose if necessary
if (!taxa_are_rows(mh_TE_PSV_transformed.ps)) {
  otu_mat <- t(otu_mat)
}

# Calculate prevalence for each taxon
taxa_prevalence <- apply(otu_mat, 1, function(x) sum(x > 0))

# Convert to a data frame
prevalence_df <- data.frame(
  Taxa = names(taxa_prevalence),  # Use names for taxa names
  Prevalence = taxa_prevalence    # Prevalence values
)

# Check the resulting data frame
head(prevalence_df)

ggplot(prevalence_df, aes(x = Prevalence)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(
    title = "Taxa Prevalence Across All Samples",
    x = "Number of Samples",
    y = "Number of Taxa"
  )
