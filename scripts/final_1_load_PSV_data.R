#### Mh++ load msweep data ####

# Read in metadata
mh_TE_map <- read.csv("metadata_swabcomp.csv", header = T, row.names = 1)
mh_TE_map$X16S_qPCR_mean_copy
# Read in count matrix
mh_TE_otu.df <- read.table("BRDnoBRD_mSWEEP_results/mSweep_results_sub_count_rel_abund_byGSV_c0_99q_count_matrix.tsv", header = T, row.names = 1, sep = "\t")
# Remove scientific notation by formatting the matrix with sprintf to control the number of decimal places
mh_TE_otu_clean <- apply(mh_TE_otu.df, 2, function(x) as.numeric(sprintf("%.3f", x)))

row.names(mh_TE_otu_clean) <- row.names(mh_TE_otu.df)
mh_TE_otu <- otu_table(mh_TE_otu_clean, taxa_are_rows = T)


# Make PS object
mh_TE_GSV.ps <- merge_phyloseq(sample_data(mh_TE_map), mh_TE_otu)
mh_TE_GSV.ps


# ## Make object with transformed data using reads #s msweep abundances ####
# # Multiply the counts by the values in the column "themisto_aligned_reads"
# multiplier <- as.numeric(sample_data(mh_TE_GSV.ps)$themisto_aligned_reads)
# 
# # Define a function to multiply counts by the corresponding multiplier
# multiply_by_sample <- function(otu_counts, multiplier) {
#   return(t(t(otu_counts) * multiplier))
# }

### Filter out low abundance GSVs ####
# Filter taxa based on a minimum threshold of 0.01
#mh_TE_GSV_filtered.ps <- prune_taxa(taxa_sums(mh_TE_GSV.ps) > 0.01, mh_TE_GSV.ps)

# another filtering option 

# # Define a function to check if a taxon has a minimum count greater than 0.01 in any sample
# filter_by_sample <- function(otu_table, threshold) {
#   # Returns TRUE for taxa that have counts > threshold in at least one sample
#   rowSums(otu_table > threshold) > 0
# }
# 
# # Apply the function and prune the taxa
# mh_TE_GSV_filtered.ps <- prune_taxa(
#   filter_by_sample(otu_table(mh_TE_GSV.ps), 0.01),
#   mh_TE_GSV.ps
# )
# 
# 
# 
# ### Apply the multiplier ####
# 
# # Get the OTU table as a matrix
# otu_counts <- as(otu_table(mh_TE_GSV_filtered.ps), "matrix")
# 
# # Apply the custom function to multiply each column by the corresponding multiplier
# otu_counts_transformed <- multiply_by_sample(otu_counts, multiplier)
# 
# # Remove scientific notation by formatting the matrix with sprintf to control the number of decimal places
# otu_counts_transformed <- apply(otu_counts_transformed, 2, function(x) as.numeric(sprintf("%.0f", x)))
# 
# # Reassign correct GSV ID's
# rownames(otu_counts_transformed) <- rownames(otu_counts)
# 
# # Create a new OTU table with the transformed counts
# otu_table_transformed <- otu_table(otu_counts_transformed, taxa_are_rows = TRUE)
# 
# # Replace the OTU table in the original phyloseq object
# mh_TE_GSV_transformed.ps <- merge_phyloseq(sample_data(mh_TE_map), otu_table_transformed)
# 
# 
# 
# 
# # Extract the OTU table
# otu_table_transformed <- otu_table(mh_TE_GSV_transformed.ps)
# 
# # Ensure that the OTU table is in matrix format and taxa are rows
# otu_table_matrix <- as(otu_table_transformed, "matrix")
# 
# # Calculate the number of taxa with counts greater than 0.0 in each sample
# taxa_per_sample <- apply(otu_table_matrix, 2, function(x) sum(x > 0))
# 
# # Convert the result to a data frame for easy viewing
# taxa_per_sample_df <- data.frame(Sample = names(taxa_per_sample), TaxaCount = taxa_per_sample)
# 
# # Print the data frame
# print(taxa_per_sample_df)

#### Subset samples #####
# Subset samples where the "Concentration" column has the value "R"
mh_TE_GSV_R_samples <- subset_samples(mh_TE_GSV.ps, Concentration == "R")

# Prune taxa to remove those with zero counts across all samples
mh_TE_GSV_transformed.ps <- prune_taxa(taxa_sums(mh_TE_GSV_R_samples) > 0, mh_TE_GSV_R_samples)

mh_TE_map <- as(sample_data(mh_TE_GSV_transformed.ps), "data.frame")

