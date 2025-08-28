# Load required libraries
library(reshape2)
library(dplyr)
library(tidyr)


# Read ANI data
MhANI <- read.table("tree_results/mh_2435_fastani_output.txt")
MhANI.mat <- acast(MhANI, V1~V2, value.var = "V3")


# Change to hclust format
MhANI.hclust <- hclust(dist(MhANI.mat), method = "complete")

# Initialize an empty list to store the cluster columns
cluster_list <- list()

# Get the row names from the original data to ensure consistent ordering
row_names <- rownames(MhANI.mat)

# Loop over k values from 2 to 30
for (k in 2:30) {
  # Perform hierarchical clustering and cut the tree at k clusters
  my_gene_col <- cutree(tree = MhANI.hclust, k = k)
  
  # Convert the cluster assignments to a dataframe with the correct row names
  my_gene_df <- data.frame(cluster = as.factor(my_gene_col), row.names = names(my_gene_col))
  
  # Ensure the dataframe has the correct order of row names
  my_gene_df <- my_gene_df[row_names, , drop = FALSE]
  
  # Add the cluster assignments to the list
  cluster_list[[paste0("k_", k)]] <- my_gene_df$cluster
}

# Combine all cluster assignments into a single dataframe
final_df <- as.data.frame(cluster_list)

# Set the row names to ensure they match the original data
rownames(final_df) <- row_names

# Inspect the final dataframe
print(final_df)

# Write the dataframe to a TSV file
## Remember that this order of genomes might be different than their order in the themisto database.
# N.B. If you don't change the order of this spreadsheet to match the database all your results will be wrong!
write.table(final_df, file = "final_2025_ANI_groupings_2-30.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)


## Calculate within group ANI ####
# Extract the k=9 cluster assignments
# Extract the k=8 cluster assignments as a named vector
clusters_k8 <- final_df$k_8
names(clusters_k8) <- rownames(final_df)  # assign genome names

# Initialize a vector to store average within-cluster ANI
avg_within_ANI <- numeric(length(unique(clusters_k8)))

# Loop over each cluster
for (clust in unique(clusters_k8)) {
  # Get genomes in this cluster
  genomes_in_cluster <- names(clusters_k8[clusters_k8 == clust])
  
  # Subset the ANI matrix for these genomes
  sub_ANI <- MhANI.mat[genomes_in_cluster, genomes_in_cluster]
  
  # Get the upper triangle values (excluding diagonal)
  upper_vals <- sub_ANI[upper.tri(sub_ANI)]
  
  # Compute the mean
  avg_within_ANI[as.character(clust)] <- mean(upper_vals)
}

# View results
avg_within_ANI

### Alternative way to calculate within GSV ANI averages ####
# Make sure cluster assignments are named by genome
clusters_k8 <- final_df$k_8
names(clusters_k8) <- rownames(final_df)

# Convert the ANI matrix to long format
ani_long <- as.data.frame(as.table(as.matrix(MhANI.mat))) %>%
  rename(Genome1 = Var1, Genome2 = Var2, ANI = Freq)

# Add cluster assignment for each genome
ani_long <- ani_long %>%
  mutate(Cluster1 = clusters_k8[Genome1],
         Cluster2 = clusters_k8[Genome2]) %>%
  # Keep only pairs within the same cluster, excluding self-comparisons
  filter(Cluster1 == Cluster2 & Genome1 != Genome2)

# Compute average within-cluster ANI
avg_within_ANI <- ani_long %>%
  group_by(Cluster1) %>%
  summarise(avg_ANI = mean(ANI), .groups = "drop")

# View results
avg_within_ANI


## Calculate ANI distances between representative genomes ####
# List of genomes of interest
genomes_of_interest <- c(
  "/scratch/user/enriquedoster/2025_Mh_validation/20250110_Mh_Genbank_minus_outliers/GCA_007965105.1_ASM796510v1_genomic.fna",
  "/scratch/user/enriquedoster/2025_Mh_validation/20250110_Mh_Genbank_minus_outliers/GCA_039617675.1_ASM3961767v1_genomic.fna",
  "/scratch/user/enriquedoster/2025_Mh_validation/20250110_Mh_Genbank_minus_outliers/GCA_033335735.1_ASM3333573v1_genomic.fna",
  "/scratch/user/enriquedoster/2025_Mh_validation/20250110_Mh_Genbank_minus_outliers/GCA_030370375.1_ASM3037037v1_genomic.fna",
  "/scratch/user/enriquedoster/2025_Mh_validation/20250110_Mh_Genbank_minus_outliers/GCA_023516315.1_ASM2351631v1_genomic.fna",
  "/scratch/user/enriquedoster/2025_Mh_validation/20250110_Mh_Genbank_minus_outliers/GCA_900638155.1_56527_B01_genomic.fna",
  "/scratch/user/enriquedoster/2025_Mh_validation/20250110_Mh_Genbank_minus_outliers/GCA_046115575.1_ASM4611557v1_genomic.fna",
  "/scratch/user/enriquedoster/2025_Mh_validation/20250110_Mh_Genbank_minus_outliers/GCA_038145235.1_ASM3814523v1_genomic.fna"
)

# Subset the ANI matrix for these genomes (rows AND columns)
sub_ANI <- MhANI.mat[genomes_of_interest, genomes_of_interest]

# Convert to a tidy dataframe of pairwise comparisons
pairwise_ANI <- as.data.frame(as.table(sub_ANI))

# Rename columns for clarity
colnames(pairwise_ANI) <- c("Genome1", "Genome2", "ANI")

# Remove self-comparisons if you want only true pairwise comparisons
pairwise_ANI <- pairwise_ANI[pairwise_ANI$Genome1 != pairwise_ANI$Genome2, ]

# Optional: keep only one direction (Genome1 < Genome2) to avoid duplicates
pairwise_ANI <- pairwise_ANI[!duplicated(t(apply(pairwise_ANI[,1:2], 1, sort))), ]

# Inspect results
pairwise_ANI

#write.table(pairwise_ANI, file = "final_2025_ANI_by_representative.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

