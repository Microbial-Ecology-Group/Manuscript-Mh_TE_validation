library(ggtree)
library(ggplot2)
library(dendextend)
library(ape)
library(tibble)
library(reshape2)
library(tidyverse)

# Read ANI data
MhANI <- read.table("/scratch/user/enriquedoster/2025_Mh_validation/ANI_analysis/mh_2435_fastani_output.txt")
MhANI.mat <- acast(MhANI, V1 ~ V2, value.var = "V3")

# Read metadata (keep only `k_10` and `Assembly.Submitter`)
tree_metadata <- read.table("updated_metadata_mh_genbank.txt", sep = "\t", header = TRUE, fill = NA, stringsAsFactors = FALSE)
tree_metadata <- tree_metadata[, c("Filename", "k_10")] #"Assembly.Submitter"

# Match the 'Filename' column from metadata with column names of the matrix
matching_columns <- colnames(MhANI.mat) %in% tree_metadata$Filename
MhANI.mat_filtered <- MhANI.mat[, matching_columns]

# Extract labels and annotations based on matching 'Filename' from the filtered metadata
col_annotations <- tree_metadata[match(colnames(MhANI.mat_filtered), tree_metadata$Filename), c("k_10", "Filename")]

# Prepare the hierarchical clustering tree (dendrogram) from the filtered matrix
hc <- hclust(dist(MhANI.mat_filtered), method = "complete") # 
tree <- as.phylo(hc)
