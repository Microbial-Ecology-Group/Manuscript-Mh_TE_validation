# Load necessary libraries
library(ape)
library(ggtree)
library(dplyr)

# Read the tree file
mh_tree <- ggtree::read.tree("tree_results/Final_tree_min0.dnd")

# Read the metadata file
metadata <- read.table("tree_results/metadata_genbank_MH_contig_20240528.txt", sep = "\t", header = TRUE)

# Fill missing infraspecific_name with empty strings
metadata$infraspecific_name[is.na(metadata$infraspecific_name)] <- ""

# Ensure the metadata has a column named 'tip_label' that matches the tip labels in the tree
# Assuming the column in metadata to match with the tree tip labels is named 'id'
# Adjust the column name 'id' to the correct column if necessary
metadata <- metadata %>%
  rename(tip_label = Name)

# Create a data frame with tree tip labels
tip_labels_df <- data.frame(tip_label = mh_tree$tip.label)


# Merge the tree tip labels with the metadata
#merged_data <- left_join(tip_labels_df, metadata, by = c("tip_label","Name"))

# Merge the tree tip labels with the metadata by 'tip_label' and 'Name'
merged_data <- full_join(as.data.frame(tip_labels_df), metadata, by = c( "tip_label"))

merged_data$infraspecific_name[is.na(merged_data$infraspecific_name)] <- ""

# Ensure the merged_data is in the correct order for ggtree
merged_data <- merged_data[match(mh_tree$tip.label, merged_data$tip_label), ]

# Add infraspecific_name to mh_tree
mh_tree$tip.label <- merged_data$infraspecific_name

# Plot the tree with ggtree and label tips using 'infraspecific_name', adjusting label size
p <- ggtree(mh_tree, layout = "fan") +
  geom_tiplab(aes(label = label), size = 3, open.angle = 15, hjust = 1, yscale = 20) # Adjust the size value as needed


# Print the plot
print(p)
