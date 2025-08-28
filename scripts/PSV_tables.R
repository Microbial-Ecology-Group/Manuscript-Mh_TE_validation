# Extract metadata from the phyloseq object
metadata <- as.data.frame(sample_data(mh_TE_PSV.ps))

# Create a combined 2x2x2 table
combined_table <- table(metadata$AnyCulPos, metadata$AnyPCRPos, metadata$TEPos)

# Convert the 3-dimensional table to a 2-dimensional data frame for better visualization
combined_df <- as.data.frame(combined_table)
colnames(combined_df) <- c("AnyCulPos", "AnyPCRPos", "TEPos", "Count")

# Print the combined table
print("Combined 2x2x2 table:")
print(combined_df)
