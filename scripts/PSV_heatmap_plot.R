# Complex heatmap

library(ComplexHeatmap)
library(circlize)
library(metagMisc)

# Extract the OTU table as a matrix, and ensure it's in the right orientation
otu_matrix <- as(otu_table(mh_TE_PSV_transformed.ps), "matrix")


# If taxa are in rows, you can transpose the matrix if needed
if (taxa_are_rows(mh_TE_PSV_transformed.ps)) {
  otu_matrix <- t(otu_matrix)
}

# Log-transform, adding 1 to avoid log(0)
log_otu_matrix <- log(otu_matrix + 1)

# Flatten the matrix to calculate the median
log_otu_vector <- as.vector(log_otu_matrix)

# Calculate the median
median_value <- median(log_otu_vector)

# Define custom color mapping, ensuring 0 maps to white
col_fun <- colorRamp2(c(0, 5, max(log_otu_vector)), c("white", "blue", "red"))

# Create the heatmap with 0s clearly marked as white
Heatmap(log_otu_matrix, name = "LogAbundance", col = col_fun)


# Extract sample data (annotation data)
sample_annotation <- as(sample_data(mh_TE_PSV_transformed.ps), "data.frame")

# Ensure that the order of rows in the OTU matrix matches the annotation data
log_otu_matrix <- log(otu_matrix + 1)
log_otu_matrix <- log_otu_matrix[order(sample_annotation$Mh_level), ]  # Order rows by Mh_level
sample_annotation <- sample_annotation[order(sample_annotation$Mh_level), ]  # Sort sample annotation to match

# Define the annotation for percent_mh_of_raw as a barplot
row_anno <- rowAnnotation(
  percent_mh_of_raw = anno_barplot(sample_annotation$percent_mh_of_raw)
)

# Define custom color mapping with zero clearly white
col_fun <- colorRamp2(c(0, 10, max(log_otu_matrix)), c("white", "blue", "red"))

# # Create the heatmap with annotations and ordered by Mh_level
# Heatmap(log_otu_matrix, name = "LogAbundance", col = col_fun, 
#         right_annotation = row_anno, 
#         show_row_names = TRUE, cluster_rows = FALSE)
# 

# Ensure the OTU matrix is numeric
log_otu_matrix <- as.matrix(log_otu_matrix)
log_otu_matrix <- apply(log_otu_matrix, 2, as.numeric)

# Ensure sample_annotation has no NA values
sample_annotation <- sample_annotation[complete.cases(sample_annotation$percent_mh_of_raw, sample_annotation$Mh_level), ]

# Ensure row names of log_otu_matrix match those of sample_annotation
rownames(log_otu_matrix) <- rownames(sample_annotation)


# Define the Mh_level colors
mh_level_colors <- c("zero" = "yellow", 
                     "low" = "light gray", 
                     "medium" = "#6497B1", 
                     "high" = "#005B96", 
                     "highest" = "#03396C")

# Define the BRD colors
brd_colors <- c("Y" = "red", "N" = "grey")

# Create the left annotation for percent_mh_of_raw
left_anno <- rowAnnotation(
  percent_mh_of_raw = anno_barplot(sample_annotation$percent_mh_of_raw, 
                                   annotation_name_side = "top")
)

# Create the right annotation with correctly rotated labels
right_anno <- rowAnnotation(
  Mh_level = sample_annotation$Mh_level,
  BRD = sample_annotation$BRD,
  col = list(Mh_level = mh_level_colors, BRD = brd_colors),
  annotation_name_gp = gpar(fontsize = 14),  # Adjust font size for annotation labels
  annotation_name_rot = 90  # Rotate the annotation names by 90 degrees
)


ht <- Heatmap(log_otu_matrix, 
              name = "LogAbundance", 
              col = col_fun, 
              left_annotation = left_anno,  # Left annotation
              right_annotation = right_anno,  # Right annotations
              show_row_names = FALSE,  # Remove row names
              column_names_side = "top",  # Column labels on top
              column_names_rot = 45,  # Rotate column labels to 45 degrees
              column_names_gp = gpar(fontsize = 14),  # Font size for column names
              cluster_rows = FALSE,  # Maintain row order
              rect_gp = gpar(col = "black", lwd = 0.5, alpha = .9)  # Add borders (lines) between cells
)
# # Custom function to draw only vertical lines in the heatmap
# draw_vertical_borders <- function(j, i, x, y, width, height, fill) {
#   # Draw only the vertical lines (on the left side of each cell)
#   grid.segments(
#     x0 = x - width / 2, x1 = x - width / 2, 
#     y0 = y - height / 2, y1 = y + height / 2, 
#     gp = gpar(col = "black", lwd = 1)
#   )
# }
# 
# # Heatmap with custom border drawing function
# ht <- Heatmap(log_otu_matrix, 
#               name = "LogAbundance", 
#               col = col_fun, 
#               left_annotation = left_anno,  # Left annotation
#               right_annotation = right_anno,  # Right annotations
#               show_row_names = FALSE,  # Remove row names
#               column_names_side = "top",  # Column labels on top
#               column_names_rot = 45,  # Rotate column labels to 45 degrees
#               column_names_gp = gpar(fontsize = 14),  # Font size for column names
#               cluster_rows = FALSE,  # Maintain row order
#               rect_gp = gpar(col = NA),  # Remove default borders
#               cell_fun = draw_vertical_borders  # Add custom vertical borders
# )
draw(ht, padding = unit(c(2, 14, 2, 2), "mm"))  # Increase the left margin (last value)

# Specify the output file and resolution
png("../Writing/Figures/heatmap_PSV_abund.png", width = 2000, height = 1500, res = 300)

draw(ht, padding = unit(c(2, 14, 2, 2), "mm"))  # Increase the left margin (last value)

# Turn off the device to save the file
dev.off()
