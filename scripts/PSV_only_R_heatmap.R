# Load necessary libraries
library(phyloseq)
library(ComplexHeatmap)
library(circlize)

sample_data(mh_TE_GSV_transformed.ps)$Mh_level <- factor(sample_data(mh_TE_GSV_transformed.ps)$Mh_level , levels = c("zero","low","medium","high","highest"))

# Extract the OTU table as a matrix, and ensure it's in the right orientation
otu_matrix <- as(otu_table(mh_TE_GSV_transformed.ps), "matrix")

# If taxa are in rows, transpose the matrix if needed
if (taxa_are_rows(mh_TE_GSV_transformed.ps)) {
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

# Extract sample data (annotation data)
sample_annotation <- as(sample_data(mh_TE_GSV_transformed.ps), "data.frame")
mh_TE_map <- as(sample_data(mh_TE_GSV_transformed.ps), "data.frame")

str(sample_annotation)

# Ensure that the order of rows in the OTU matrix matches the annotation data
log_otu_matrix <- log_otu_matrix[order(mh_TE_map$Mh_level), ]  # Order rows by Mh_level
mh_TE_map <- mh_TE_map[order(mh_TE_map$Mh_level), ]  # Sort sample annotation to match

# Ensure mh_TE_map has no NA values in the relevant columns
mh_TE_map <- mh_TE_map[complete.cases(mh_TE_map$percent_mh_of_nonhost, 
                                                      mh_TE_map$Mh_level, 
                                                      mh_TE_map$X16S_qPCR_mean_copy), ]

# Update row names of log_otu_matrix to match those of mh_TE_map
rownames(log_otu_matrix) <- rownames(mh_TE_map)

# Define the Mh_level and BRD colors
mh_level_colors <- c("zero" = "yellow", 
                     "low" = "light gray", 
                     "medium" = "#6497B1", 
                     "high" = "#005B96", 
                     "highest" = "#03396C")
brd_colors <- c("Y" = "red", "N" = "green3")

# Define the annotation for two barplots on the left with rotated labels
# left_anno <- rowAnnotation(
#   `16S qPCR \n copy #` = anno_barplot(mh_TE_map$X16S_qPCR_mean_copy),  # Adjust font size
#   `Mh % of \n raw reads`= anno_barplot(mh_TE_map$percent_mh_of_nonhost),
#   gap = unit(1.5, "mm"), # puts a gap between the bar plots
#   width = unit(3,"cm"), # adjusts the width of each bar plot
#   annotation_name_side = "top", # where the anno label goes
#   annotation_name_gp = gpar(fontsize = 10, just = "right"), # making the anno label fit
#   annotation_name_rot = 90, # angle of the annotation label (0deg in this case)
#   annotation_name_align = TRUE, # this is supposed to align the label to the plot
#   annotation_name_offset = unit(1.5,"mm") # provides separation from the annotation and label 
# )
# Find max values for each bar to reverse the axis
left_anno <- rowAnnotation(
  `16S qPCR \n copy #` = anno_barplot(
    mh_TE_map$X16S_qPCR_mean_copy,
    gp = gpar(fill = "steelblue"),
    bar_width = 0.9,
    flip = TRUE        # <-- this inverts the bar scale
  ),
  `Mh % of \n nonhost reads` = anno_barplot(
    mh_TE_map$percent_mh_of_nonhost,
    gp = gpar(fill = "orange"),
    bar_width = 0.9,
    flip = TRUE        # <-- flip this one too
  ),
  gap = unit(1.5, "mm"),
  width = unit(3,"cm"),
  annotation_name_side = "top",
  annotation_name_gp = gpar(fontsize = 10, just = "right"),
  annotation_name_rot = 90,
  annotation_name_align = TRUE,
  annotation_name_offset = unit(1.5,"mm")
)


# Right annotation with vertically rotated labels for Mh_level and BRD
right_anno <- rowAnnotation(
  `Mh level` = mh_TE_map$Mh_level,
  BRD        = mh_TE_map$BRD,
  col = list(
    `Mh level` = mh_level_colors,  # <- key matches exactly
    BRD        = brd_colors
  ),
  annotation_name_side = "top",
  annotation_name_gp = gpar(fontsize = 10),  # Font size for annotation labels
  annotation_name_rot = c(90, 90)  # Rotate the labels for each annotation vertically
)

# Define color mapping for the heatmap with 0s as white
col_fun <- colorRamp2(c(0, 10, max(log_otu_matrix)), c("white", "blue", "red"))


# Create the heatmap with proper margins
(heat_map_PSV <- Heatmap(log_otu_matrix, name = "Log Abundance", col = col_fun, 
                         left_annotation = left_anno,
                         right_annotation = right_anno,  
                         rect_gp = gpar(col = rgb(0,0,0,alpha=.3), lwd = .5),  # add borders (lines) between cells
                         show_row_names = FALSE,  # remove row names
                         column_names_side = "top",  # column labels on top
                         column_names_rot = 45,  # rotate column labels to 45 degrees
                         column_names_gp = gpar(fontsize = 10),  # font size for column names
                         column_names_centered = TRUE,
                         cluster_rows = FALSE,   # maintain row order
                         show_column_dend = TRUE) #  dendrogram
)


# Specify the output file and resolution
png("../Writing/June_2025_draft/figures/Fig_3_GSV_Heatmap.png", width = 2000, height = 1500, res = 300)
heat_map_PSV
# Turn off the device to save the file
dev.off()


svg("../Writing/June_2025_draft/figures/Fig_3_GSV_Heatmap.svg", width = 16, height = 10)  # open SVG device
print(heat_map_PSV)                             # draw the plot
dev.off() 

