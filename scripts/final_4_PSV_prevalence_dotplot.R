library(phyloseq)
library(ggplot2)
library(stringr)
# Assume mh_TE_GSV_transformed.ps is your phyloseq object
# Extract the OTU table from the phyloseq object
mh_TE_GSV_transformed_RA.ps <- transform_sample_counts(mh_TE_GSV_transformed.ps, function(x) {x/sum(x)}*100)

otu_table <- otu_table(mh_TE_GSV_transformed_RA.ps)

# Convert OTU table to a data frame
otu_df <- as.data.frame(as.matrix(otu_table))
otu_df$GSV <- rownames(otu_df)

# Convert to long format
otu_long <- melt(otu_df, id.vars = "GSV", variable.name = "Sample", value.name = "RelativeAbundance")

# Save to CSV file
#write.csv(otu_long, "GSV_relative_abundance_long_format.csv", row.names = FALSE)
# Round RelativeAbundance to two significant figures

otu_long <- otu_long[otu_long$RelativeAbundance > 0, ]

otu_long$RelativeAbundance <- round(otu_long$RelativeAbundance, 3)


# Plot the scatter plot
ggplot(otu_long, aes(x = GSV, y = RelativeAbundance, color = GSV)) +
  geom_point(alpha = 0.6, size = 3) +
  labs(title = "Relative Abundance of GSVs Across Samples", x = "GSV", y = "Relative Abundance", color = "GSV") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0,  size = 10))
  #scale_y_log10()  # Optional: Log scale for better visualization if data has a wide range


otu_long$GSV = factor(otu_long$GSV, levels = c("1","2","5","8","3","7"))

# Define the colors and assign them to the factor levels
tableau8_colors <- c(
  "3" = "#1F77B4", 
  "2" = "#AEC7E8", 
  "1" = "#FF7F0E", 
  #"4" = "#F7B6D2", 
  "5" = "#DBDB8D", 
  "7" = "#98DF8A", 
  "8" = "#F7B6D2"
)
# --- Order GSVs numerically (works whether GSV is "1","2",... or "GSV_1","GSV 2", etc.)
levs <- otu_long %>%
  mutate(GSV_num = as.numeric(str_extract(as.character(GSV), "\\d+"))) %>%
  arrange(GSV_num) %>%
  distinct(GSV) %>%
  pull(GSV)

otu_long <- otu_long %>%
  mutate(GSV = factor(GSV, levels = levs))

# --- Plot
GSV_rel_abund.plot <- ggplot(otu_long, aes(x = GSV, y = RelativeAbundance, color = GSV)) +
  geom_jitter(alpha = 0.8, size = 3, width = 0.1, height = 0) +
  labs( x = "GSV", y = "Relative Abundance", color = "GSV") +
  scale_color_manual(values = tableau8_colors) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, size = 10),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8), # box around plot
    panel.grid.major.x = element_blank(),   # remove vertical major grid lines
    panel.grid.minor.x = element_blank(),   # remove vertical minor grid lines
    panel.grid.minor = element_blank()
  )

# Specify the output file and resolution
png("../Writing/June_2025_draft/figures/Fig_4_GSV_Rel_abund_across_samples.png", width = 2000, height = 1500, res = 300)
GSV_rel_abund.plot
# Turn off the device to save the file
dev.off()


svg("../Writing/June_2025_draft/figures/Fig_4_GSV_Rel_abund_across_samples.svg", width = 16, height = 10)  # open SVG device
print(GSV_rel_abund.plot)                             # draw the plot
dev.off() 

# ### Plot with prevalence line ####
# 
# # Calculate the frequency of each GSV
# GSV_counts <- aggregate(RelativeAbundance ~ GSV, data = otu_long, FUN = function(x) sum(x > 0))
# GSV_counts$RelativeAbundance <- paste0("n=", GSV_counts$RelativeAbundance)
# 
# # Plot the scatter plot with jitter and text annotations
# ggplot(otu_long, aes(x = reorder(GSV, -RelativeAbundance), y = RelativeAbundance, color = GSV)) +
#   geom_jitter(alpha = 0.6, size = 3, width = 0.2, height = 0) +
#   geom_text(data = GSV_counts, aes(x = GSV, y = 0.7, label = RelativeAbundance), color = "black", vjust = 1.5) +
#   labs(title = "Relative Abundance of GSVs Across Samples", x = "GSV", y = "Relative Abundance", color = "GSV") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 0, size = 10))
#   #scale_y_log10()  # Optional: Log scale for better visualization if data has a wide range
# 
# 
# 
# 
# # Calculate the frequency of each GSV
# GSV_counts <- aggregate(RelativeAbundance ~ GSV, data = otu_long, FUN = function(x) sum(x > 0))
# 
# # Create a secondary y-axis scale
# sec_y_max <- max(GSV_counts$RelativeAbundance)
# primary_y_max <- max(otu_long$RelativeAbundance)
# scale_factor <- primary_y_max / sec_y_max
# 
# ggplot(otu_long, aes(x = reorder(GSV, -RelativeAbundance), y = RelativeAbundance, color = GSV)) +
#   geom_jitter(alpha = 0.6, size = 3, width = 0.2, height = 0) +
#   geom_line(data = GSV_counts, aes(x = GSV, y = RelativeAbundance * scale_factor, group = 1), color = "black") +
#   geom_text(data = GSV_counts, aes(x = GSV, y = RelativeAbundance * scale_factor, label = paste0("n=", RelativeAbundance)), color = "black", vjust = -0.6, hjust = -0.03) +
#   labs(title = "Relative Abundance of GSVs Across Samples", x = "GSV", y = "Relative Abundance", color = "GSV") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 0, size = 10)) +
#   scale_y_continuous(sec.axis = sec_axis(~ . / scale_factor, name = "GSV prevalence"))
# 
# #### CHATGPT option ####
# library(ggplot2)
# 
# # Define Tableau 20 colors
# tableau20_colors <- c(
#   "#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C", "#98DF8A", 
#   "#D62728", "#FF9896", "#9467BD", "#C5B0D5", "#8C564B", "#C49C94", 
#   "#E377C2", "#F7B6D2", "#7F7F7F", "#C7C7C7", "#BCBD22", "#DBDB8D", 
#   "#17BECF", "#9EDAE5"
# )
# 
# tableau9_colors <- c(
#   "#1F77B4", "#AEC7E8", "#FF7F0E","#F7B6D2","#DBDB8D","#9EDAE5","#98DF8A","#C7C7C7","#C5B0D5")
# 
# as.factor(otu_long$GSV)
# 
# otu_long$GSV = factor(otu_long$GSV, levels = c("sp3","sp8","sp4","sp1","sp2","sp7","sp6","sp5"))
# 
# 
# ggplot(otu_long, aes(x = reorder(GSV, -RelativeAbundance), y = RelativeAbundance, color = GSV)) +
#   geom_jitter(alpha = 0.6, size = 3, width = 0.2, height = 0) +
#   geom_line(data = GSV_counts, aes(x = GSV, y = RelativeAbundance * scale_factor, group = 1), color = "black") +
#   geom_text(data = GSV_counts, aes(x = GSV, y = RelativeAbundance * scale_factor, label = paste0("n=", RelativeAbundance)), color = "black", vjust = -0.6, hjust = -0.03) +
#   labs( x = "GSV", y = "Relative Abundance", color = "GSV") +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 0, size = 10),
#     panel.grid = element_blank(),
#     panel.border = element_rect(color = "black", fill = NA, size = 1)
#   ) +
#   scale_color_manual(values = tableau9_colors) +
#   scale_y_continuous(sec.axis = sec_axis(~ . / scale_factor, name = "GSV prevalence"))
# 
# 
# ggplot(otu_long, aes(x = reorder(GSV, -RelativeAbundance), y = RelativeAbundance, color = GSV)) +
#   geom_jitter(alpha = 0.6, size = 3, width = 0.2, height = 0) +
#   geom_line(data = GSV_counts, aes(x = GSV, y = RelativeAbundance * scale_factor, group = 1), color = "black") +
#   geom_text(data = GSV_counts, aes(x = GSV, y = RelativeAbundance * scale_factor, label = paste0("n=", RelativeAbundance)), color = "black", vjust = -0.65, hjust = -0.03) +
#   labs(x = "GSV", y = "Relative Abundance", color = "GSV") +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 0, size = 10),
#     panel.grid.major.y = element_line(color = "lightgray", size = 0.5, linetype = "solid"),
#     panel.grid.minor = element_blank(),
#     panel.grid.major.x = element_blank(),
#     panel.border = element_rect(color = "black", fill = NA, size = 1)
#   ) +
#   scale_color_manual(values = tableau9_colors) +
#   scale_y_continuous(sec.axis = sec_axis(~ . / scale_factor, name = "GSV prevalence"))
