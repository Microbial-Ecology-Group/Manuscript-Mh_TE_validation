# Relative abundance plots for GSV data
# The only thing you need as input for this script is a phyloseq object with the GSV data (mhpp_GSV.ps)

# For the part of the plot showing the percent of counts classified as GSVs, you'll need a column in your
# metadata with this percentage already calculated and listed as a proportion (percent_mh_of_raw).

# Load library
library(ggsci) # if not using this package, make sure to switch out the palatte in the figures
library(metagMisc)
library(dplyr)
library(ggdendro)
library(cowplot) # Used for making the grid plots


#####  CSS TRANSFORM    ####
###
# mh_TE_GSV_transformed.ps <- prune_taxa(taxa_sums(mh_TE_GSV_transformed.ps) > 0, mh_TE_GSV_transformed.ps)
# any(taxa_sums(mh_TE_GSV_transformed.ps)==0) # QUADRUPLE CHECKING - nope good.
# 
# mh_TE_GSV_transformed.ps.df <- as(sample_data(mh_TE_GSV_transformed.ps), "data.frame")

mh_TE_GSV_transformed.ps.df   <- as(sample_data(mh_TE_GSV_transformed.ps), "data.frame")
###
#####  RELATIVE ABUNDANCE   #####
###

# Transform CSS counts to rel abundance
rel_abund <- transform_sample_counts(mh_TE_GSV_transformed.ps, function(x) {x/sum(x)}*100)

# melt data
rel_abund_melt <- psmelt(rel_abund)


###### Re-label low abundance taxa into single category #####
# rel_abund_melt_agg <- rel_abund_melt %>%
#   group_by(Sample, GSV) %>%  # Group by both Sample and GSV
#   summarise(Abundance = sum(Abundance)) %>%  # Sum the Abundance for each group
#   ungroup()

# By pass labeling
rel_abund_melt_agg <- rel_abund_melt 



##
###### Clustering at ASV level ####
##

# Cluster using the function hclust() and default settings
ps_mhpp_GSV.dist <- vegdist(t(otu_table(mh_TE_GSV_transformed.ps)), method = "bray")
ps_mhpp_GSV.hclust <- hclust(ps_mhpp_GSV.dist)
plot(ps_mhpp_GSV.hclust) # example plot

# Extract data as dendrogram
ps_mhpp_GSV.dendro <- as.dendrogram(ps_mhpp_GSV.hclust)
ps_mhpp_GSV.dendro.data <- dendro_data(ps_mhpp_GSV.dendro, type = "rectangle")
ps_mhpp_GSV.dendro.data #  this object contains $segments and $labels

# Sample names in order based on clustering
sample_names_ps_mhpp_GSV <- ps_mhpp_GSV.dendro.data$labels$label


# Add metadata #### 
# make sure that the column (or change name) in your phyloseq object matches with this column: ps_mhpp_GSV.dendro.data$labels
ps_mhpp_GSV.dendro.data$labels <- left_join(ps_mhpp_GSV.dendro.data$labels, sample_data(mh_TE_GSV_transformed.ps), by = c("label" = "sample_id"))

# Setup the data, so that the layout is inverted (this is more 
# "clear" than simply using coord_flip())
segment_data_ps_mhpp_GSV <- with(
  segment(ps_mhpp_GSV.dendro.data),
  data.frame(x = y, y = x, xend = yend, yend = xend))

# Use the dendrogram label data to position the gene labels
gene_pos_table_ps_mhpp_GSV <- with(
  ps_mhpp_GSV.dendro.data$labels, 
  data.frame(y_center = x, gene = as.character(label), x = y , Level = as.character(Mh_level), 
             Concentration = as.character(Concentration), percent_mh_of_raw = percent_mh_of_raw,  
             Nasal.Mh.Culture = as.character(Nasal.Mh.Culture),NS.Mh.PCR = as.character(NS.Mh.PCR),
             BRD = as.character(BRD), DGNP.Mh.PCR = as.character(DGNP.Mh.PCR),
             height = 1))


# Table to position the samples
sample_pos_table_ps_mhpp_GSV <- data.frame(sample = sample_names_ps_mhpp_GSV) %>%
  dplyr::mutate(x_center = (1:dplyr::n()),  width = 1)

##
######## Relative abundance bar plot #########
##

# Use class melted data and add gene locations
joined_ra_GSV_ps_GSV_melt <- rel_abund_melt_agg %>%
  left_join(gene_pos_table_ps_mhpp_GSV, by = c("Sample" = "gene")) %>%
  left_join(sample_pos_table_ps_mhpp_GSV, by = c("Sample" = "sample")) 

# Calculate the mean relative abundance of each class taxa, sort by most abundant to least
factor_by_abund <- rel_abund_melt %>%
  dplyr::group_by(OTU) %>%
  dplyr::summarize(median_GSV = median(Abundance)) %>%
  arrange(-median_GSV)

# Use the sorted class levels to clean up the order of taxa in the relative abundance plots
joined_ra_GSV_ps_GSV_melt$GSV <- factor(joined_ra_GSV_ps_GSV_melt$OTU, levels = as.character(factor_by_abund$OTU))

# Limits for the vertical axes
gene_axis_limits_ps_class <- with(
  gene_pos_table_ps_mhpp_GSV, 
  c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))) +   0.1 * c(-1, 1) # extra spacing: 0.1


# Useful 20 colors + black and white (https://sashamaps.net/docs/resources/20-colors/)
#'#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#a9a9a9', '#ffffff', '#000000'

# Option to change name
#joined_ra_GSV_ps_GSV_melt <- joined_ra_GSV_ps_GSV_melt %>%
#  mutate(GSV = str_replace_all(GSV, "Mannheimia_haemolytica_strain_", "Mh "))

factor_by_abund <- joined_ra_GSV_ps_GSV_melt %>%
  dplyr::group_by(GSV) %>%
  dplyr::summarize(median_GSV = median(Abundance)) %>%
  arrange(-median_GSV)

# Use the sorted class levels to clean up the order of taxa in the relative abundance plots
joined_ra_GSV_ps_GSV_melt$GSV <- factor(joined_ra_GSV_ps_GSV_melt$GSV, levels = as.character(factor_by_abund$GSV))

# Change order of Mh levels in the "Level" column
gene_pos_table_ps_mhpp_GSV$Level = factor(gene_pos_table_ps_mhpp_GSV$Level, levels = c("zero","low","medium","high","highest"))

# Dendrogram plot ####
plt_dendr_GSV <- ggplot(segment_data_ps_mhpp_GSV) + 
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend), size=.75, lineend="round", linejoin="round") +
  geom_point(data=gene_pos_table_ps_mhpp_GSV, aes(x, y_center, colour=Level , shape=Level ),
             size=8, stroke=0.75, position=position_nudge(x=-0.03)) +
  scale_y_continuous(breaks=gene_pos_table_ps_mhpp_GSV$y_center, 
                     labels=gene_pos_table_ps_mhpp_GSV$gene, 
                     limits=gene_axis_limits_ps_class, 
                     expand=c(0, 0)) + 
  labs(x="Ward's Distance", y="", title="") +
  scale_x_reverse() + 
  #scale_colour_manual(name="BRD diagnosis", values = c("dark gray", "red","#6497B1", "#005B96", "#03396C", "#011F4B")) +  # Set legend title and colors
  scale_colour_manual(name="16S Mh level", values = c("yellow","light gray", "#6497B1", "#005B96", "#03396C", "#011F4B")) +  # Set legend title and colors
  scale_shape_manual(values=c(15,15,15,15,15,15), guide="none") +  
  theme_bw() + 
  theme(legend.position="right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16),
        panel.border=element_blank(),
        plot.title=element_text(size=30),
        panel.grid=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.x=element_line(linewidth=0.75),
        axis.title.x=element_text(size=14),
        axis.text.x=element_text(size=12, colour="black"))
plt_dendr_GSV

#
## Relative abundance plot ####
#
plt_rel_GSV <- ggplot(joined_ra_GSV_ps_GSV_melt, 
                      aes(x = x_center, y = Abundance, fill = GSV, 
                          height = height, width = width)) + 
  coord_flip() +
  geom_bar(stat = "identity", colour = "black") +
  #scale_fill_manual(values = col_vector) + #use this if not using color palette from ggthemes of ggsci
  scale_fill_tableau(palette = "Tableau 20")+
  scale_x_continuous(breaks = sample_pos_table_ps_mhpp_GSV$x_center, 
                     labels = sample_pos_table_ps_mhpp_GSV$sample, 
                     expand = c(0, 0)) + 
  # For the y axis, alternatively set the labels as: gene_position_table$gene
  #scale_y_continuous(breaks = gene_pos_table[, "y_center"], labels = rep("", nrow(gene_pos_table)),limits = gene_axis_limits, expand = c(0, 0)) + 
  labs(x = "", y = "GSV Relative abundance") +
  theme_bw() +
  theme(legend.position = "right",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.75),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(linewidth = 0.75, lineend = "square", colour = "black"),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12, colour = "black"),
        plot.margin = unit(c(1, 0.2, 0.2, -0.7), "cm"), 
        panel.grid.minor = element_blank())
plt_rel_GSV


## Make plot showing relative abundance of GSV hits to nonhost counts ####

# Add a column for the remaining percentage to make up 100%
mh_TE_GSV_transformed.ps.df$remainder_percentage <- 100 - mh_TE_GSV_transformed.ps.df$percent_mh_of_raw

# Reshape data to long format for ggplot2
long_data <- tidyr::pivot_longer(mh_TE_GSV_transformed.ps.df, cols = c("percent_mh_of_raw", "remainder_percentage"), 
                                 names_to = "category", values_to = "value")

long_data$category <- factor(long_data$category, levels = c("remainder_percentage","percent_mh_of_raw"))

long_data$sample_id <- factor(long_data$sample_id, levels = as.factor(gene_pos_table_ps_mhpp_GSV$gene))

plt_small_bar <- ggplot(long_data, aes(x = sample_id, y = value, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = ifelse(category == "percent_mh_of_raw", sprintf("%.1f%%", value), "")), 
            position = position_stack(vjust = 0.5), hjust = -0.2,
            size = 3,  # Adjust this value for size
            fontface = "bold") +
  coord_flip() +
  scale_y_continuous(
    breaks = c(0, 100),
    labels = scales::percent_format(scale = 1)(c(0, 100))
  ) +
  scale_fill_manual(values = c("percent_mh_of_raw" = "darkred", "remainder_percentage" = "darkgrey")) +
  labs(y = "Mh %", fill = "") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"))

# Display the plot
plt_small_bar
#plots_dendro_rel_GSV <- plot_grid(plt_dendr_GSV, plt_rel_GSV, align = 'h', rel_widths = c(0.5, 1.5))



# Make the combined plot #####

# Extract the legends from each plot
legend_dendr_GSV <- get_legend(plt_dendr_GSV)
legend_rel_GSV <- get_legend(plt_rel_GSV)

# Remove the legends from the original plots
plt_dendr_GSV <- plt_dendr_GSV + theme(legend.position = "none")
plt_rel_GSV <- plt_rel_GSV + theme(legend.position = "none")

combined_legend <- ggdraw() +
  draw_plot(legend_dendr_GSV, x = .072, y = 0.55, width = 0.9, height = 0.3) +  # Adjust the x value to shift the top legend to the left
  draw_plot(legend_rel_GSV, x = 0, y = 0, width = 1, height = 0.7)

combined_legend

plots_dendro_rel_GSV <- plot_grid(
  plt_dendr_GSV, 
  plt_rel_GSV, 
  combined_legend,
  align = 'h', 
  ncol = 3,
  rel_widths = c(0.3,1.1, 0.3)
)

plots_dendro_rel_GSV

### Combined plot with mh relative abundance
plots_dendro_rel_GSV_wbar <- plot_grid(
  plt_dendr_GSV, 
  plt_small_bar,
  plt_rel_GSV, 
  combined_legend,
  align = 'h', 
  ncol = 4,
  rel_widths = c(0.15,0.07 , 0.85, 0.1)
)

plots_dendro_rel_GSV_wbar


#ggsave("mh_GSV_dendrogram_rel_abun_by_Sample_type.jpg", plot=last_plot(), device="jpeg", path="figures/", scale = .5, units = "in", dpi = 1600)

svg("figures/mh_GSV_dendrogram_rel_abun_by_BRD.svg", width = 16, height = 10)  # open SVG device
print(plots_dendro_rel_GSV_wbar)                             # draw the plot
dev.off() 


# # This requires a directory called figures to be in your working directory
# # All 3 lines have to be run for the file to be created.
png("figures/mh_GSV_dendrogram_rel_abun_by_BRD.png", width = 1400, height = 800)
plots_dendro_rel_GSV_wbar
dev.off()

