# PSV 
library(cowplot)

# Calculate alpha diversity and update metadata file in phyloseq object with new results
mh_TE_PSV.div <- estimate_richness(mh_TE_PSV_transformed.ps, measures = c("Observed","Shannon","Simpson","InvSimpson"))

# Making sure row names match with phyloseq objects
rownames(mh_TE_PSV.div) <- sample_names(mh_TE_PSV_transformed.ps)

# Extract the existing sample data from the phyloseq object
existing_sample_data <- sample_data(mh_TE_PSV_transformed.ps)

# Merge alpha diversity data with the existing sample data
new_sample_data <- merge(existing_sample_data, mh_TE_PSV.div, by = "row.names", all = TRUE)
row.names(new_sample_data) <- new_sample_data$Row.names
new_sample_data$Row.names <- NULL

# Update the sample data in the phyloseq object
sample_data(mh_TE_PSV_transformed.ps) <- sample_data(new_sample_data)

sample_data(mh_TE_PSV_transformed.ps)$Mh_level <- factor(sample_data(mh_TE_PSV_transformed.ps)$Mh_level , levels = c("zero","low","medium","high","highest"))


# Plot richness and diversity
# Richness by sample type
plot1 <- ggplot(sample_data(mh_TE_PSV_transformed.ps), mapping=aes(x= Mh_level , y = Observed, fill = BRD , color = BRD )) +
  theme_bw() + 
  labs(y = "Observed PSVs") +
  geom_boxplot(alpha=0.35, position=position_dodge(width=0.8)) +
  geom_jitter(position=position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8)) +
  geom_smooth() +
  scale_fill_manual(values = c("darkred", "darkcyan")) +
  scale_color_manual(values = c("darkred", "darkcyan")) +
  theme(legend.position = "right",
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        legend.title = element_text("Sample Type", size = 10),
        strip.background = element_rect(fill= "black", size = 1.0),
        strip.text = element_text(size =24, colour = "white"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.7),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 18),
        panel.border = element_rect(colour = "black", linewidth = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

plot1
#ggsave("figures/mh_PSV_richness_barplot.jpg", plot=last_plot(), device="jpeg", path="figures/", scale = 1, units = "in", dpi = 600)


# PW wilcoxon tests of richness by type
mhpp_pw_obs_Mh_level <- pairwise.wilcox.test(sample_data(mh_TE_PSV_transformed.ps)$Observed, sample_data(mh_TE_PSV_transformed.ps)$Mh_level, p.adjust.method = "BH")
mhpp_pw_obs_Mh_level$p.value # p-values for Observed PSVs

mhpp_pw_obs_BRD <- pairwise.wilcox.test(sample_data(mh_TE_PSV_transformed.ps)$Observed, sample_data(mh_TE_PSV_transformed.ps)$BRD, p.adjust.method = "BH")
mhpp_pw_obs_BRD$p.value # p-values for Shannons



## Shannon's ####

# Plot richness and diversity
# Richness by sample type
plot2 <-ggplot(sample_data(mh_TE_PSV_transformed.ps), mapping=aes(x= Mh_level , y = Shannon, fill = BRD , color = BRD )) +
  theme_bw() + 
  labs(y = "Shannon's index") +
  geom_boxplot(alpha=0.35, position=position_dodge(width=0.8)) +
  geom_jitter(position=position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8)) +
  geom_smooth() +
  scale_fill_manual(values = c("darkred", "darkcyan")) +
  scale_color_manual(values = c("darkred", "darkcyan")) +
  theme(legend.position = "right",
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        legend.title = element_text("", size = 10),
        strip.background = element_rect(fill= "black", size = 1.0),
        strip.text = element_text(size =24, colour = "white"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.7),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 18),
        panel.border = element_rect(colour = "black", linewidth = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
plot2
#ggsave("figures/mh_PSV_shannon_barplot.jpg", plot=last_plot(), device="jpeg", path="figures/", scale = 1, units = "in", dpi = 600)

# PW wilcoxon tests of richness by type
mhpp_pw_shann_Mh_level <- pairwise.wilcox.test(sample_data(mh_TE_PSV_transformed.ps)$Shannon, sample_data(mh_TE_PSV_transformed.ps)$Mh_level, p.adjust.method = "BH")
mhpp_pw_shann_Mh_level$p.value # p-values for Shannons

mhpp_pw_shann_BRD <- pairwise.wilcox.test(sample_data(mh_TE_PSV_transformed.ps)$Shannon, sample_data(mh_TE_PSV_transformed.ps)$BRD, p.adjust.method = "BH")
mhpp_pw_shann_BRD$p.value # p-values for Shannons


# Combine the plots
plot_grid(plot1, plot2, ncol = 2, align = "v")

# Extract metadata
metadata <- data.frame(sample_data(mh_TE_PSV_transformed.ps))

# Reshape the data
long_data <- metadata %>%
  pivot_longer(cols = c(Observed, Shannon),
               names_to = "Variable",
               values_to = "Value")

# Create the plot
# Create the combined plot with vertical facets
ggplot(long_data, aes(x = Mh_level, y = Value, fill = Concentration, color = Concentration)) +
  geom_boxplot(alpha=0.35, position=position_dodge(width=0.8)) +
  geom_jitter(position=position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8)) +
  facet_grid(Variable ~ ., scales = "free", space = "free_x") +
  scale_fill_manual(values = c("darkred", "darkcyan")) +
  scale_color_manual(values = c("darkred", "darkcyan")) +
  theme_bw() +
  labs(y = NULL) +
  theme(legend.position = "right",
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        legend.title = element_text("Sample Type", size = 10),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.7),
        panel.border = element_rect(colour = "black", linewidth = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_rect(fill= "black", size = 1.0),
        strip.text = element_text(size = 14, colour = "white"),
        # Remove x-axis labels and ticks for the top facet
        strip.placement = "outside")
