# PCV beta diversity


#############################################################################################
####################################   CSS TRANSFORM    ###################################

sample_data(mh_TE_PSV_transformed.ps)$Mh_level <- factor(sample_data(mh_TE_PSV_transformed.ps)$Mh_level , levels = c("zero","low","medium","high","highest"))

mh_TE_PSV_transformed.ps <- prune_taxa(taxa_sums(mh_TE_PSV_transformed.ps) > 0, mh_TE_PSV_transformed.ps)
any(taxa_sums(mh_TE_PSV_transformed.ps)==0) # QUADRUPLE CHECKING - nope good.

mh_TE_PSV_transformed.ps.css <- phyloseq_transform_css(mh_TE_PSV_transformed.ps, log = F)
mh_TE_PSV_transformed.ps.css.df <- as(sample_data(mh_TE_PSV_transformed.ps.css), "data.frame")


###
####
##### All samples- Beta diversity #######
####
### 

# Bray curtis distance matrix did converge with all samples
# Had to use a different beta diversity index, euclidean, calculated on hellinger transformed counts
data.dist <- vegdist(decostand(t(otu_table(mh_TE_PSV_transformed.ps.css)), "hell"), "euclidean") 
data.ord <- vegan::metaMDS(comm = t(data.dist), distance = "none", try = 10, trymax = 999, autotransform = F)

plt_ord_by_Mh_level <- plot_ordination(mh_TE_PSV_transformed.ps.css, data.ord, color = "Mh_level") +
  theme_bw() +
  #labs(title ="Sample type") +
  stat_ellipse(aes(fill= Mh_level), geom="polygon", alpha = 0.25) +
  theme(legend.position = "right",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "grey91", linewidth = 1.0),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", linewidth = 0.7),
        plot.title = element_text(size = 30),
        panel.border = element_rect(colour = "black", linewidth = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
plt_ord_by_Mh_level

## pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
Mh_level.permanova <- pairwise.adonis(data.dist, mh_TE_PSV_transformed.ps.css.df$Mh_level, perm = 9999, p.adjust.m = "BH")
Mh_level.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
Mh_level.disper <- betadisper(data.dist, mh_TE_PSV_transformed.ps.css.df$Mh_level)
Mh_level.permdisp <- permutest(Mh_level.disper, permutations = 9999, pairwise = T)
Mh_level.permdisp # looks like a few are significant

## pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
Concentration.permanova <- pairwise.adonis(data.dist, mh_TE_PSV_transformed.ps.css.df$Concentration, perm = 9999, p.adjust.m = "BH")
Concentration.permanova #

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
Concentration.disper <- betadisper(data.dist, mh_TE_PSV_transformed.ps.css.df$Concentration)
Concentration.permdisp <- permutest(Concentration.disper, permutations = 9999, pairwise = T)
Concentration.permdisp #

# Assuming Mh_level.disper is your betadisper object
dispersion_data <- as.data.frame(Mh_level.disper$distances)
dispersion_data$Mh_level <- mh_TE_PSV_transformed.ps.css.df$Mh_level

## Dispersion plot ####
ggplot(dispersion_data, aes(x = Mh_level, y = Mh_level.disper$distances, color = Mh_level)) + 
  geom_boxplot() + 
  labs(y = "Dispersion", x = "Mh level") +
  theme_bw() +
  scale_color_flatui(dir="-1") +
  theme(legend.position = "right",
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        legend.title = element_text(size = 10),
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

#
##
#### Subset without "Zero" ####
##
#

mhpp_PSV_nonzero.ps.css <- subset_samples(mh_TE_PSV_transformed.ps.css, Mh_level != "zero")
mhpp_PSV_nonzero.ps.css.df <- as(sample_data(mhpp_PSV_nonzero.ps.css), "data.frame")



# Had to use a different beta diversity index, euclidean, calculated on hellinger transformed counts
data_nonzero.dist <- vegdist(decostand(t(otu_table(mhpp_PSV_nonzero.ps.css)), "hell"), "euclidean") 
data_nonzero.ord <- vegan::metaMDS(comm = t(data_nonzero.dist), distance = "none", try = 10, trymax = 999, autotransform = F)

plt_ord_by_Mh_level <- plot_ordination(mhpp_PSV_nonzero.ps.css, data_nonzero.ord, color = "Mh_level") +
  theme_bw() +
  #labs(title ="Sample type") +
  #stat_ellipse(aes(fill= Mh_level), geom="polygon", alpha = 0.25) +
  theme(legend.position = "right",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "grey91", size = 1.0),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", linewidth = 0.7),
        plot.title = element_text(size = 30),
        panel.border = element_rect(colour = "black", linewidth = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
plt_ord_by_Mh_level

## pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
Mh_level_nonzero.permanova <- pairwise.adonis(data_nonzero.dist, mhpp_PSV_nonzero.ps.css.df$Mh_level, perm = 9999, p.adjust.m = "BH")
Mh_level_nonzero.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
Mh_level_nonzero.disper <- betadisper(data_nonzero.dist, mhpp_PSV_nonzero.ps.css.df$Mh_level)
Mh_level_nonzero.permdisp <- permutest(Mh_level_nonzero.disper, permutations = 9999, pairwise = T)
Mh_level_nonzero.permdisp # looks like a few are significant


# Plot richness and diversity
# Richness by sample type
ggplot(sample_data(mh_TE_PSV_transformed.ps), mapping=aes(x= Mh_level , y = Observed, fill = Mh_level , color = Mh_level )) +
  theme_bw() + 
  labs(y = "Observed PSVs") +
  geom_boxplot(alpha=0.35, position=position_dodge(width=0.8)) +
  geom_jitter(position=position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8)) +
  geom_smooth() +
  #scale_fill_manual(values = c("darkred", "darkcyan")) +
  #scale_color_manual(values = c("darkred", "darkcyan")) +
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

