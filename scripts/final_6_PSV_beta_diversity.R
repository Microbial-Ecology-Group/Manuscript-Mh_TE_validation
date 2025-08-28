# PCV beta diversity
library(pairwiseAdonis)

#############################################################################################
####################################   CSS TRANSFORM    ###################################

sample_data(mh_TE_PSV_transformed.ps)$Mh_level <- factor(sample_data(mh_TE_PSV_transformed.ps)$Mh_level , levels = c("zero","low","medium","high","highest"))

mh_TE_PSV_transformed.ps <- prune_taxa(taxa_sums(mh_TE_PSV_transformed.ps) > 0, mh_TE_PSV_transformed.ps)
any(taxa_sums(mh_TE_PSV_transformed.ps)==0) # QUADRUPLE CHECKING - nope good.

#mh_TE_PSV_transformed.ps.css <- phyloseq_transform_css(mh_TE_PSV_transformed.ps, log = F)
mh_TE_PSV_transformed.ps.df <- as(sample_data(mh_TE_PSV_transformed.ps), "data.frame")


###
####
##### All samples- Beta diversity #######
####
### 

# Bray curtis distance matrix did converge with all samples
# Had to use a different beta diversity index, euclidean, calculated on hellinger transformed counts
data.df <- as(mh_TE_PSV_transformed.ps@sam_data,"data.frame") # make dataframe from metadata
data.dist <- vegdist(t(otu_table(mh_TE_PSV_transformed.ps)), "jaccard") 
data.ord <- vegan::metaMDS(comm = t(data.dist), distance = "none", try = 10, trymax = 999, autotransform = F)

plt_ord_by_Mh_level <- plot_ordination(mh_TE_PSV_transformed.ps, data.ord, color = "Mh_level") +
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

# Adonis2 ####
# 1.  Make a metadata frame that contains ONLY the grouping variable
meta  <- mh_TE_PSV_transformed.ps.df %>% 
  select(Mh_level) %>%                 # keep the column you’ll test
  mutate(Mh_level = factor(Mh_level))  # be sure it’s a factor

# 2.  Re-order rows so they line up with the distance matrix
meta  <- meta[ match(labels(data.dist), rownames(meta)), , drop = FALSE ]

# 3.  Run PERMANOVA with 9 999 permutations
adonis2_res <- adonis2(
  data.dist ~ Mh_level,          # distance object on the left of “~”
  data        = meta,            # one-row-per-sample data frame
  permutations = 9999,
  by          = "terms"          # gives a marginal test for each term
)

# BRD adonis2

meta  <- mh_TE_PSV_transformed.ps.df %>% 
  select(BRD) %>%                 # keep the column you’ll test
  mutate(BRD = factor(BRD))  # be sure it’s a factor

# 2.  Re-order rows so they line up with the distance matrix
meta  <- meta[ match(labels(data.dist), rownames(meta)), , drop = FALSE ]

# 3.  Run PERMANOVA with 9 999 permutations
adonis2_res <- adonis2(
  data.dist ~ BRD,          # distance object on the left of “~”
  data        = meta,            # one-row-per-sample data frame
  permutations = 9999,
  by          = "terms"          # gives a marginal test for each term
)

print(adonis2_res)
## pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
Mh_level.permanova <- pairwise.adonis(data.dist, mh_TE_PSV_transformed.ps.df$Mh_level, perm = 9999, p.adjust.m = "BH")
Mh_level.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
Mh_level.disper <- betadisper(data.dist, mh_TE_PSV_transformed.ps.df$Mh_level)
Mh_level.permdisp <- permutest(Mh_level.disper, permutations = 9999, pairwise = T)
Mh_level.permdisp # looks like a few are significant

## pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
BRD.permanova <- pairwise.adonis(data.dist, mh_TE_PSV_transformed.ps.df$BRD, perm = 9999, p.adjust.m = "BH")
BRD.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
BRD.disper <- betadisper(data.dist, mh_TE_PSV_transformed.ps.df$BRD)
BRD.permdisp <- permutest(BRD.disper, permutations = 9999, pairwise = T)
BRD.permdisp # looks like a few are significant




## pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
Concentration.permanova <- pairwise.adonis(data.dist, mh_TE_PSV_transformed.ps.df$Concentration, perm = 9999, p.adjust.m = "BH")
Concentration.permanova #

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
Concentration.disper <- betadisper(data.dist, mh_TE_PSV_transformed.ps.df$Concentration)
Concentration.permdisp <- permutest(Concentration.disper, permutations = 9999, pairwise = T)
Concentration.permdisp #

# Assuming Mh_level.disper is your betadisper object
dispersion_data <- as.data.frame(Mh_level.disper$distances)
dispersion_data$Mh_level <- mh_TE_PSV_transformed.ps.df$Mh_level

## Dispersion plot ####
ggplot(dispersion_data, aes(x = Mh_level, y = Mh_level.disper$distances, color = Mh_level)) + 
  geom_boxplot() + 
  labs(y = "Dispersion", x = "Mh level") +
  theme_bw() +
  scale_color_flatui() +
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

data_nonzero.dist <- vegdist(t(otu_table(mhpp_PSV_nonzero.ps.css)), "jaccard") 

data_nonzero.ord <- vegan::metaMDS(comm = t(data_nonzero.dist), distance = "none", try = 10, trymax = 999, autotransform = F)

plt_ord_by_Mh_level <- plot_ordination(mhpp_PSV_nonzero.ps.css, data_nonzero.ord, color = "Mh_level") +
  theme_bw() +
  #labs(title ="Sample type") +
  stat_ellipse(aes(fill= Mh_level), geom="polygon", alpha = 0.25) +
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




## pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
BRD_nonzero.permanova <- pairwise.adonis(data_nonzero.dist, mhpp_PSV_nonzero.ps.css.df$BRD, perm = 9999, p.adjust.m = "BH")
BRD_nonzero.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
BRD_nonzero.disper <- betadisper(data_nonzero.dist, mhpp_PSV_nonzero.ps.css.df$BRD)
BRD_nonzero.permdisp <- permutest(BRD_nonzero.disper, permutations = 9999, pairwise = T)
BRD_nonzero.permdisp # looks like a few are significant


## Final Figure #####


mh_level_colors <- c("zero" = "yellow", 
                     "low" = "light gray", 
                     "medium" = "#6497B1", 
                     "high" = "#005B96", 
                     "highest" = "#03396C")

data.df <- as(mh_TE_PSV_transformed.ps@sam_data,"data.frame") # make dataframe from metadata
data.dist <- vegdist(t(otu_table(mh_TE_PSV_transformed.ps)), "jaccard") 
data.ord <- vegan::metaMDS(comm = t(data.dist), distance = "none", try = 10, trymax = 999, autotransform = F)



data.ord.plot <- ordiplot(data.ord$points)
data.ord.plot.scrs <- scores(data.ord.plot, display = "sites")
data.ord.plot.scrs <- cbind(as.data.frame(data.ord.plot.scrs), Mh_level= data.df$Mh_level)
data.ord.plot.cent <- aggregate(cbind(MDS1,MDS2) ~ Mh_level, data = data.ord.plot.scrs, FUN = mean)
data.ord.plot.segs <- merge(data.ord.plot.scrs, setNames(data.ord.plot.cent, c("Mh_level","cMDS1","cMDS2")), by = 'Mh_level', sort = F)

## plot it
data.ord.ggplot <- ggplot(data.ord.plot.segs, aes(fill = Mh_level, colour = Mh_level)) + theme_bw() +
  labs(x="NMDS1", y="NMDS2", fill = "Mh level",  # Add legend title for fill aesthetic
       colour = "Mh level") + # Add legend title for colour aesthetic) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  #scale_x_continuous(limits = xlims) + scale_y_continuous(limits = ylims) +
  geom_point(aes(x=MDS1,y=MDS2, colour= Mh_level), alpha = 0.5, shape = 18, size = 5) +
  stat_ellipse(geom = "polygon", aes(x=MDS1,y=MDS2), alpha = c(0.05), lty=2, level = 0.95, linewidth = 0.2) +
  geom_segment(aes(x=MDS1,y=MDS2, xend=cMDS1, yend=cMDS2), alpha = 0.1, linewidth = 1) +
  geom_point(aes(x=cMDS1, y=cMDS2), size = 18, shape = 18) +
  geom_text(aes(x=cMDS1, y=cMDS2, label = Mh_level), colour = "white", size = 5) +
  scale_fill_manual(values = mh_level_colors) +
  scale_colour_manual(values = mh_level_colors) +
  theme(legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 12, colour = "black"))

data.ord.ggplot

# Specify the output file and resolution
png("figures/PSV_beta_diversity_by_Mh_level.png", width = 2000, height = 1500, res = 150)
data.ord.ggplot
# Turn off the device to save the file
dev.off()


svg("figures/PSV_beta_diversity_by_Mh_level.svg", width = 16, height = 10)  # open SVG device
print(data.ord.ggplot)                             # draw the plot
dev.off() 
