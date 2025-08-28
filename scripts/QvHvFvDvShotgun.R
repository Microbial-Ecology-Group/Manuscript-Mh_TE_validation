######################################### 1. PREP/IMPORT ####
### SET WD
#setwd("/scratch/user/u.lp124119/20250609_MhValidation_QvHvFvDvShotgun/")

#### LOAD LIBRARIES
library(ggplot2)
library(ggsignif)
library(cowplot)
library(dplyr)
library(scales)
library(grid)
library(gridExtra)
library(multcompView)
library(rcompanion)

#### IMPORT DATA
data <- read.csv("BRDnoBRD_mSWEEP_results/20250709_MhVal_QvHvFvDvShotgun.csv", row.names = 1)
data # noice

# reorder variable levels Group to the way i want it to look
data$Group <- factor(data$Group, levels = c("Shotgun","TE Single-capture (Quarter)","TE Single-capture (Half)","TE Single-capture (Full)","TE Double-capture (Quarter)","TE Double-capture (Half)"))
data$Group # double noice

## palette
group_palette <- c("salmon1","dodgerblue4","dodgerblue2","skyblue1","darkgreen","olivedrab3")

######################################### 2. PLOTS AND STATS ####
### QF_Merged Reads Plot (w/ legend so I can extract it later)
# number_merged_QF_reads_plot_legend <- ggplot(data, aes(x=Group, y= QF_NumReads_PostMerge, fill= Group, colour = Group)) + theme_bw() +
#   geom_point(size = 3) +
#   labs(y= "Raw Reads") +
#   geom_boxplot(alpha=0.3, linewidth = 1.25) +
#   scale_fill_manual(values = group_palette) +
#   scale_colour_manual(values = group_palette) +
#   scale_y_continuous(labels = label_comma()) +
#   geom_text(aes(label = SampleNameRep), position = position_nudge(x=0.28, y= 1.5), show.legend = F,
#             size = 5) +
#   theme(plot.margin = unit(c(1,1,1,1),"lines"),
#         panel.border = element_rect(colour = "black", linewidth = 2),
#         legend.position = "right",
#         legend.title = element_blank(),
#         legend.key.size = unit(4,"lines"),
#         legend.key.height = unit(6,"lines"),
#         legend.text = element_text(size = 36),
#         panel.grid.major.x = element_blank(),
#         axis.title.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.title.y = element_text(size = 44),
#         axis.text.y = element_text(size = 24, colour = "black"),
#         axis.ticks.y = element_line(colour = "black", linewidth = 1))
# number_merged_QF_reads_plot_legend
# 
# ### QF_MergedReads Plot (no legend for actual figure)
# number_merged_QF_reads_plot <- ggplot(data, aes(x=Group, y= QF_NumReads_PostMerge, fill= Group, colour = Group)) + theme_bw() +
#   geom_point(size = 3) +
#   labs(y= "Qaulity Filtered Merged Reads") +
#   geom_boxplot(alpha=0.3, linewidth = 1.25) +
#   scale_fill_manual(values = group_palette) +
#   scale_colour_manual(values = group_palette) +
#   scale_y_continuous(labels = label_comma()) +
#   scale_x_discrete(expand = c(0,1)) +
#   geom_text(aes(label = SampleNameRep), position = position_nudge(x=0.48), show.legend = F,
#             size = 5) +
#   theme(plot.margin = unit(c(1,1,1,1),"lines"),
#         panel.border = element_rect(colour = "black", linewidth = 2),
#         legend.position = "none",
#         legend.title = element_blank(),
#         legend.key.size = unit(4,"lines"),
#         legend.key.height = unit(6,"lines"),
#         legend.text = element_text(size = 36),
#         panel.grid.major.x = element_blank(),
#         axis.title.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.title.y = element_text(size = 44),
#         axis.text.y = element_text(size = 24, colour = "black"),
#         axis.ticks.y = element_line(colour = "black", linewidth = 1)) +
#   geom_signif(colour ="black", textsize = 12, size = 0.75, vjust = 0.5, position = position_nudge(y=0.5))
# number_merged_QF_reads_plot
# 
# number_merged_QF_reads.stats <- pairwise.wilcox.test(data$QF_NumReads_PostMerge, data$Group, p.adjust.method = "BH")
# number_merged_QF_reads.Pvalues <- number_merged_QF_reads.stats$p.value
# number_merged_QF_reads.Pvalues %>%
#   fullPTable() %>%
#   multcompLetters()
# 
# ### Deduped Reads Plot 
# deduped_merged_QF_reads_plot <- ggplot(data, aes(x=Group, y= PostMerge_Dedup_Num, fill= Group, colour = Group)) + theme_bw() +
#   geom_point(size = 3) +
#   labs(y= "Deduplicated Merged Reads") +
#   geom_boxplot(alpha=0.3, linewidth = 1.25) +
#   scale_fill_manual(values = group_palette) +
#   scale_colour_manual(values = group_palette) +
#   scale_y_continuous(labels = label_comma()) +
#   scale_x_discrete(expand = c(0,1)) +
#   geom_text(aes(label = SampleNameRep), position = position_nudge(x=0.48), show.legend = F,
#             size = 5) +
#   theme(plot.margin = unit(c(1,1,1,1),"lines"),
#         panel.border = element_rect(colour = "black", linewidth = 2),
#         legend.position = "none",
#         legend.title = element_blank(),
#         legend.key.size = unit(4,"lines"),
#         legend.key.height = unit(6,"lines"),
#         legend.text = element_text(size = 36),
#         panel.grid.major.x = element_blank(),
#         axis.title.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.title.y = element_text(size = 44),
#         axis.text.y = element_text(size = 24, colour = "black"),
#         axis.ticks.y = element_line(colour = "black", linewidth = 1))
# deduped_merged_QF_reads_plot
# deduped_merged_QF.stats <- pairwise.wilcox.test(data$PostMerge_Dedup_Num, data$Group, p.adjust.method = "BH") 
# 
# deduped_merged_QF.Pvalues <- deduped_merged_QF.stats$p.value
# deduped_merged_QF.Pvalues %>%
#   fullPTable() %>%
#   multcompLetters()
# 
# 
# ### Non-host reads
# nonhost_reads_plot <- ggplot(data, aes(x=Group, y= NonHost_Num, fill= Group, colour = Group)) + theme_bw() +
#   geom_point(size = 3) +
#   labs(y= "Non-host Reads") +
#   geom_boxplot(alpha=0.3, linewidth = 1.25) +
#   scale_fill_manual(values = group_palette) +
#   scale_colour_manual(values = group_palette) +
#   scale_y_continuous(labels = label_comma()) +
#   scale_x_discrete(expand = c(0,1)) +
#   geom_text(aes(label = SampleNameRep), position = position_nudge(x=0.5), show.legend = F,
#             size = 5) +
#   theme(plot.margin = unit(c(1,1,1,1),"lines"),
#         panel.border = element_rect(colour = "black", linewidth = 2),
#         legend.position = "none",
#         legend.title = element_blank(),
#         legend.key.size = unit(4,"lines"),
#         legend.key.height = unit(6,"lines"),
#         legend.text = element_text(size = 36),
#         panel.grid.major.x = element_blank(),
#         axis.title.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.title.y = element_text(size = 44),
#         axis.text.y = element_text(size = 24, colour = "black"),
#         axis.ticks.y = element_line(colour = "black", linewidth = 1))
# nonhost_reads_plot
# nonhost_reads.stats <- pairwise.wilcox.test(data$NonHost_Num, data$Group, p.adjust.method = "BH") 
# 
# nonhost_reads.Pvalues <- nonhost_reads.stats$p.value
# nonhost_reads.Pvalues %>%
#   fullPTable() %>%
#   multcompLetters()
# 
# ### Classified Mh
# Mh_reads_plot <- ggplot(data, aes(x=Group, y= MhReads_Num, fill= Group, colour = Group)) + theme_bw() +
#   geom_point(size = 3) +
#   labs(y = expression("Reads Classified as "*italic("Mh"))) +
#   geom_boxplot(alpha=0.3, linewidth = 1.25) +
#   scale_fill_manual(values = group_palette) +
#   scale_colour_manual(values = group_palette) +
#   scale_y_continuous(labels = label_comma()) +
#   scale_x_discrete(expand = c(0,1)) +
#   geom_text(aes(label = SampleNameRep), position = position_nudge(x=0.5), show.legend = F,
#             size = 5) +
#   theme(plot.margin = unit(c(1,1,1,1),"lines"),
#         panel.border = element_rect(colour = "black", linewidth = 2),
#         legend.position = "none",
#         legend.title = element_blank(),
#         legend.key.size = unit(4,"lines"),
#         legend.key.height = unit(6,"lines"),
#         legend.text = element_text(size = 36),
#         panel.grid.major.x = element_blank(),
#         axis.title.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.title.y = element_text(size = 44),
#         axis.text.y = element_text(size = 24, colour = "black"),
#         axis.ticks.y = element_line(colour = "black", linewidth = 1))
# Mh_reads_plot
# Mh_reads.stats <- pairwise.wilcox.test(data$MhReads_Num, data$Group, p.adjust.method = "BH") 
# 
# Mh_reads.Pvalues <- Mh_reads.stats$p.value
# Mh_reads.Pvalues
# Mh_reads.Pvalues %>%
#   fullPTable() %>%
#   multcompLetters()

# 
# #### combining plots with the legend, making all purty like, and adding significance letters
# plot_legend <- get_legend(number_merged_QF_reads_plot_legend) # grab legend
# ggdraw() + draw_plot(plot_legend) # legend looks gouda
# 
# reads_plot1 <- plot_grid(number_merged_QF_reads_plot, deduped_merged_QF_reads_plot, nonhost_reads_plot, Mh_reads_plot, plot_legend, ncol = 5, rel_widths = c(1,1,0.95,0.9)) # base model combined plot
# reads_plot2 <- ggdraw() + 
#   draw_plot(reads_plot1) +
#   draw_text("a", size = 28, y= c(0.965,0.965, 0.117, 0.1), x= c(0.094, 0.3005,0.501,0.684)) + # all A
#   draw_text("b", size = 28, y=c(0.54,0.706,0.602,0.635,0.595), x=c(0.114,0.132,0.149,0.167,0.184)) + # b's for QF
#   draw_text("b", size = 28, y=c(0.172,0.21,0.18,0.193), x=c(0.320,0.3375,0.373,0.390)) + # b's for Deduped
#   draw_text("c", size = 28, y=c(0.27), x=c(0.354)) + # c for deduped
#   draw_text("b", size = 28, y=c(0.68), x=c(0.587)) + # b for nonhost
#   draw_text("ab", size = 28, y=c(0.665,0.17,0.169,0.86), x=c(0.518,0.535,0.552,0.568)) + # ab's nonhost
#   draw_text("ab", size = 28, y=c(0.163), x=c(0.7188)) + # ab for Mh
#   draw_text("abc", size = 28, y=c(0.694,0.82), x=c(0.702,0.754)) + # abcs for Mh
#   draw_text("bc", size = 28, y=c(0.161), x=c(0.736)) + # bc for Mh
#   draw_text("c", size = 28, y=c(0.68), x=c(0.771)) # c for Mh
# 
# reads_plot2 # looks great (Figure 2); dimensions to exactly match Figure 2: 3450px by 1125px
# 

# Enrique's figure ####
plot_wilcox_letters <- function(df, yvar, ylab, palette,
                                test_var      = "Group",
                                alpha         = 0.05,
                                letter_nudge  = 0.05,
                                jitter_width  = 0.08,
                                letter_size   = 14) {
  
  df <- dplyr::filter(df, !is.na(.data[[yvar]]))
  y_sym   <- rlang::ensym(yvar)
  grp_sym <- rlang::sym(test_var)
  
  df[[test_var]] <- droplevels(factor(df[[test_var]]))
  groups <- levels(df[[test_var]])
  k      <- length(groups)
  if (k < 2) stop("Need ≥ 2 groups in ", test_var, ".")
  
  ## 1) pairwise Wilcoxon (BH)
  pw <- pairwise.wilcox.test(df[[yvar]], df[[test_var]],
                             p.adjust.method = "BH")
  
  ## 2) full symmetric p-value matrix
  p_full <- matrix(1, k, k, dimnames = list(groups, groups))
  rn <- rownames(pw$p.value); cn <- colnames(pw$p.value)
  for (r in rn) for (c in cn) p_full[r, c] <- p_full[c, r] <- pw$p.value[r, c]
  diag(p_full) <- 1
  
  print(pw)
  ## 3) fix undefined cells only
  nonfin <- !is.finite(p_full)
  if (any(nonfin)) {
    warning("Pairwise Wilcoxon had undefined p-values; setting those cells to 1 (no difference).")
    p_full[nonfin] <- 1
  }
  
  ## 4) compact-letter display
  letters <- multcompView::multcompLetters(p_full, threshold = alpha)$Letters
  letters <- letters[groups]  # keep panel order
  
  ## 5) y-positions
  letter_df <- tibble::tibble(!!grp_sym := groups, Letters = letters) |>
    dplyr::left_join(
      df |>
        dplyr::group_by(!!grp_sym) |>
        dplyr::summarise(ypos = max(.data[[yvar]], na.rm = TRUE) * (1 + letter_nudge),
                         .groups = "drop"),
      by = test_var
    )
  print(letter_df)
  ## 6) plot
  ggplot2::ggplot(
    df,
    ggplot2::aes(x = !!grp_sym, y = !!y_sym, fill = !!grp_sym, colour = !!grp_sym)
  ) +
    ggplot2::geom_boxplot(alpha = .30, linewidth = 1.25, outlier.shape = NA) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_text(data = letter_df,
                       ggplot2::aes(y = ypos, label = Letters),
                       size = letter_size, vjust = 0, colour = "black") +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_colour_manual(values = palette) +
    ggplot2::labs(y = ylab) +
    ggplot2::scale_y_continuous(labels = scales::label_comma()) +
    ggplot2::scale_x_discrete(expand = c(0, 1)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.margin  = ggplot2::margin(1, 1, 1, 1, "lines"),
      panel.border = ggplot2::element_rect(colour = "black", linewidth = 2),
      legend.position = "none",
      panel.grid.major.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.text.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(size = 44),
      axis.text.y  = ggplot2::element_text(size = 24, colour = "black"),
      axis.ticks.y = ggplot2::element_line(colour = "black", linewidth = 1)
    )
}




###############################################################################
##  Build the four panels
###############################################################################
group_palette <- c("salmon1","dodgerblue4","dodgerblue2","skyblue1","darkgreen","olivedrab3")

opt_number_merged_QF_reads_plot <- plot_wilcox_letters(
  data,
  yvar = "QF_NumReads_PostMerge",
  ylab = "Quality-filtered merged reads",
  palette = group_palette,
  test_var = "Group"
)

opt_deduped_merged_QF_reads_plot <- plot_wilcox_letters(
  data,
  yvar = "PostMerge_Dedup_Num",
  ylab = "Deduplicated merged reads",
  palette = group_palette,
  test_var = "Group"
)

opt_nonhost_reads_plot <- plot_wilcox_letters(
  data,
  yvar = "NonHost_Num",
  ylab = "Non-host reads",
  palette = group_palette,
  test_var = "Group"
)

opt_Mh_reads_plot <- plot_wilcox_letters(
  data,
  yvar = "MhReads_Num",
  ylab = expression("Reads classified as "*italic("M. haemolytica")),
  palette = group_palette,
  test_var = "Group"
)

opt_GSV_reads_plot <- plot_wilcox_letters(
  data,
  yvar = "GSV_counts",
  ylab = expression("Reads classified as GSVs"),
  palette = group_palette,
  test_var = "Group"
)

###############################################################################
##  Shared legend and combined figure
###############################################################################
# extract legend once
opt_shared_legend   <- opt_number_merged_QF_reads_plot +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(4,"lines"),
        legend.key.height = unit(6,"lines"),
        legend.text = element_text(size = 36))
opt_shared_legend <- cowplot::get_legend(opt_shared_legend)

# assemble the grid
optimization_reads_plot <- cowplot::plot_grid(
  opt_number_merged_QF_reads_plot,
  opt_deduped_merged_QF_reads_plot,
  opt_nonhost_reads_plot,
  opt_Mh_reads_plot,
  opt_GSV_reads_plot,
  opt_shared_legend,
  ncol       = 6,
  align      = "h",
  axis       = "tb",
  rel_widths = c(1, 1, 1, 1, 1.1)   # tweak as desired
)

###############################################################################
##  Display / save
###############################################################################
print(optimization_reads_plot)


# png("../Writing/June_2025_draft/figures/Fig_BRDnoBRD_TEvsShotgun_counts.png", width = 3450, height = 1125)
# reads_plot
# dev.off()




png("../Writing/June_2025_draft/figures/Fig_QvHvFvDvShotgun_counts.png", width = 3450, height = 1125)
optimization_reads_plot
dev.off()

######################################### 3. GETTING SOME SUMMARY STATS FOR THE TEXT ####
## sample sizes
count(data, Technique) # 39 TE, 5 shotgun
count(data, Group) # 5-9 by group

#### split into TE and shotgun for raw read cluster and Merged/QF read numbers
TE_Data <- data[which(data$Technique=="TE"),]
Shotgun_Data <- data[which(data$Technique=="Shotgun"),]
## TE numbers
# PE Clusters
mean(TE_Data$PE_Clusters) # 71,846,997
min(TE_Data$PE_Clusters) # 1,170,835
max(TE_Data$PE_Clusters) # 153,793,168

# QF & Merged Reads
mean(TE_Data$QF_NumReads_PostMerge) # 99,690,079

## Shotgun numbies
# PE Clusters
mean(Shotgun_Data$PE_Clusters) # 131,233,940
min(Shotgun_Data$PE_Clusters) # 120,365,961
max(Shotgun_Data$PE_Clusters) # 141,073,067

## Merged & QF
mean(Shotgun_Data$QF_NumReads_PostMerge) # 99,690,079

######## Mh READS
mean(Shotgun_Data$MhReads_Num) # 1,102.6 in shotgun
sd(Shotgun_Data$MhReads_Num)/(sqrt(5)) # 804.72 SEM
mean(Shotgun_Data$MhReads_Num)/mean(Shotgun_Data$QF_NumReads_PostMerge)*100 # 0.0005 % of QF/merged reads

### split into TE modified libraries
SC_Full <- TE_Data[which(TE_Data$Group=="TE Single-capture (Full)"),]
SC_Half <- TE_Data[which(TE_Data$Group=="TE Single-capture (Half)"),]
SC_Quarter <- TE_Data[which(TE_Data$Group=="TE Single-capture (Quarter)"),]
DC_Half <- TE_Data[which(TE_Data$Group=="TE Double-capture (Half)"),]
DC_Qaurter<- TE_Data[which(TE_Data$Group=="TE Double-capture (Quarter)"),]

mean(SC_Full$MhReads_Num) # 123,541.00
sd(SC_Full$MhReads_Num)/sqrt(length(SC_Full$MhReads_Num)) # 71,655.40

mean(SC_Half$MhReads_Num) # 57,009.25
sd(SC_Half$MhReads_Num)/sqrt(length(SC_Half$MhReads_Num)) # 47,010.45

mean(SC_Quarter$MhReads_Num) # 101,318.00
sd(SC_Quarter$MhReads_Num)/sqrt(length(SC_Quarter$MhReads_Num)) # 64,214.77

mean(DC_Half$MhReads_Num) # 269,888.80
sd(DC_Half$MhReads_Num)/sqrt(length(DC_Half$MhReads_Num)) #47,090.53

mean(DC_Half$MhReads_Num)/mean(DC_Half$QF_NumReads_PostMerge)*100 # 0.2297 % of QF/merged reads
(mean(DC_Half$MhReads_Num)/mean(DC_Half$QF_NumReads_PostMerge)*100)/(mean(Shotgun_Data$MhReads_Num)/mean(Shotgun_Data$QF_NumReads_PostMerge)*100) # 455.02-fold increase v. shotgun

mean(DC_Qaurter$MhReads_Num) # 205,801.90
sd(DC_Qaurter$MhReads_Num)/sqrt(length(DC_Qaurter$MhReads_Num)) #55,444.48

mean(DC_Qaurter$MhReads_Num)/mean(DC_Qaurter$QF_NumReads_PostMerge)*100 # 0.2346 % of QF/merged reads
(mean(DC_Qaurter$MhReads_Num)/mean(DC_Qaurter$QF_NumReads_PostMerge)*100)/(mean(Shotgun_Data$MhReads_Num)/mean(Shotgun_Data$QF_NumReads_PostMerge)*100) # 463.72-fold increase v. shotgun



