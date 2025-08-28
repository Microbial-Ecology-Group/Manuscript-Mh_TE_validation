# Count summaries compared to shotgun sequencing
library(scales)
library(ggsignif)
library(cowplot)
library(dplyr)
library(scales)
library(grid)
library(gridExtra)
library(multcompView)
library(rcompanion)

seq_data.df <- read.table("BRDnoBRD_mSWEEP_results/BRDnoBRD_count_summaries.txt", header = T, row.names = 1, sep = "\t")

seq_data.df <- seq_data.df %>%
  filter(Concentration != "F") %>%
  select(-Concentration)

str(seq_data.df)
#Total_QC_reads          
# Total_QC_Deduped_reads  
# Total_Nonhost_reads     
# Total_Mh_Extracted_reads

seq_data.df$Group <- factor(seq_data.df$Group, levels = c("Shotgun","TE Double-capture (Half)"))
seq_data.df$Group # double noice

## palette
group_palette <- c("salmon1","olivedrab3")

## Wilcox testing
wilcox_list <- list()

y_vars <- c("Total_QC_reads",
            "Total_QC_Deduped_reads",
            "Total_Nonhost_reads",
            "Total_Mh_Extracted_reads")


for (y in y_vars) {
  w <- pairwise.wilcox.test(
    seq_data.df[[y]],
    seq_data.df$Group,
    p.adjust.method = "BH",
    exact = FALSE
  )
  wilcox_list[[y]] <- w
}

wilcox_list



# ### Classified Mh
# Mh_reads_plot <- ggplot(seq_data.df, aes(x=Group, y= Total_Mh_Extracted_reads, fill= Group, colour = Group)) + theme_bw() +
#   geom_point(size = 3) +
#   labs(y = expression("Reads Classified as "*italic("Mh"))) +
#   geom_boxplot(alpha=0.3, linewidth = 1.25) +
#   scale_fill_manual(values = group_palette) +
#   scale_colour_manual(values = group_palette) +
#   #scale_y_continuous(labels = label_comma()) +
#   scale_x_discrete(expand = c(0,1)) +
#   #geom_text(aes(label = SampleNameRep), position = position_nudge(x=0.5), show.legend = F,
#   #          size = 5) +
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
# Mh_reads.stats <- pairwise.wilcox.test(seq_data.df$Total_Mh_Extracted_reads, seq_data.df$Group, p.adjust.method = "BH") 
# 
# Mh_reads.Pvalues <- Mh_reads.stats$p.value
# Mh_reads.Pvalues %>%
#   fullPTable() %>%
#   multcompLetters()


# svg("figures/mh_GSV_dendrogram_rel_abun_by_BRD.svg", width = 16, height = 10)  # open SVG device
# print(plots_dendro_rel_GSV_wbar)                             # draw the plot
# dev.off() 
# 
# 
# # # This requires a directory called figures to be in your working directory
# # # All 3 lines have to be run for the file to be created.
# png("figures/mh_GSV_dendrogram_rel_abun_by_BRD.png", width = 1400, height = 800)
# plots_dendro_rel_GSV_wbar
# dev.off()


## Trying new way to make figures ####
###############################################################################
##  Helper: box-and-dot plot + Wilcoxon letters
###############################################################################
# plot_wilcox_letters <- function(df, yvar, ylab, palette,
#                                 test_var      = "Group",
#                                 alpha         = 0.05,
#                                 letter_nudge  = 0.05,
#                                 jitter_width  = 0.08,
#                                 letter_size   = 16) {
#   
#   ## ------------------------------------------------------------------
#   ## 0.  clean & symbols
#   ## ------------------------------------------------------------------
#   df <- dplyr::filter(df, !is.na(.data[[yvar]]))
#   y_sym   <- rlang::ensym(yvar)
#   grp_sym <- rlang::sym(test_var)
#   
#   df[[test_var]] <- droplevels(factor(df[[test_var]]))
#   groups <- levels(df[[test_var]])
#   k      <- length(groups)
#   
#   ## ------------------------------------------------------------------
#   ## 1.  pair-wise Wilcoxon (BH)
#   ## ------------------------------------------------------------------
#   pw <- pairwise.wilcox.test(df[[yvar]], df[[test_var]],
#                              p.adjust.method = "BH")
#   print(pw)
#   ## ------------------------------------------------------------------
#   ## 2.  rebuild FULL symmetric p-value matrix
#   ## ------------------------------------------------------------------
#   p_full <- matrix(1, k, k, dimnames = list(groups, groups))
#   
#   # pw$p.value rows are groups[2:k], columns are groups[1:(k-1)]
#   rn <- rownames(pw$p.value)
#   cn <- colnames(pw$p.value)
#   
#   for (r in rn) {
#     for (c in cn) {
#       p_full[r, c] <- p_full[c, r] <- pw$p.value[r, c]
#     }
#   }
#   diag(p_full) <- 1                   # self-comparison
#   
#   ## ------------------------------------------------------------------
#   ## 3.  compact-letter display
#   ## ------------------------------------------------------------------
#   print(p_full)
#   letters <- multcompView::multcompLetters(p_full,
#                                            threshold = alpha)$Letters
#   letters <- letters[groups]          # ensure exact order
#   print(letters)
#   ## ------------------------------------------------------------------
#   ## 4.  y-positions for the letters
#   ## ------------------------------------------------------------------
#   letter_df <- tibble::tibble(!!grp_sym := groups,
#                               Letters   = letters) |>
#     dplyr::left_join(
#       df |>
#         dplyr::group_by(!!grp_sym) |>
#         dplyr::summarise(ypos = max(.data[[yvar]], na.rm = TRUE) *
#                            (1 + letter_nudge),
#                          .groups = "drop"),
#       by = test_var)
#   
#   ## ------------------------------------------------------------------
#   ## 5.  build plot
#   ## ------------------------------------------------------------------
#   ggplot2::ggplot(df,
#                   ggplot2::aes(x = !!grp_sym, y = !!y_sym,
#                                fill = !!grp_sym, colour = !!grp_sym)) +
#     ggplot2::geom_boxplot(alpha = .30, linewidth = 1.25,
#                           outlier.shape = NA) +
#     #ggplot2::geom_point(size = 3,position = ggplot2::position_jitter(width = jitter_width)) +
#     ggplot2::geom_text(data = letter_df,
#                        ggplot2::aes(y = ypos, label = Letters),
#                        size = letter_size, vjust = 0, colour = "black") +
#     ggplot2::scale_fill_manual(values = palette) +
#     ggplot2::scale_colour_manual(values = palette) +
#     ggplot2::labs(y = ylab) +
#     ggplot2::scale_y_continuous(labels = scales::label_comma()) +
#     ggplot2::scale_x_discrete(expand = c(0, 1)) +
#     ggplot2::theme_bw() +
#     ggplot2::theme(
#       plot.margin  = ggplot2::margin(1, 1, 1, 1, "lines"),
#       panel.border = ggplot2::element_rect(colour = "black", linewidth = 2),
#       legend.position = "none",
#       #legend.title = element_blank(),
#       panel.grid.major.x = ggplot2::element_blank(),
#       axis.title.x = ggplot2::element_blank(),
#       axis.text.x  = ggplot2::element_blank(),
#       axis.ticks.x = ggplot2::element_blank(),
#       axis.title.y = ggplot2::element_text(size = 44),
#       axis.text.y  = ggplot2::element_text(size = 24, colour = "black"),
#       axis.ticks.y = ggplot2::element_line(colour = "black", linewidth = 1)
#     )
# }
###############################################################################
##  Build the four panels
###############################################################################
group_palette <- c("salmon1","olivedrab3")

seq_data.df$Log_nonhost <- log(seq_data.df$Total_Nonhost_reads)
seq_data.df$Log_mh <- log(seq_data.df$Total_Mh_Extracted_reads)
seq_data.df$Log_gsv <- log(seq_data.df$GSV_counts)

Total_Mh_Extracted_reads
GSV_counts

number_merged_QF_reads_plot <- plot_wilcox_letters(
  seq_data.df,
  yvar = "Total_QC_reads",
  ylab = "Quality-filtered merged reads",
  palette = group_palette,
  test_var = "Group"
)

deduped_merged_QF_reads_plot <- plot_wilcox_letters(
  seq_data.df,
  yvar = "Total_QC_Deduped_reads",
  ylab = "Deduplicated merged reads",
  palette = group_palette,
  test_var = "Group"
)

nonhost_reads_plot <- plot_wilcox_letters(
  seq_data.df,
  yvar = "Log_nonhost",
  ylab = "Log non-host reads",
  palette = group_palette,
  test_var = "Group"
)

Mh_reads_plot <- plot_wilcox_letters(
  seq_data.df,
  yvar = "Log_mh",
  ylab = expression("Log reads classified as "*italic("M. haemolytica")),
  palette = group_palette,
  test_var = "Group"
)

GSV_reads_plot <- plot_wilcox_letters(
  seq_data.df,
  yvar = "Log_gsv",
  ylab = expression("Log reads classified as GSVs"),
  palette = group_palette,
  test_var = "Group"
)

###############################################################################
##  Shared legend and combined figure
###############################################################################
# extract legend once
legend_plot   <- number_merged_QF_reads_plot +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(4,"lines"),
        legend.key.height = unit(6,"lines"),
        legend.text = element_text(size = 36))
shared_legend <- cowplot::get_legend(legend_plot)

# assemble the grid
reads_plot <- cowplot::plot_grid(
  number_merged_QF_reads_plot,
  deduped_merged_QF_reads_plot,
  nonhost_reads_plot,
  Mh_reads_plot,
  GSV_reads_plot,
  shared_legend,
  ncol       = 6,
  align      = "h",
  axis       = "tb",
  rel_widths = c(1, 1, 1, 1, 1.1)   # tweak as desired
)

###############################################################################
##  Display / save
###############################################################################
print(reads_plot)


png("../Writing/June_2025_draft/figures/Fig_BRDnoBRD_TEvsShotgun_counts.png", width = 3450, height = 1125)
reads_plot
dev.off()
#3450px by 1125px


combined_reads_plot <- plot_grid(optimization_reads_plot,reads_plot , nrow = 2,
                                 align      = "h", labels = "AUTO",
                                 rel_widths = c(1, 1), label_size = 40)
png("../Writing/June_2025_draft/figures/log_Fig_Combined_counts.png", width = 3450, height = 2000)
combined_reads_plot
dev.off()

## Combine both count types ####


## ───────────────────────────────────────────────────────────────
## 1.  make the four optimisation panels  (NO legend)
## ───────────────────────────────────────────────────────────────
opt_row_core <- plot_grid(
  opt_number_merged_QF_reads_plot + theme(legend.position = "none"),
  opt_deduped_merged_QF_reads_plot + theme(legend.position = "none"),
  opt_nonhost_reads_plot + theme(legend.position = "none"),
  opt_Mh_reads_plot + theme(legend.position = "none"),
  opt_GSV_reads_plot + theme(legend.position = "none"),
  ncol       = 5,
  align      = "h",
  axis       = "tb",
  labels = c("A)","B)","C)","D)","E)"),
  label_size = 38
)



## 2.  make the four “reads” panels        (NO legend)
reads_row_core <- plot_grid(
  number_merged_QF_reads_plot + theme(legend.position = "none"),
  deduped_merged_QF_reads_plot + theme(legend.position = "none"),
  nonhost_reads_plot + theme(legend.position = "none"),
  Mh_reads_plot + theme(legend.position = "none"),
  GSV_reads_plot + theme(legend.position = "none"),
  ncol       = 5,
  align      = "h",
  axis       = "tb",
  labels = c("F)","G)","H)","I)","J)"),
  label_size = 38
)

## ───────────────────────────────────────────────────────────────
## 3.  extract a legend for EACH row
## ───────────────────────────────────────────────────────────────
opt_legend   <- get_legend(opt_number_merged_QF_reads_plot   + theme(legend.position = "right",
                                                                     legend.title = element_blank(),
                                                                     legend.key.size = unit(4,"lines"),
                                                                     legend.key.height = unit(6,"lines"),
                                                                     legend.text = element_text(size = 36)))

reads_legend <- get_legend(number_merged_QF_reads_plot  + theme(legend.position = "right",
                                                                legend.title = element_blank(),
                                                                legend.key.size = unit(4,"lines"),
                                                                legend.key.height = unit(6,"lines"),
                                                                legend.text = element_text(size = 36)))

## choose one width large enough for BOTH legends
legend_width <- max(sum(opt_legend$widths), sum(reads_legend$widths))

## 1️⃣  find a common width vector
common_w <- unit.pmax(opt_legend$widths, reads_legend$widths)

## 2️⃣  give BOTH legends that exact vector
opt_legend$widths   <- common_w
reads_legend$widths <- common_w
## ───────────────────────────────────────────────────────────────
## 4.  build FULL rows: 5 plots + legend (same rel_widths)
## ───────────────────────────────────────────────────────────────
legend_rel <- 0.28          # relative width for the legend column

opt_row <- plot_grid(
  opt_row_core,  opt_legend,
  ncol       = 2,
  rel_widths = c(1.5, legend_rel)
)

reads_row <- plot_grid(
  reads_row_core,  reads_legend,
  ncol       = 2,
  rel_widths = c(1.5, legend_rel)
)

## ───────────────────────────────────────────────────────────────
## 5.  stack the two rows
## ───────────────────────────────────────────────────────────────
combined_reads_plot <- plot_grid(
  opt_row,
  reads_row,
  ncol        = 1,         # vertical stacking
  align       = "v",
  axis        = "lr",
  rel_heights = c(1, 1)
)

print(combined_reads_plot)

png("../Writing/June_2025_draft/figures/log_Fig_Combined_counts.png", width = 3450, height = 2000)
combined_reads_plot
dev.off()
