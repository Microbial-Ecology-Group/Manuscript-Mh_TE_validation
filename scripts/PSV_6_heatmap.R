# MicroViz annotated heatmaps
library(microViz)

mhpp_PSV_nonzero.ps <- subset_samples(mhpp_PSV.ps, Mh_level != "zero")
mhpp_PSV_nonzero.ps <- prune_taxa(taxa_sums(mhpp_PSV_nonzero.ps) > 0, mhpp_PSV_nonzero.ps)
any(taxa_sums(mhpp_PSV_nonzero.ps)==0) # QUADRUPLE CHECKING - nope good.

mhpp_PSV_nonzero.ps <- prune_taxa(taxa_sums(mhpp_PSV_nonzero.ps) > 400000, mhpp_PSV_nonzero.ps)


cols1 <- distinct_palette(n = 2, add = NA)
names(cols1) <- unique(sample_data(mhpp_PSV_nonzero.ps)$Concentration)


cols2 <- distinct_palette(n = 4, add = NA)
names(cols2) <- unique(sample_data(mhpp_PSV_nonzero.ps)$Mh_level)

mhpp_PSV_nonzero.ps %>%
  tax_transform("compositional") %>%
  comp_heatmap(
    tax_anno = taxAnnotation(
      Prev. = anno_tax_prev(bar_width = 0.3, size = grid::unit(1, "cm"))
    ),
    sample_anno = sampleAnnotation(
      State1 = anno_sample_cat("Concentration", col = cols1),
      State2 = anno_sample_cat("Mh_level", col = cols2)
    )
  )
