# Load libraries
library(phyloseq)
library(microViz)
library(dplyr)
library(ggplot2)
library(data.table)
library(tidyr)
#library(palmerpenguins)
library(ggthemes)
library(forcats)
library(vegan)
library(microbiome)
library(metagenomeSeq)
library(metagMisc)
library(ggsci)
library(pairwise)
#library(devtools)
library(pairwiseAdonis)
library(ggtree)
#library(speedyseq)
#library(seqateurs)
#library(ggstatsplot)
#library(lme4) # for mixed effects regression models
#library(oddsratio)
#library(epiDisplay)
#library(jsmodule)

# Phyloseq has issue with the phyloseq_to_metagenomeSeq() function. I made the change suggested here:
# https://github.com/joey711/phyloseq/issues/1118

phyloseq_to_metagenomeSeq = function(physeq, ...){
  # Enforce orientation. Samples are columns
  if( !taxa_are_rows(physeq) ){ physeq <- t(physeq)}
  # Coerce count data to vanilla matrix of integers
  countData = round(as(otu_table(physeq), "matrix"), digits=0)
  # Create sample annotation if possible
  if(!is.null(sample_data(physeq,FALSE))){
    ADF = AnnotatedDataFrame(data.frame(sample_data(physeq)))  
  } else { 
    ADF = NULL 
  }
  # Create taxa annotation if possible
  if(!is.null(tax_table(physeq,FALSE))){
    TDF = AnnotatedDataFrame(data.frame(OTUname = taxa_names(physeq),
                                        physeq@tax_table@.Data,row.names = taxa_names(physeq)))
  } else {
    TDF = AnnotatedDataFrame(data.frame(OTUname = taxa_names(physeq),
                                        row.names = taxa_names(physeq)))
  }
  # Create MRexperiment
  if(requireNamespace("metagenomeSeq")){
    mrobj = metagenomeSeq::newMRexperiment(counts = countData, phenoData = ADF, featureData = TDF,...)
    # Calculate normalization factor
    if (sum(colSums(countData > 0) > 1) < ncol(countData)) {
      p = suppressMessages(metagenomeSeq::cumNormStat(mrobj))
    }
    else {
      p = suppressMessages(metagenomeSeq::cumNormStatFast(mrobj))
    }
    mrobj = metagenomeSeq::cumNorm(mrobj, p = p)
    return(mrobj)
  }
}

#setDTthreads(percent=60)
#getDTthreads(verbose = TRUE)
#Cstack_info()

# install.packages("remotes")
# remotes::install_github("mikemc/speedyseq")
#remotes::install_github("alexpiper/seqateurs")