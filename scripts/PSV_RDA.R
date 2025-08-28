# RDA
# RDA


cleanset <- data.frame(sample_data(mh_TE_PSV_transformed.ps))
cleanset <- cleanset[, c("Concentration","Pen","Mh_level","BRD","Nasal.Mh.Culture","DGNP.Mh.Culture",
                         "Proct.Mh.Culture","X16S_Mean_RA_per","X16S_RA_nasal_per",
                         "NS.Mh.Ct","DGNP.Mh.Ct","MALDI.Genotype")]
#"DGNP.Mh.PCR","NS.Mh.PCR","PS.Mh.PCR","PS.Mh.Ct","Order",

############################################### RDA MODEL ###################################################

#### CSS transformation of counts
mh_TE_PSV_transformed.ps.css <- phyloseq_transform_css(mh_TE_PSV_transformed.ps, log = F)

norm.mt<-as(otu_table(mh_TE_PSV_transformed.ps.css),"matrix")
### This matrix has rows as otus and columns as samples (all our work before used columns as otus)
## we first remove all rows (otus) that have count zero by summing per row and keeping those rows with
## sum greater than zero
norm.rs = apply(norm.mt,1,sum)
norm.mt = norm.mt[norm.rs>0,]
## We then transform the normalized matrix so that rows are samples and columns are otus (for vegan to work correctly)
norm.mt = t(norm.mt)
#norm.mt<-norm.mt[cleanset_samples,]
#Using function decostand() in vegan to transform using hellinger standardizations
hell.norm.mt <- decostand(norm.mt,"hell")
### We define the full model that includes all possible variables in the metadata file
mod1 <- rda(hell.norm.mt ~ ., data =  cleanset)
anova(mod1, by = "term", perm = 1000)

### We define the null model that only includes an intercept (1)
mod0 <- rda(hell.norm.mt ~ 1, data =  cleanset)
### We use the step function to find the best model by removing one term at a time.
### We use permutation test (perMANOVA) to do so
mod <- step(mod0, scope = formula(mod1), test = "perm", direction="both")
anova(mod, by = "term", perm = 1000)
anova(mod,perm = 1000)

# > anova(mod, by = "term", perm = 1000)
# Permutation test for rda under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# Model: rda(formula = hell.norm.mt ~ DGNP.Mh.Ct + NS.Mh.Ct + Pen + X16S_RA_nasal_per + Concentration, data = cleanset)
# Df  Variance      F Pr(>F)  
# DGNP.Mh.Ct         1 0.0067237 6.4781  0.077 .
# NS.Mh.Ct           1 0.0059802 5.7618  0.041 *
#   Pen                1 0.0035836 3.4527  0.095 .
# X16S_RA_nasal_per 11 0.0199878 1.7507  0.268  
# Concentration      1 0.0018989 1.8295  0.031 *
#   Residual           5 0.0051895                
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# > anova(mod,perm = 1000)
# Permutation test for rda under reduced model
# Permutation: free
# Number of permutations: 999
# 
# Model: rda(formula = hell.norm.mt ~ DGNP.Mh.Ct + NS.Mh.Ct + Pen + X16S_RA_nasal_per + Concentration, data = cleanset)
# Df Variance     F Pr(>F)
# Model    15 0.038174 2.452  0.158
# Residual  5 0.005190


