# Staging script for mh validation study

source("scripts/load_libraries.R")

# Count summary code
source("scripts/final_1_load_PSV_data.R")

source("scripts/final_3")


# Rel abundance of PSVs
source("scripts/final_4_PSV_prevalence_dotplot.R")



# Off target benchmarking ####

# this script is for viewing raw results from benchmarking and determine trends in FP rates
source("scripts/off_target_FP_by_GSV.R")

# THis is used to parse the various off target benchmarking results with different filters
source("scripts/off_target_FP_by_Conf.R")