# Script for cluster to run assembly with FF16FvCB Tleaf DIST 5
library(plant)
library(plant.assembly)
library(loggr)

args <- commandArgs(trailingOnly=TRUE)
if(length(args)==0){
    print("No arguments supplied.")
}else{
    for(i in 1:length(args)){
         eval(parse(text=args[[i]]))
    }
}

source("R/assembly.R")
dir.create("output_cluster", showWarnings = FALSE)
vpd <- vpd/10
print( vpd)
res <- run_assembly_FvCB_narea_lma_Narea_LTR(vpd = vpd,
                                         name_param_NareaLTR = "output/data_slope_narea.csv",
                                         disturbance_mean_interval = 5)
saveRDS(res, file = file.path("output_cluster", paste0("FvCB_lma_Narea_NvLTR_dist5_vpd",vpd , ".rds")))
