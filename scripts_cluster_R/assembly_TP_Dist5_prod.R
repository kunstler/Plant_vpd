# Script for cluster to run assembly with TP DIST 5
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
prod <- prod/10
print(prod)
res <- run_assembly_elev_slope(site_prod = prod, name_data_param = "output/data_slope_TP.csv",
                               disturbance_mean_interval = 5)
saveRDS(res, file = file.path("output_cluster", paste0("TP_lma_dist5_prod",prod , ".rds")))
