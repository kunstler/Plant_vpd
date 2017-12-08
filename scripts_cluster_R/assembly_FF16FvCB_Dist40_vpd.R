# Script for cluster to run assembly with FF16FvCB DIST 40
library(plant)
library(plant.assembly)
library(loggr)

source("R/assembly.R")

vpd <- commandArgs(trailingOnly=TRUE)/10
res <- run_assembly_FvCB_narea_lma(vpd = vpd, disturbance_mean_interval = 40)
saveRDS(res, file = file.path("output_cluster", paste0("FvCB_lma_Narea_dist5_vpd", , ".rds")))

