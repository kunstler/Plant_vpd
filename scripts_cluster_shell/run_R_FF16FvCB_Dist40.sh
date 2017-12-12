#!/bin/sh
#SBATCH --job-name=run_R_vpd
#SBATCH -p q-irstea
#SBATCH --time=00:05:00                                    
#SBATCH -A u_emgr                                          
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH --mail-type=ALL                                   
#SBATCH --mail-user=georges.kunstler@irstea.fr                                                                                                    
#SBATCH --mem=80000
  
module unload gcc
module load R/gcc
module load lapack/gcc/64/3.8.0-with-blas
module load nlopt/2.4.2                                                         

echo ${SLURM_JOB_NODELIST}
echo ${SLURM_ARRAY_TASK_ID}


# on lance le script R en precisant bien son chemin 
R CMD BATCH "--args vpd=${SLURM_ARRAY_TASK_ID}" scripts_cluster_R/assembly_FF16FvCB_Dist40_vpd.R FF16FvCB_Dist40_vpd${SLURM_ARRAY_TASK_ID}.out
