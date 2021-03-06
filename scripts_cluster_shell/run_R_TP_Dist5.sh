#!/bin/sh
#SBATCH --job-name=run_R_Dist5_TP_prod%j
#SBATCH -p q-irstea
#SBATCH --time=95:00:00                                    
#SBATCH -A u_emgr                                          
#SBATCH -o run_R_Dist5_TP_prod%j.out
#SBATCH -e run_R_Dist5_TP_prod%j.err
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
R CMD BATCH "--args prod=${SLURM_ARRAY_TASK_ID}" scripts_cluster_R/assembly_TP_Dist5_prod.R TP_Dist5_prod${SLURM_ARRAY_TASK_ID}.out
