#!/bin/bash
#SBATCH --job-name=one_d
#SBATCH --output=res_omp.txt
#SBATCH --error=err_omp.txt
##SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=haibo.ni02@gmail.com
#SBATCH --mail-type=ALL
##SBATCH --begin=18:30
##SBATCH --time=3-0  # set time limit
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# python run_AF_1Hz.py 3


# sh get_all.sh

sh run.sh