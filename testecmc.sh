#!/bin/bash

#SBATCH --job-name=ecmc_test
#SBATCH --output=%x.o
#SBATCH --time=00:20:00
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --partition=cpu_short

#Source necessary modules
source modules_load.sh

# Run MPI script
srun build/gauge_ecmc_cb inputs/ecmc_cb.txt
