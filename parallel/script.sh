#!/bin/bash
#SBATCH --job-name=scaling
#SBATCH --chdir=.
#SBATCH --output=omp_%j.out
#SBATCH --error=omp_%j.err
#SBATCH --time=00:10:00
#SBATCH --qos=gp_debug

N=${SLURM_CPUS_PER_TASK}
echo $N
mpirun -n ${N} ./MD.exe

