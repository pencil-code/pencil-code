#!/bin/sh
#SBATCH --job-name=vortex_braid
#SBATCH --ntasks=64
#SBATCH --nodes=2
#SBATCH --partition=staff
#SBATCH --output=run_%j.out
#SBATCH --error=run_%j.err
#SBATCH --time=14:00:00
#SBATCH --exclusive

echo $SLURM_JOB_NUM_NODES
echo $SLURM_NTASKS
echo $SLURM_TASKS_PER_NODE

#mpirun -np 8 src/start.x
#mpirun -np 8 src/run.x
#start.csh
run.csh
