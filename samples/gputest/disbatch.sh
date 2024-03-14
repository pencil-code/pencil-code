#!/bin/bash -l
#SBATCH --output=test.out
#SBATCH --partition=dev-g  # Partition (queue) name
#SBATCH --nodes=1 # Total number of nodes 
#SBATCH --ntasks-per-node=1    # 8 MPI ranks per node, 16 total (2x8)
#SBATCH --gpus-per-node=1       # Allocate one gpu per MPI rank
#SBATCH --time=00:05:00       # Run time (d-hh:mm:ss)
#SBATCH --account=project_462000523 # Project for billing
#SBATCH --cpus-per-task=2
#
#
#export MPICH_GPU_SUPPORT_ENABLED=1
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK
#export SRUN_CPUS_PER_TASK=4
export OMP_MAX_ACTIVE_LEVELS=2
export OMP_NUM_THREADS=2
export OMP_PROC_BIND=TRUE

./start.csh
#export MPICH_GPU_SUPPORT_ENABLED=1
#srun ./src/start.x
#srun ./src/run.x
./run.csh

