#!/bin/bash -l
#SBATCH --output=test.out
#SBATCH --partition=gpu  # Partition (queue) name
#SBATCH --nodes=2 # Total number of nodes 
#SBATCH --ntasks-per-node=8
#SBATCH --time=00:15:00       # Run time (d-hh:mm:ss)
#SBATCH --account=naiss2023-3-23
#SBATCH --cpus-per-task=7
#SBATCH --exclusive

#module load rocm/5.0.2 craype-accel-amd-gfx90a cmake/3.20.1
#git submodule update --remote --merge

source src/.moduleinfo

export OMP_NUM_THREADS=7
export OMP_PROC_BIND=true
export OMP_MAX_ACTIVE_LEVELS=2
export OMP_WAIT_POLICY=PASSIVE
./start.csh
export MPICH_GPU_SUPPORT_ENABLED=1
#srun ./src/run.x
./run.csh


