#!/bin/bash -l
#SBATCH --output=test.out
#SBATCH --partition=gpu  # Partition (queue) name
#SBATCH --nodes=4 # Total number of nodes 
#SBATCH --ntasks-per-node=4
#SBATCH --time=00:15:00       # Run time (d-hh:mm:ss)
#SBATCH --account=naiss2023-3-23
#SBATCH --cpus-per-task=16
#SBATCH --exclusive

#module load PDC/23.03 rocm/5.0.2 craype-accel-amd-gfx90a cmake/3.20.1
source src/.moduleinfo

export OMP_NUM_THREADS=16
export OMP_PROC_BIND=true
export OMP_MAX_ACTIVE_LEVELS=2
export OMP_WAIT_POLICY=PASSIVE
./start.csh
export MPICH_GPU_SUPPORT_ENABLED=1
#srun ./src/run.x
./run.csh


