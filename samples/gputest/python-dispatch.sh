#!/bin/bash -l
#SBATCH --output=test.out
#SBATCH --partition=dev-g  # Partition (queue) name
#SBATCH --nodes=2 # Total number of nodes 
#SBATCH --ntasks-per-node=8
#SBATCH --gpus-per-node=8       # Allocate one gpu per MPI rank
#SBATCH --time=00:05:00       # Run time (d-hh:mm:ss)
#SBATCH --account=project_462000523 # Project for billing
#SBATCH --cpus-per-task=7
#SBATCH --exclusive
#SBATCH --mem=0

module load cray-python

export OMP_NUM_THREADS=7
export OMP_PROC_BIND=close,spread
export OMP_MAX_ACTIVE_LEVELS=2
export OMP_WAIT_POLICY=PASSIVE
./start.csh
export MPICH_GPU_SUPPORT_ENABLED=1
export LD_PRELOAD=./src/libPC.so
#srun ./src/run.x
#./run.csh
export ODOP_PATH=~/monitoring/odop/odop/odop_obs/
source ~/monitoring/venv/bin/active
srun python call.py



