#!/bin/bash -l

#SBATCH --account=project_462000295

#SBATCH -N 1

#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-node=1
#SBATCH -p dev-g

#SBATCH -t 00:05:00
#SBATCH --cpus-per-task=1


module purge
module load CrayEnv
module load PrgEnv-cray
module load craype-accel-amd-gfx90a
module load rocm
module load buildtools
module load cray-python
#./start.csh
./run.csh
