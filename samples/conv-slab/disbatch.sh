#!/bin/bash -l

#SBATCH --account=ituomine

#SBATCH -N 1

#SBATCH -n 2
#SBATCH --ntasks-per-node=2
#SBATCH -p test

#SBATCH -t 00:05:00

export OMP_NUM_THREADS=1

module load \

gcc/11.3.0 \

intel-oneapi-mkl/2022.1.0 \

openmpi/4.1.4 \

csc-tools \

StdEnv \

cuda/11.7.0
rm ./data/time_series.dat
./start.csh
touch data/jobid.dat
./run.csh


