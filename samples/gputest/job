#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p gputest
#SBATCH -t 00:05:00
#SBATCH --gres=gpu:k80:1
#SBATCH --ntasks-per-node=1

#module load gcc/4.9.0 cuda/7.5 openmpi/2.1.2 cmake/3.5.2
#srun --cpus-per-task=1 -N 1 -n 1 --ntasks-per-node=1 --gres=gpu:k80:1 ./run.csh

rm -f LOCK
./start.csh
touch data/jobid.dat
./run.csh
