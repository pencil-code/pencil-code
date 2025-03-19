#!/bin/bash -l
#SBATCH --output=test.out
#SBATCH --partition=debug  # Partition (queue) name
#SBATCH --nodes=2 # Total number of nodes 
#SBATCH --ntasks-per-node=8
#SBATCH --time=00:15:00       # Run time (d-hh:mm:ss)
#SBATCH --account=project_462000523 # Project for billing
#SBATCH --exclusive
#SBATCH --mem=0

./start.csh
./run.csh


