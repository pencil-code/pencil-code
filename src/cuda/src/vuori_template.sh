#!/bin/bash
#SBATCH -J astaroth_test 
#SBATCH -N 1 
#SBATCH -n 1
#SBATCH -t 00:15:00
#SBATCH --mail-type=ALL
#SBATCH -o output_astaroth_test_%j.out
#SBATCH -e output_astaroth_test_%j.err
#SBATCH -p gpu
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1

./runme 
