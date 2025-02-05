#!/bin/bash
#SBATCH --job-name=test_torchfort
#SBATCH --account=ituomine
#SBATCH --partition=gputest
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:v100:1
#SBATCH --mem=32G
#SBATCH -o slurm-%x_%J.out
####SBATCH --mail-type=ALL

sleep 150000

##apptainer shell --nv -B pencil-code:$PENCIL_HOME -B scripts:scripts torchfort_0.2.0.sif
