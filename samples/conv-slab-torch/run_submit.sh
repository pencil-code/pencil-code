#!/bin/bash
#SBATCH --job-name=test_torchfort
#SBATCH --account=ituomine
#SBATCH --partition=gputest
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:a100:4
#SBATCH --mem=32G
#SBATCH -o slurm-%x_%J.out
####SBATCH --mail-type=ALL
#SBATCH --output=train.out

apptainer exec --nv -B /users/mreinhar/pencil-code:/users/mreinhar/pencil-code -B scripts:/opt/scripts -B /scratch/project_2001062/mreinhar/pencil-code/samples/conv-slab-torch/data:/scratch/project_2001062/mreinhar/pencil-code/samples/conv-slab-torch/data torchfort_0.2.0.sif bash /opt/scripts/run.sh
