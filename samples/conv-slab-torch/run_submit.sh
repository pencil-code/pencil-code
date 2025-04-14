#!/bin/bash
#SBATCH --job-name=test_torchfort
#SBATCH --account=project_2000403
#SBATCH --partition=gputest
#SBATCH --time=00:03:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --gres=gpu:v100:1
##SBATCH --mem=32G
#SBATCH -o slurm-%x_%J.out
####SBATCH --mail-type=ALL
#SBATCH --output=train.out

#TP: for god knows for what reason a blank data folder gets generated on top of the correct data folder after running
#TP: and even more weirdly if one deletes the blank folder the correct one reappears!
apptainer exec --nv -B $PENCIL_HOME:$PENCIL_HOME -B scripts:/opt/scripts -B /scratch/project_2000403/toukopur/data:$PWD/data -B /scratch/project_2000403/toukopur/data:/data torchfort_bisonflex.sif bash rm -rf data

apptainer exec --nv -B $PENCIL_HOME:$PENCIL_HOME -B scripts:/opt/scripts -B /scratch/project_2000403/toukopur/data:$PWD/data -B /scratch/project_2000403/toukopur/data:/data torchfort_bisonflex.sif bash /opt/scripts/run.sh

