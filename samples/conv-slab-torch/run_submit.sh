#!/bin/bash
#SBATCH --job-name=test_torchfort
#SBATCH --account=project_2000403
#SBATCH --partition=gpu
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --gres=gpu:v100:1
##SBATCH --mem=32G
#SBATCH -o slurm-%x_%J.out
####SBATCH --mail-type=ALL
#SBATCH --output=train.out


rm /scratch/project_2000403/sgiridha/data_conv-slab-torch/data/training/stationary.pt

#TP: for god knows for what reason a blank data folder gets generated on top of the correct data folder after running
#TP: and even more weirdly if one deletes the blank folder the correct one reappears!
#apptainer exec --nv -B $PENCIL_HOME:$PENCIL_HOME -B scripts:/opt/scripts -B /scratch/project_2000403/sgiridha/data:$PWD/data -B /scratch/project_2000403/sgiridha/data_conv-slab-torch:/data torchfort_bisonflex.sif bash rm -rf data

#bash clean.sh

#apptainer exec --nv -B $PENCIL_HOME:$PENCIL_HOME -B scripts:/opt/scripts -B /scratch/project_2000403/toukopur/data:$PWD/data -B /scratch/project_2000403/toukopur/data:/data torchfort_bisonflex.sif bash
#apptainer exec --nv -B $PENCIL_HOME:$PENCIL_HOME -B scripts:/opt/scripts -B /scratch/project_2000403/sgiridha/data:$PWD/data -B /scratch/project_2000403/sgiridha/data:/data torchfort_bisonflex.sif bash /opt/scripts/run.sh



apptainer exec --nv -B $PENCIL_HOME:$PENCIL_HOME -B scripts:/opt/scripts -B /scratch/project_2000403/sgiridha/data_conv-slab-torch/data:$PENCIL_HOME/samples/conv-slab-torch/data -B /scratch/project_2000403/sgiridha/data_conv-slab-torch/data:/data -B /scratch/project_2000403/sgiridha/snapshots:/$PENCIL_HOME/samples/conv-slab-torch/snapshots -B /scratch/project_2000403/sgiridha/snapshots:/snapshots torchfort_bisonflex.sif bash rm -rf data bash rm -rf snapshots


apptainer exec --nv -B $PENCIL_HOME:$PENCIL_HOME -B scripts:/opt/scripts -B /scratch/project_2000403/sgiridha/data_conv-slab-torch/data:$PENCIL_HOME/samples/conv-slab-torch/data -B /scratch/project_2000403/sgiridha/data_conv-slab-torch/data:/data -B /scratch/project_2000403/sgiridha/snapshots:/$PENCIL_HOME/samples/conv-slab-torch/snapshots -B /scratch/project_2000403/sgiridha/snapshots:/snapshots torchfort_bisonflex.sif bash



apptainer exec --nv -B $PENCIL_HOME:$PENCIL_HOME -B /appl:/appl -B scripts:/opt/scripts -B /scratch/project_2000403/sgiridha/data_conv-slab-torch/data:$PENCIL_HOME/samples/conv-slab-torch/data -B /scratch/project_2000403/sgiridha/data_conv-slab-torch/data:/data -B /scratch/project_2000403/sgiridha/snapshots:/$PENCIL_HOME/samples/conv-slab-torch/snapshots -B /scratch/project_2000403/sgiridha/snapshots:/snapshots torchfort_bisonflex.sif bash /opt/scripts/run.sh
