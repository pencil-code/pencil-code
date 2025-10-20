#!/bin/bash
#SBATCH --job-name=test_torchfort
#SBATCH --account=project_2000403
#SBATCH --partition=gputest
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --gres=gpu:v100:1
##SBATCH --mem=64G
#SBATCH -o slurm-%x_%J.out
####SBATCH --mail-type=ALL
#SBATCH --output=train.out


#rm /scratch/project_2000403/sgiridha/data_conv-slab-torch/data/training/stationary.pt

#TP: for god knows for what reason a blank data folder gets generated on top of the correct data folder after running
#TP: and even more weirdly if one deletes the blank folder the correct one reappears!
#apptainer exec --nv -B $PENCIL_HOME:$PENCIL_HOME -B scripts:/opt/scripts -B /scratch/project_2000403/sgiridha/data:$PWD/data -B /scratch/project_2000403/sgiridha/data_conv-slab-torch:/data torchfort_bisonflex.sif bash rm -rf data

#bash clean.sh

#apptainer exec --nv -B $PENCIL_HOME:$PENCIL_HOME -B scripts:/opt/scripts -B /scratch/project_2000403/toukopur/data:$PWD/data -B /scratch/project_2000403/toukopur/data:/data torchfort_bisonflex.sif bash
#apptainer exec --nv -B $PENCIL_HOME:$PENCIL_HOME -B scripts:/opt/scripts -B /scratch/project_2000403/$USER/data:$PWD/data -B /scratch/project_2000403/$USER/data:/data torchfort_bisonflex.sif bash /opt/scripts/run.sh

apptainer exec --nv -B $PENCIL_HOME:$PENCIL_HOME -B scripts:/opt/scripts -B /scratch/project_2000403/sgiridha/data_conv-slab-torch/data:$PENCIL_HOME/samples/conv-slab-torch/data -B /scratch/project_2000403/sgiridha/data_conv-slab-torch/data:/data -B /scratch/project_2000403/sgiridha/snapshots_conv-slab-torch/snapshots:/$PENCIL_HOME/samples/conv-slab-torch/snapshots -B /scratch/project_2000403/sgiridha/snapshots_conv-slab-torch/snapshots:/snapshots torchfort_latest_withdf5.sif bash rm -rf data bash rm -rf snapshots


apptainer exec --nv -B $PENCIL_HOME:$PENCIL_HOME -B scripts:/opt/scripts -B /scratch/project_2000403/sgiridha/data_conv-slab-torch/data:$PENCIL_HOME/samples/conv-slab-torch/data -B /scratch/project_2000403/sgiridha/data_conv-slab-torch/data:/data -B /scratch/project_2000403/sgiridha/snapshots_conv-slab-torch/snapshots:/$PENCIL_HOME/samples/conv-slab-torch/snapshots -B /scratch/project_2000403/sgiridha/snapshots_conv-slab-torch/snapshots:/snapshots torchfort_latest_withdf5.sif bash



apptainer exec --nv -B $PENCIL_HOME:$PENCIL_HOME -B /appl:/appl -B scripts:/opt/scripts -B /scratch/project_2000403/sgiridha/data_conv-slab-torch/data:$PENCIL_HOME/samples/conv-slab-torch/data -B /scratch/project_2000403/sgiridha/data_conv-slab-torch/data:/data -B /scratch/project_2000403/sgiridha/snapshots_conv-slab-torch/snapshots:/$PENCIL_HOME/samples/conv-slab-torch/snapshots -B /scratch/project_2000403/sgiridha/snapshots_conv-slab-torch/snapshots:/snapshots torchfort_latest_withdf5.sif bash /opt/scripts/run.sh
