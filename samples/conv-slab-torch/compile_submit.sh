#!/bin/bash
#SBATCH --job-name=test_torchfort
#SBATCH --account=project_2000403
#SBATCH --partition=test
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
##SBATCH -o slurm-%x_%J.out
####SBATCH --mail-type=ALL
#SBATCH --output=comp.out

#TP: have to bind /appl to get the right CUDA compiler for PC-A
#    the one in the container is not the correct one

apptainer exec --nv -B $PENCIL_HOME:$PENCIL_HOME -B /appl:/appl -B scripts:/opt/scripts -B /users/$USER/tmpdir:/opt/tmpdir torchfort_latest_withdf5.sif bash /opt/scripts/compile.sh

#apptainer exec --nv -B pencil-code:/opt/pencil-code -B scripts:/opt/scripts torchfort_0.2.0.sif bash
#source $PENCIL_HOME/sourceme.sh
#pc_build
