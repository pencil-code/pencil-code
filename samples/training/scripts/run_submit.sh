#!/bin/bash
#SBATCH --job-name=test_torchfort
#SBATCH --account=project_2016901
#SBATCH --partition=gpu
#SBATCH --time=05:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --gres=gpu:v100:1
#SBATCH --mem=128G
#SBATCH -o slurm-%x_%J.out
####SBATCH --mail-type=ALL
#SBATCH --output=train.out


# current sample name
sample_name=conv-slab

# training or inference
mode=training

data_src="/scratch/project_2000403/$USER/data_$sample_name/data"
snap_src="/scratch/project_2000403/$USER/snapshots_$sample_name/snapshots"
sample_src="$PENCIL_HOME/samples/training/$sample_name/$mode"


apptainer exec --nv \
		-B $PENCIL_HOME:$PENCIL_HOME \
		-B /appl:/appl \
		-B $sample_src/scripts/run.sh:/opt/scripts/run.sh \
		-B $data_src:$sample_src/data \
		-B $snap_src:$sample_src/snapshots \
		torchfort_latest_withdf5.sif bash /opt/scripts/run.sh

