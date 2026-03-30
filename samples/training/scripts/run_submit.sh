#!/bin/bash
#SBATCH --job-name=test_torchfort
#SBATCH --account=project_2016901
#SBATCH --partition=gputest
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --gres=gpu:v100:1
#SBATCH --mem=128G
#SBATCH -o slurm-%x_%J.out
####SBATCH --mail-type=ALL
#SBATCH --output=train.out


# current sample name
sample_name=helical-MHDturb

# training or inference
#mode=training_multi
mode=inference

ml_model=unet


#data_src="/scratch/project_2000403/$USER/data_$sample_name-multi-gpu/data"
data_src="/scratch/project_2000403/$USER/data_$sample_name/data"
snap_src="/scratch/project_2000403/$USER/snapshots_$sample_name/snapshots"
sample_src="$PENCIL_HOME/samples/training/$sample_name/$mode"

run_dir="$sample_src/run.in"
dt_train_min=0.14
dt_train_max=0.175


module load python-data
module load pytorch/2.4

export PYTHONPATH=$PYTHONPATH:$PENCIL_HOME/samples/training/models





if [[ "$mode" == "inference" ]]; then
	
		if [[ ! -f "$data_src/training/stationary.pt" || ! -f  "$data_src/training/stats_current_output.pt" ]]; then
    	echo "Error: Either model file: stationary.pt or stats file: stats_current_output.pt is not found."
    	exit 1
		fi

    srun -n 1 python3 -c "from build_files import ptTObin; ptTObin('$data_src/training')"

else
	srun -n 1 python3 -c "from build_files import build_model; build_model('$data_src/training', '$data_src/training', '$ml_model')"

	srun -n 1 python3 -c "from build_files import build_loss; build_loss('$data_src/training', '$data_src/training')"

	srun -n 1 python3 -c "from build_files import rand_dt_train; rand_dt_train('$run_dir', $dt_train_min, $dt_train_max)"

	srun -n 1 python3 -c "from build_files import build_torchfort_config; build_torchfort_config('$data_src/training', '$ml_model')"
fi

apptainer exec --nv \
		-B $PENCIL_HOME:$PENCIL_HOME \
		-B /appl:/appl \
		-B $sample_src/scripts/run.sh:/opt/scripts/run.sh \
		-B $data_src:$sample_src/data \
		-B $snap_src:$sample_src/snapshots \
		torchfort_latest_withdf5.sif bash /opt/scripts/run.sh


