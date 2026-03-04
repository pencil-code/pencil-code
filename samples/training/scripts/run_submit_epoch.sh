#!/bin/bash
#SBATCH --job-name=test_torchfort
#SBATCH --account=project_2016901
#SBATCH --partition=gpu
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --gres=gpu:v100:1
#SBATCH --mem=32G
#SBATCH -o slurm-%x_%J.out
####SBATCH --mail-type=ALL
#SBATCH --output=train.out


# current sample name
sample_name=conv-slab

# training or inference
mode=training

# unet or fno
ml_model=unet

dt_train_min=0.000233
dt_train_max=0.000311


data_src="/scratch/project_2000403/$USER/data_$sample_name/data"
snap_src="/scratch/project_2000403/$USER/snapshots_$sample_name/snapshots"
sample_src="$PENCIL_HOME/samples/training/$sample_name/$mode"

run_dir="$sample_src/run.in"


module load python-data
module load pytorch/2.4

export PYTHONPATH=$PYTHONPATH:$PENCIL_HOME/samples/training/models

EPOCHS=5

for ((i=1; i<=EPOCHS; i++)); do
    echo "=== Starting Epoch $i / $EPOCHS ==="

		if [ $i -le 2 ]; then
			# fresh start with known mean and std deviation instead of cold start (which might tend to overfit to the largest values)
			srun python3 -c "from build_files import build_model; build_model('$data_src/training', '$data_src/training', '$ml_model')"
		fi

		
		# resume the running norm values instead of starting from scratch
		srun python3 -c "from build_files import build_loss; build_loss('$data_src/training', '$data_src/training')"


		# add some randomization for better generalization
		srun python3 -c "from build_files import rand_dt_train; rand_dt_train('$run_dir', $dt_train_min, $dt_train_max)"

		srun python3 -c "from build_files import build_torchfort_config; build_torchfort_config('$data_src/training', '$ml_model')"

    
		apptainer exec --nv \
				-B $PENCIL_HOME:$PENCIL_HOME \
				-B /appl:/appl \
				-B $sample_src/scripts/run.sh:/opt/scripts/run.sh \
				-B $data_src:$sample_src/data \
				-B $sample_src/scripts:$sample_src/scripts \
				-B $snap_src:$sample_src/snapshots \
				torchfort_latest_withdf5.sif bash /opt/scripts/run.sh



    echo "=== Finished Epoch $i ==="
		
		if [ $i -ne 1 ]; then
    	# Copy model checkpoints and logs for this epoch
    	cp "${data_src}/training/stationary.pt" \
      	 "${data_src}/training/${ml_model}/${ml_model}_epoch_256_${i}.pt"

    	cp train_loss_0.csv \
      	 "${data_src}/training/${ml_model}/train_loss_rank_0_epoch_${i}.csv"

    	cp val_loss_0.csv \
      	 "${data_src}/training/$ml_model/val_loss_rank_0_epoch_${i}.csv"
		fi
		
		# discard the first model because it has already tried to optimize for largest values
		if [ $i -eq 1 ]; then
			rm -f "${data_src}/training/stationary.pt"
		fi
		
		cp stats_current_output.pt \
			"${data_src}/training/stats_current_output.pt"

		cp stats_current_input.pt \
			"${data_src}/training/stats_current_input.pt"
			
done

cp train.out \
	"${data_src}/training/${ml_model}/train_epoch_${EPOCHS}.out"



cp stats_current_output.pt \
"${data_src}/training/${ml_model}/stats_current_output.pt"

cp stats_current_input.pt \
	"${data_src}/training/${ml_model}/stats_current_input.pt"

