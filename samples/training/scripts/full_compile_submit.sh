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
#sample_name=helical-MHDturb

sample_name=cartesian-convection-kramers

# training or inference
mode=training

data_src="/scratch/project_2000403/$USER/data_$sample_name/data"
snap_src="/scratch/project_2000403/$USER/snapshots_$sample_name/snapshots"
sample_src="$PENCIL_HOME/samples/training/$sample_name/$mode"

#TP: have to bind /appl to get the right CUDA compiler for PC-A
#    the one in the container is not the correct one


apptainer exec --nv \
		-B $PENCIL_HOME:$PENCIL_HOME \
		-B /appl:/appl \
		-B $sample_src/scripts/compile_p1.sh:/opt/scripts/compile_p1.sh \
		-B /users/$USER/tmpdir:/opt/tmpdir \
		$sample_src/torchfort_latest_withdf5.sif bash /opt/scripts/compile_p1.sh


WATCH_FILE="src/astaroth/submodule/build/runtime_build/overrides.h"

rm -r "$WATCH_FILE"

apptainer exec --nv \
		-B $PENCIL_HOME:$PENCIL_HOME \
		-B /appl:/appl \
		-B $sample_src/scripts/run.sh:/opt/scripts/run.sh \
		-B $data_src:$sample_src/data \
		-B $snap_src:$sample_src/snapshots \
		torchfort_latest_withdf5.sif bash /opt/scripts/run.sh &

APP_PID=$!

echo "Checking for creation of $WATCH_FILE..."
while [ ! -f "$WATCH_FILE" ]; do
    sleep 2  
done

if [ -f "$WATCH_FILE" ]; then
    echo "overrides file created! Stopping Apptainer..."
    kill $APP_PID
fi

apptainer exec --nv \
		-B $PENCIL_HOME:$PENCIL_HOME \
		-B /appl:/appl \
		-B $sample_src/scripts/compile_p2.sh:/opt/scripts/compile_p2.sh \
		-B /users/$USER/tmpdir:/opt/tmpdir \
		$sample_src/torchfort_latest_withdf5.sif bash /opt/scripts/compile_p2.sh

