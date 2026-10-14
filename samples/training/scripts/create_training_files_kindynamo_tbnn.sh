#!/bin/bash
#SBATCH --job-name=create_ml_model
#SBATCH --account=project_2016901
#SBATCH --partition=gputest
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:gh200:1
#SBATCH --mem=16G
#SBATCH --output=model_com.out

module load python-data
module load python-pytorch/2.10

# Project name
project_name=project_2016901

# Sample name
sample_name=kin-dynamo

# Physics setup name
physics_setup_name=default-physics-setup

# Model name
model_name=kindynamo_tbnn

# Experiment name
experiment_name=DEBUGGING

data_src="/scratch/$project_name/$USER/$sample_name/$physics_setup_name/$physics_setup_src/$model_name/$experiment_name/data"
training_src="$data_src/training/"

mkdir -p $data_src
mkdir -p $training_src

echo $data_src
echo $training_src

if [ ! -L data ]; then
   ln -s $data_src data
fi

srun python3 $PENCIL_HOME/samples/training/models/kindynamo_build_training_files.py $training_src $model_name
cp $PENCIL_HOME/samples/training/models/ptTobin.py $training_src
