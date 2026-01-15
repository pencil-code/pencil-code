#!/bin/bash
#SBATCH --job-name=create_ml_model
#SBATCH --account=project_2000403
#SBATCH --partition=gputest
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:v100:1
#SBATCH --mem=16G
#SBATCH --output=model_com.out

module load python-data
module load pytorch/2.4

# current sample name
sample_name=helical-MHDturb 

data_src="/scratch/project_2000403/$USER/data_$sample_name/data/training/"

srun python3 $PENCIL_HOME/samples/training/models/build_training_files.py $data_src
