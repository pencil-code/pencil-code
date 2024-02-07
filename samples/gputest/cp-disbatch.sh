#!/bin/bash -l
#SBATCH --output=test.out # Name of stdout output file
#SBATCH --partition=debug  # Partition (queue) name
#SBATCH --nodes=1               # Total number of nodes 
#SBATCH --ntasks-per-node=1     # 8 MPI ranks per node, 16 total (2x8)
#SBATCH --time=00:05:00       # Run time (d-hh:mm:ss)
#SBATCH --account=project_462000448  # Project for billing
#SBATCH --cpus-per-task=4



module load CrayEnv
module load rocm
#export MPICH_GPU_SUPPORT_ENABLED=1
#export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK
export SRUN_CPUS_PER_TASK=4
export OMP_MAX_ACTIVE_LEVELS=2
export OMP_NUM_THREADS=4

#srun start.csh
srun $PENCIL_HOME/samples/gputest/src/start.x
rm $PENCIL_HOME/samples/gputest/data/timeseries.dat
srun $PENCIL_HOME/samples/gputest/src/run.x

