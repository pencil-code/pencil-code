#!/bin/bash -l
#SBATCH --output=test.out
#SBATCH --partition=dev-g  # Partition (queue) name
#SBATCH --gpus-per-node=1       # Allocate one gpu per MPI rank
##SBATCH --partition=debug  # Partition (queue) name
#SBATCH --nodes=2 # Total number of nodes 
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:15:00       # Run time (d-hh:mm:ss)
#SBATCH --account=project_462000613 # Project for billing
#SBATCH --cpus-per-task=7
##SBATCH --exclusive
#SBATCH --mem=0

export OMP_NUM_THREADS=7
export OMP_PROC_BIND=close,spread
export OMP_MAX_ACTIVE_LEVELS=2
export OMP_WAIT_POLICY=PASSIVE
#./start.csh
export MPICH_GPU_SUPPORT_ENABLED=1
#srun ./src/run.x
#./run.csh
#
#To Ondrej modify this to source odop and the python stuff it requires
source ~/new-monitoring/env/bin/activate

#To Ondrej modify this to point your odop installation
export ODOP_REAL_PATH=~/new-monitoring/odop/odop/odop_obs
export ODOP_PATH=$ODOP_REAL_PATH/
rm  -f $ODOP_REAL_PATH/metric_database/*.csv
rm  -f $ODOP_REAL_PATH/metric_database/*.png
rm  -f $ODOP_REAL_PATH/logs/*.txt
#export LD_PRELOAD=./src/libPC.so:/users/toukopur/pencil-code/pencil-private/projects/PC-A/cpu/src/libPCStart.so
export LD_PRELOAD=./src/libPC.so
#export LD_PRELOAD=/users/toukopur/pencil-code/pencil-private/projects/PC-A/cpu/src/libPCStart.so
#
#To Ondrej: change this to point the the sample you want to run alongside on the CPU with this one
export PC_CPU_SAMPLE=/users/toukopur/pencil-code/pencil-private/projects/PC-A/cpu
#TP: should be same as --nodes in here
export PC_CPU_SAMPLE_NODES=1
#TP: for one-to-one mapping between GPUs and CPUs
export PC_CPU_SAMPLE_RANKS_PER_NODE=8

#TP: this does not matter if one is not using Allas.
export OS_AUTH_URL=https://pouta.csc.fi:5001/v3
# With the addition of Keystone we have standardized on the term **project**
# as the entity that owns the resources.
export OS_PROJECT_ID=9665d280d1e34e8699056c8dd24e7938
export OS_PROJECT_NAME="project_2000403"
export OS_USER_DOMAIN_NAME="Default"
if [ -z "$OS_USER_DOMAIN_NAME" ]; then unset OS_USER_DOMAIN_NAME; fi
export OS_PROJECT_DOMAIN_ID="default"
if [ -z "$OS_PROJECT_DOMAIN_ID" ]; then unset OS_PROJECT_DOMAIN_ID; fi
# unset v2.0 items in case set
unset OS_TENANT_ID
unset OS_TENANT_NAME
# In addition to the owning entity (tenant), OpenStack stores the entity
# performing the action as the **user**.
export OS_USERNAME="USERNAME"
# With Keystone you pass the keystone password.
export OS_PASSWORD="PASSWORD"
export OS_IDENTITY_API_VERSION="3"


srun python call.py
#srun python op-tasks/reduce.py
#srun python /users/toukopur/pencil-code/samples/gputest/op-tasks/reduce.py


