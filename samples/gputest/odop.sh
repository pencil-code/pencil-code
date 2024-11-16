#!/bin/bash -l
#SBATCH --output=test.out
#SBATCH --partition=dev-g  # Partition (queue) name
#SBATCH --nodes=1 # Total number of nodes 
#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-node=1       # Allocate one gpu per MPI rank
#SBATCH --time=00:15:00       # Run time (d-hh:mm:ss)
#SBATCH --account=project_462000613 # Project for billing
#SBATCH --cpus-per-task=7
#SBATCH --exclusive
#SBATCH --mem=0

export OMP_NUM_THREADS=7
export OMP_PROC_BIND=close,spread
export OMP_MAX_ACTIVE_LEVELS=2
export OMP_WAIT_POLICY=PASSIVE
./start.csh
export MPICH_GPU_SUPPORT_ENABLED=1
#srun ./src/run.x
#./run.csh
source ~/new-monitoring/venv/bin/activate
export ODOP_REAL_PATH=~/new-monitoring/odop/odop/odop_obs
export ODOP_PATH=$ODOP_REAL_PATH/
rm  -f $ODOP_REAL_PATH/metric_database/*.csv
rm  -f $ODOP_REAL_PATH/metric_database/*.png
rm  -f $ODOP_REAL_PATH/logs/*.txt
export LD_PRELOAD=./src/libPC.so

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
export OS_USERNAME="USER_NAME_HERE"
# With Keystone you pass the keystone password.
export OS_PASSWORD="PASS_WORD_HERE"
export OS_IDENTITY_API_VERSION="3"


srun python call.py


