#!/bin/bash
#SBATCH --job-name=torch
#SBATCH --account=project_2016901
#SBATCH --partition=gpumedium
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=72
#SBATCH --gres=gpu:gh200:1
####SBATCH --mem=12G
#SBATCH -o slurm-%x_%J.out
####SBATCH --mail-type=ALL
#SBATCH --output=train.out

echo "***********************************************"
echo "The Job ID: $SLURM_JOB_ID"
echo "***********************************************"


unset CMAKE_PREFIX_PATH
unset CMAKE_INCLUDE_PATH
unset CMAKE_LIBRARY_PATH

#module purge
module load nvhpc
module load cmake


export CC=mpicc
export CXX=mpicxx
export FC=mpif90a

PROJ_DIR="/projappl/project_2016901/$USER"


HDF5_PATH="$PROJ_DIR/hdf5_nvhpc_parallel"
TORCHFORT_PATH="$PROJ_DIR/torchfort-install"
LIBTORCH_PATH="$PROJ_DIR/libtorch"
YAMLCPP_PATH="$PROJ_DIR/yaml-cpp-install"

export PYTHONUSERBASE="$PROJ_DIR/my-python-env"

CUDA_LIB="$NVHPC/cuda/lib64"
NCCL_LIB=$(ls -d $NVHPC/Linux_aarch64/*/comm_libs/*/nccl/lib 2>/dev/null | tail -n 1)
export LD_LIBRARY_PATH="$:$LD_LIBRARY_PATH:$HDF5_PATH/lib:$TORCHFORT_PATH/lib:$LIBTORCH_PATH/lib:$YAMLCPP_PATH/lib64:$CUDA_LIB:$NCCL_LIB"

export NVCC_PREPEND_FLAGS="--expt-relaxed-constexpr ${NVCC_PREPEND_FLAGS}"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_PROC_BIND=spread
export OMP_MAX_ACTIVE_LEVELS=2
export OMP_WAIT_POLICY=PASSIVE

export OMPI_MCA_coll_hcoll_enable=0


rm data/training/stationary.pt
rm wandb_helper.log 


export WANDB_API_KEY="xxx"
export WANDB_LOGGING_DIR="data/training/wandb_log/wandb_${SLURM_JOB_ID}"
mkdir -p "${WANDB_LOGGING_DIR}"

export TORCHFORT_LOGDIR="${WANDB_LOGGING_DIR}"

mkdir -p "${TORCHFORT_LOGDIR}"

#mkdir -p "${TORCHFORT_LOGDIR}" "${WANDB_LOGGING_DIR}"

TORCHFORT_INSTALL_DIR="/projappl/project_2016901/sgiridha/torchfort-install"

python -u ${TORCHFORT_INSTALL_DIR}/bin/python/wandb_helper.py \
      --wandb_dir=${WANDB_LOGGING_DIR} \
      --wandb_group="roihutest" \
      --wandb_project="torchfort" \
      --wandb_entity="girishreyas2000" \
      --run_tag="job-${SLURM_JOB_ID}" \
      --timeout=100 > wandb_helper.log 2>&1 &
WANDB_PID=$!

pc_run start
pc_run run

#sleep 10
#kill -TERM $WANDB_PID 2>/dev/null || true
#wait $WANDB_PID 2>/dev/null || true

wait $WANDB_PID
