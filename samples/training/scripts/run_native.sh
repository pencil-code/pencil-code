#!/bin/bash
#SBATCH --job-name=torch
#SBATCH --account=project_2016901
#SBATCH --partition=gputest
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=7
#SBATCH --gres=gpu:v100:1
#SBATCH --mem=128G
#SBATCH -o slurm-%x_%J.out
####SBATCH --mail-type=ALL
#SBATCH --output=train.out

echo "***********************************************"
echo "The Job ID: $SLURM_JOB_ID"
echo "***********************************************"


unset CMAKE_PREFIX_PATH
unset CMAKE_INCLUDE_PATH
unset CMAKE_LIBRARY_PATH
unset CC
unset CXX
unset FC

module purge
module use /appl/opt/nvhpc/modulefiles
module load nvhpc-hpcx-cuda11/24.11
module load gcc/11.3.0
module load cmake

export NVCC_PREPEND_FLAGS="--compiler-bindir=$(which g++) -Xcompiler -gz=none"
GCC11_LIBDIR=$(gcc --print-file-name=libgcc_s.so.1 | xargs dirname)
export LD_LIBRARY_PATH=$GCC11_LIBDIR:$LD_LIBRARY_PATH

export NVHPC_CUDA_HOME=/appl/opt/nvhpc/Linux_x86_64/24.11/cuda/11.8
GCC11_BIN=$(dirname $(which gcc))
export PATH=$GCC11_BIN:$PATH

export CC=nvc
export CXX=nvc++
export FC=nvfortran

export CXXFLAGS="-std=c++17"
export CMAKE_CXX_FLAGS="-std=c++17"

export LD_LIBRARY_PATH=/scratch/project_2016901/sgiridha/hdf5_nvhpc_parallel/lib:$LD_LIBRARY_PATH

export LD_LIBRARY_PATH=/appl/opt/nvhpc/Linux_x86_64/24.11/cuda/11.8/lib64:$LD_LIBRARY_PATH

export CUFILE_LIB_DIR=/appl/opt/nvhpc/Linux_x86_64/24.11/cuda/11.8/targets/x86_64-linux/lib
export HDF5_PATH=/scratch/project_2016901/sgiridha/hdf5_nvhpc_parallel/lib
export CUDA_PATH=/appl/opt/nvhpc/Linux_x86_64/24.11/cuda/11.8/lib64
export LD_LIBRARY_PATH=$HDF5_PATH:$CUDA_PATH:$LD_LIBRARY_PATH:$CUFILE_LIB_DIR

export OMP_NUM_THREADS=7
export OMP_PROC_BIND=close,spread
export OMP_MAX_ACTIVE_LEVELS=2
export OMP_WAIT_POLICY=PASSIVE

export OMPI_MCA_coll_hcoll_enable=0

export NCCL_DEBUG=INFO
export TORCH_DISTRIBUTED_DEBUG=INFO


rm data/training/stationary.pt

export CUDA_HOME=$NVHPC_BASE/cuda/11.8
export CUDA_PATH=$NVHPC_BASE/cuda/11.8
export CUDAToolkit_ROOT=$NVHPC_BASE/cuda/11.8
export CUDA_TOOLKIT_ROOT_DIR=$NVHPC_BASE/cuda/11.8
export CPATH=$NVHPC_BASE/cuda/11.8/include:$CPATH
export PATH=$NVHPC_BASE/cuda/11.8/bin:$PATH

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$NVHPC_BASE/cuda/11.8/lib64:\
$NVHPC_BASE/cuda/11.8/targets/x86_64-linux/lib:\

#pc_run start
pc_run run
