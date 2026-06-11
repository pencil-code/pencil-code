#!/bin/bash
#SBATCH --job-name=build_torchfort
#SBATCH --account=project_2016901
#SBATCH --partition=gpu
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --gres=gpu:v100:1
#SBATCH --mem=128G
#SBATCH -o out.out

set -e

BASE_DIR=$(pwd)
INSTALL_PREFIX=$BASE_DIR/torchfort-install
YAML_CPP_PREFIX=$BASE_DIR/yaml-cpp-install
YAML_CPP_SRC=$BASE_DIR/yaml-cpp
TORCHFORT_SRC=$BASE_DIR/TorchFort
BUILD_DIR=$BASE_DIR/torchfort-build
LIBTORCH_DIR=$BASE_DIR/libtorch         
CUDNN_DIR=$BASE_DIR/cudnn-linux-x86_64-8.7.0.84_cuda11-archive
CUSPARSELT_DIR=$BASE_DIR/libcusparse_lt-linux-x86_64-0.4.0.7-archive
HDF5_SRC=$BASE_DIR/hdf5-1.12.2
HDF5_INSTALL=$BASE_DIR/hdf5_nvhpc_parallel
JOBS=7

NVHPC_BASE=/appl/opt/nvhpc/Linux_x86_64/24.11
NVHPC_BIN=$NVHPC_BASE/compilers/bin
NVHPC_CUDA=$NVHPC_BASE/cuda/11.8
NVHPC_NCCL=$NVHPC_BASE/comm_libs/11.8/nccl
NVHPC_MPI=$NVHPC_BASE/comm_libs/11.8/openmpi4/openmpi-4.1.5

# ── Modules ───────────────────────────────────────────────────────────────────
module purge
module load gcc/11.3.0
module load git/2.35.2
module load cmake/3.23.1
module use /appl/opt/nvhpc/modulefiles
module load nvhpc-hpcx-cuda11/24.11

mkdir -p $BASE_DIR/build_tmp
export TMPDIR=$BASE_DIR/build_tmp
export TEMP=$BASE_DIR/build_tmp
export TMP=$BASE_DIR/build_tmp

export CUDA_HOME=$NVHPC_CUDA
export CUDA_PATH=$NVHPC_CUDA
export PATH=$NVHPC_BIN:$NVHPC_CUDA/bin:$PATH
export LD_LIBRARY_PATH=$CUDNN_DIR/lib:$CUSPARSELT_DIR/lib:$NVHPC_NCCL/lib:$NVHPC_MPI/lib:$NVHPC_CUDA/lib64:$LD_LIBRARY_PATH

# ── 1/7 Download dependencies if missing ──────────────────────────────────────
if [ ! -d "$CUDNN_DIR" ]; then
    echo "[1/7] Downloading cuDNN 8.7.0 for CUDA 11.8..."
    wget https://developer.download.nvidia.com/compute/cudnn/redist/cudnn/linux-x86_64/cudnn-linux-x86_64-8.7.0.84_cuda11-archive.tar.xz
    tar xf cudnn-linux-x86_64-8.7.0.84_cuda11-archive.tar.xz
    rm cudnn-linux-x86_64-8.7.0.84_cuda11-archive.tar.xz
fi

if [ ! -d "$CUSPARSELT_DIR" ]; then
    echo "[1/7] Downloading cuSPARSELt 0.4.0.7..."
    wget https://developer.download.nvidia.com/compute/cusparselt/redist/libcusparse_lt/linux-x86_64/libcusparse_lt-linux-x86_64-0.4.0.7-archive.tar.xz
    tar xf libcusparse_lt-linux-x86_64-0.4.0.7-archive.tar.xz
    rm libcusparse_lt-linux-x86_64-0.4.0.7-archive.tar.xz
fi

if [ ! -d "$LIBTORCH_DIR" ]; then
    echo "[1/7] Downloading LibTorch 2.2.2 cu118 pre-cxx11 ABI..."
    wget https://download.pytorch.org/libtorch/cu118/libtorch-shared-with-deps-2.2.2%2Bcu118.zip -O libtorch.zip
    unzip -q libtorch.zip
    rm libtorch.zip
fi

# ── 2/7 Build HDF5 (parallel, with NVHPC MPI wrappers) ───────────────────────
if [ ! -f "$HDF5_INSTALL/lib/libhdf5.so" ]; then
    echo "[2/7] Building HDF5..."
    if [ ! -d "$HDF5_SRC" ]; then
        wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.2/src/hdf5-1.12.2.tar.gz
        tar xzf hdf5-1.12.2.tar.gz
    fi
    cd $HDF5_SRC
    export CC=mpicc
    export FC=mpifort
    export CXX=mpic++
    export CFLAGS="-fPIC"
    export FCFLAGS="-fPIC"
    export CXXFLAGS="-fPIC -D_GLIBCXX_USE_CXX11_ABI=0"
    ./configure \
        --prefix=$HDF5_INSTALL \
        --enable-fortran \
        --enable-parallel \
        --enable-shared \
        --disable-tests
    make -j$JOBS
    make install
    unset CC FC CXX CFLAGS FCFLAGS CXXFLAGS
    cd $BASE_DIR
fi

# ── 3/7 Build yaml-cpp (pre-cxx11 ABI) ───────────────────────────────────────
echo "[3/7] Building yaml-cpp..."
if [ ! -d "$YAML_CPP_SRC" ]; then
    git clone https://github.com/jbeder/yaml-cpp.git $YAML_CPP_SRC
fi
rm -rf $YAML_CPP_SRC/build
mkdir -p $YAML_CPP_SRC/build && cd $YAML_CPP_SRC/build
cmake .. \
    -DCMAKE_INSTALL_PREFIX=$YAML_CPP_PREFIX \
    -DYAML_BUILD_SHARED_LIBS=OFF \
    -DCMAKE_CXX_COMPILER=g++ \
    -DCMAKE_C_COMPILER=gcc \
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
    -DCMAKE_CXX_FLAGS="-D_GLIBCXX_USE_CXX11_ABI=0"
make -j$JOBS install
cd $BASE_DIR

# ── 4/7 Clone TorchFort ───────────────────────────────────────────────────────
if [ ! -d "$TORCHFORT_SRC" ]; then
    git clone https://github.com/NVIDIA/TorchFort.git $TORCHFORT_SRC
fi

# ── 5/7 Configure TorchFort ───────────────────────────────────────────────────
echo "[5/7] Configuring TorchFort..."
rm -rf $BUILD_DIR
mkdir -p $BUILD_DIR && cd $BUILD_DIR

TORCH_CMAKE_PATH=$LIBTORCH_DIR/share/cmake/Torch

cmake $TORCHFORT_SRC \
    -DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX \
    -DTORCHFORT_YAML_CPP_ROOT=$YAML_CPP_PREFIX \
    -DYAML_CPP_LIBRARY=$YAML_CPP_PREFIX/lib64/libyaml-cpp.a \
    -DTORCHFORT_NCCL_ROOT=$NVHPC_NCCL \
    -DTORCHFORT_BUILD_EXAMPLES=OFF \
    -DCMAKE_C_COMPILER=gcc \
    -DCMAKE_CXX_COMPILER=g++ \
    -DCMAKE_Fortran_COMPILER=$NVHPC_BIN/nvfortran \
    -DCMAKE_CUDA_COMPILER=$NVHPC_BIN/nvcc \
    -DTorch_DIR=$TORCH_CMAKE_PATH \
    -DCMAKE_PREFIX_PATH="$LIBTORCH_DIR;$YAML_CPP_PREFIX;$CUDNN_DIR" \
    -DCMAKE_CUDA_FLAGS="-allow-unsupported-compiler" \
    -DCMAKE_CXX_FLAGS="-D_GLIBCXX_USE_CXX11_ABI=0" \
    -DCMAKE_C_FLAGS="-D_GLIBCXX_USE_CXX11_ABI=0" \
    -DCMAKE_EXE_LINKER_FLAGS="-L$NVHPC_CUDA/lib64/stubs -lcuda" \
    -DCMAKE_SHARED_LINKER_FLAGS="-L$NVHPC_CUDA/lib64/stubs -lcuda" \
    -DCUDNN_INCLUDE_DIR="$CUDNN_DIR/include" \
    -DCUDNN_LIBRARY="$CUDNN_DIR/lib/libcudnn.so" \
    -DCUDAToolkit_ROOT=$NVHPC_CUDA \
    -DCUDA_TOOLKIT_ROOT_DIR=$NVHPC_CUDA \
    -DMPI_HOME=$NVHPC_MPI \
    -DMPI_C_INCLUDE_PATH=$NVHPC_MPI/include \
    -DMPI_C_LIBRARIES=$NVHPC_MPI/lib/libmpi.so \
    -DMPI_CXX_INCLUDE_PATH=$NVHPC_MPI/include \
    -DMPI_CXX_LIBRARIES=$NVHPC_MPI/lib/libmpi.so \
    -DMPI_Fortran_INCLUDE_PATH=$NVHPC_MPI/include \
    -DMPI_Fortran_LIBRARIES="$NVHPC_MPI/lib/libmpi_usempif08.so;\
$NVHPC_MPI/lib/libmpi_usempi_ignore_tkr.so;\
$NVHPC_MPI/lib/libmpi_mpifh.so;\
$NVHPC_MPI/lib/libmpi.so"

# ── 6/7 Build ─────────────────────────────────────────────────────────────────
echo "[6/7] Building TorchFort..."
cd $BUILD_DIR

if [ ! -f "$YAML_CPP_PREFIX/lib/libyaml-cpp.a" ] && \
   [ ! -f "$YAML_CPP_PREFIX/lib64/libyaml-cpp.a" ]; then
    cd $YAML_CPP_SRC/build && make install && cd $BUILD_DIR
fi

make -j1

# ── 7/7 Install ───────────────────────────────────────────────────────────────
echo "[7/7] Installing..."
make install

echo "============================================="
echo "SUCCESS!"
echo "TorchFort : $INSTALL_PREFIX"
echo "HDF5      : $HDF5_INSTALL"
echo "cuDNN     : $CUDNN_DIR"
echo "cuSPARSELt: $CUSPARSELT_DIR"
echo "LibTorch  : $LIBTORCH_DIR"
echo "============================================="
