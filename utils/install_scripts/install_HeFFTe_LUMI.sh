#!/bin/bash 
export SCRATCH=/scratch/project_462001500/toukopur
export VER=2.4.1
export MPI_INC=/opt/cray/pe/mpich/9.0.1/ofi/gnu/12.3/include   # adjust if needed (§2)
export CMAKE_PREFIX_PATH="${EBROOTROCM}:${CMAKE_PREFIX_PATH}"
export BUILD_DIR=$SCRATCH/src/heffte-${VER}/build

mkdir -p $SCRATCH/src
cd $SCRATCH/src
test -f v${VER}.tar.gz || wget -q https://github.com/icl-utk-edu/heffte/archive/refs/tags/v${VER}.tar.gz
tar xf v${VER}.tar.gz

cmake -S $SCRATCH/src/heffte-${VER} -B $BUILD_DIR \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=cc \
  -DCMAKE_CXX_COMPILER=CC \
  -DHeffte_ENABLE_FFTW=OFF\
  -DHeffte_ENABLE_ROCM=ON \
  -DHeffte_ENABLE_CUDA=OFF \
  -DHeffte_ENABLE_GPU_AWARE_MPI=ON \
  -DCMAKE_HIP_ARCHITECTURES=gfx90a \
  -DCMAKE_HIP_FLAGS="-I${MPI_INC}"

cd $BUILD_DIR && make -j
#cmake --install $SCRATCH/build/heffte-${VER}-rocm
