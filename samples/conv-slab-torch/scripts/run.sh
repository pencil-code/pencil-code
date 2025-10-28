#!/bin/bash
#export PENCIL_HOME=/users/mreinha
CWD=$PWD
cd $PENCIL_HOME
source sourceme.sh
cd $CWD
export PATH=${PATH}:$PENCIL_HOME/bin:$PENCIL_HOME/utils:$PENCIL_HOME/utils/axel:$PENCIL_HOME/utils/xiangyu:$PENCIL_HOME/remesh/bin:$PENCIL_HOME/src/scripts
export PATH=${PATH}:./src/astaroth:./src/astaroth/submodule/build/src/core:./src/astaroth/submodule/build/src/core/kernels
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:./src/astaroth:./src/astaroth/submodule/build/src/core:./src/astaroth/submodule/build/src/core/kernels
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/appl/spack/v018/install-tree/gcc-11.3.0/cuda-11.7.0-zucvj4/lib64/libcufft.so.10
export PATH=${PATH}:/appl/spack/v018/install-tree/gcc-11.3.0/cuda-11.7.0-zucvj4/lib64/libcufft.so.10
export FC=mpifort
export CC=mpicc
export CXX=mpic++
echo PATH=$PATH
echo LD_LIB=$LD_LIBRARY_PATH
##export OMPI_ALLOW_RUN_AS_ROOT=1
##export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
export TORCHFORT_LOGDIR=data

export OMP_NUM_THREADS=7
export OMP_PROC_BIND=close,spread
export OMP_MAX_ACTIVE_LEVELS=2
export OMP_WAIT_POLICY=PASSIVE

#ldd  src/run.x
#bash
ls   data
ls   data/training

#pc_run start
#pc_run run
#mpiexec -n 1 --bind-to none --cpus-per-proc 7 ./src/start.x
mpiexec -n 1 --bind-to none --cpus-per-proc 7 ./src/run.x
#src/read_all_videofiles.x
