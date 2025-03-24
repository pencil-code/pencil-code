#!/bin/bash
export PATH=/opt/hdf5/bin:$PATH 
export PENCIL_HOME=/users/$USER/pencil-code
CWDIR=$PWD
cd $PENCIL_HOME
unset PENCIL_HOME
source sourceme.sh
export PATH=$PATH:$PENCIL_HOME/bin:$PENCIL_HOME/utils:$PENCIL_HOME/utils/axel:$PENCIL_HOME/utils/xiangyu:$PENCIL_HOME/remesh/bin:$PENCIL_HOME/src/scripts
#echo command path = $PATH
echo library path = $LD_LIBRARY_PATH
cd $CWDIR
export TMPDIR=/opt/tmpdir
echo $TMPDIR
#bash
#pc_build clean -f $PENCIL_HOME/config/compilers/separate/nvidia-fortran.conf TORCHFORT_PATH=/opt/torchfort HDF5_PATH=/opt/hdf5 FFLAGS+="-O0 -g -traceback"
pc_build -f $PENCIL_HOME/config/compilers/separate/nvidia-fortran.conf TORCHFORT_PATH=/opt/torchfort HDF5_PATH=/opt/hdf5 FFLAGS+="-O0 -g -traceback"
#cd src/astaroth && make
#make


#pc_build -f $PENCIL_HOME/config/compilers/separate/nvidia-fortran.conf -t read_all_videofiles

#FFLAGS+="-I/opt/hdf5/include -g -traceback" LDFLAGS="-L/opt/hdf5/lib/ -lhdf5_fortran" -s
#/users/mreinhar/pencil-code/bin/pc_build -f /users/mreinhar/pencil-code/config/compilers/GNU-GCC_MPI.conf FFLAGS+="-cuda -noswitcherror -I/opt/hdf5/include -I/opt/torchfort/include/ -g -traceback" LDFLAGS="-L/opt/torchfort/lib/ -L/opt/hdf5/lib/ -ltorchfort_fort -lhdf5_fortran -ltorchfort -cuda" -t read_all_videofiles -s
