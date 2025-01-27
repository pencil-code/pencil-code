#!/bin/bash
export PENCIL_HOME=/users/mreinhar/pencil-code
cd $PENCIL_HOME
source sourceme.sh
cd $PENCIL_HOME/samples/conv-slab-torch
export PATH=${PATH}:$PENCIL_HOME/bin:$PENCIL_HOME/utils:$PENCIL_HOME/utils/axel:$PENCIL_HOME/utils/xiangyu:$PENCIL_HOME/remesh/bin:$PENCIL_HOME/src/scripts
echo PATH=$PATH
##export OMPI_ALLOW_RUN_AS_ROOT=1
##export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
export TORCHFORT_LOGDIR=data
#pc_run start
pc_run run
#src/read_all_videofiles.x
