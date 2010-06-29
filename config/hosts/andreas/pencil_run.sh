#!/bin/sh
#This is the the executable file called by pencil.submit
#Please note that the enviroment variable might need to be 
#changed at the node for pc_run to be called.

export PENCIL_HOME=/hpc/scratch/astro/users/aos2112/pencil-code

_sourceme_quiet=1; . $PENCIL_HOME/sourceme.sh; unset _sourceme_quiet

echo $PENCIL_HOME


#start.csh

#run.csh

pc_run start
pc_run run
