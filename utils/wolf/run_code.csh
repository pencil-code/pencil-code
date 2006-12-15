#!/bin/csh

# Name:   run_code.csh
# Author: wd (Wolfgang.Dobler@ucalgary.ca)
# Date:   29-Jun-2005
# Description:
#   Start one xterm per process and run code there.
# Usage:
#   mpirun N -np 4  -x DISPLAY ${PENCIL_HOME}/samples/run_code.csh
# Keywords:
#   parallel debugging
# Limitations:
#   Only works with LAM MPI.

xterm -T "$LAMRANK" -geometry 80x40 -e 'src/run.x; sleep 600'


# End of file run_code.csh
