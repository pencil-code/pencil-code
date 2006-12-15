#!/bin/csh

# Name:   run_gdb.csh
# Author: wd (Wolfgang.Dobler@ucalgary.ca)
# Date:   29-Jun-2005
# Description:
#   Start xterm and run debugger there, starting the program src/run.x and
#   after abortion looking at the stack trace. Useful for parallel
#   debugging, where we want output from each process in a separate window
#   and want to avoid the need to type something in each window.
# Usage:
#   mpirun N -np 4  -x DISPLAY ${PENCIL_HOME}/samples/run_gdb.csh
# Keywords:
#   parallel debugging
# Limitations:
#   Only works with LAM MPI.

xterm -T "$LAMRANK" -geometry 80x60 -e '( printf "%s\n" "r" "backtrace"; cat ) | gdb src/run.x'


# End of file run_gdb.csh
