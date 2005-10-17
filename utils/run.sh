#!/bin/sh

#                       run.csh
#                      ---------
#   Run ./run.csh, but catch USR1 and USR2 signals (csh apparently cannot
#   do that) and make the code STOP in time. This assumes that SIGUSR1/2 is
#   sent ahead of SIGKILL as is the case with Gridengine's -notify option
#   (if the queue is configured accordingly) or with Gridengine's s_rt
#   resource.

# Join stderr and stout:
#$ -j y -o run.log -notify
#$ -notify
#$ -m beas
#@$-eo
#@$-o run.log
#
# Work in submit directory (SGE):
#$ -cwd

# Spawn Perl script that writes STOP file when it receives USR1 or USR2 signal
./perl_monitor &

# Protect run.csh (and thus in particular mpirun and run.x subprocesses)
# from signals. The code will stop when it sees the STOP file.
trap '' usr1
trap '' usr2


echo "run.sh: starting run.csh"
./run.csh
echo "run.sh: Done"


