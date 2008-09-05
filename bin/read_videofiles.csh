#!/bin/csh
# CVS: $Id$
#
#  Run read_videofiles.x -- to read_videofiles
#
#PBS -S /bin/csh

# examples of the calling sequence:
# qsub -l nodes=1,walltime=1:00:00 -q workq ../../../bin/read_videofiles.csh
# qsub -l ncpus=1,cput=1:00:00 -q parallel ../../../bin/read_videofiles.csh

if ($?PBS_O_WORKDIR) then
  cd $PBS_O_WORKDIR
endif

date

# Run
src/read_videofiles.x <<EOF
bx
EOF

echo ""
date
