#!/bin/csh
#
# run.csh -- driver for time stepping
#
#PBS -S /bin/csh -W group_list=UK03007 -q UK03007
#PBS -l ncpus=16,mem=8gb,walltime=90:00:00
#
if ($?PBS_O_WORKDIR) then
  cd $PBS_O_WORKDIR
endif

# Determine whether this is MPI, how many CPUS etc.
source getconf.csh

# Clean up control and data files
rm -f STOP RELOAD fort.20

# Run run.x
date
#
echo "$mpirun $mpirunops $npops src/run.x"
time $mpirun $mpirunops $npops src/run.x
#
date

cat fort.20 >> tmp/n.dat

# bsub -n 4 -q 4cpu4h -o irro.`timestr` -e irro.`timestr` irro.csh
