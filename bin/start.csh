#!/bin/csh
###                       start.csh
###                      -----------
### Run src/start.x (initialising f for src/run.x) with certain
### parameters.
#
# run.csh -- driver for time stepping
#
#PBS -S /bin/csh

# cut & paste for job submission for PBS
# qsub -l ncpus=64,mem=32gb,walltime=1:00:00 -W group_list=UK07001 -q UK07001 start.csh
# qsub -l nodes=nq1+nq2+nq3+nq4,mem=1gb,cput=0:10:00 -q p-long start.csh
# qsub -l nodes=4 -q p-long start.csh
# qsub -l nodes=nq1+nq2+nq3+nq4 -q p-long start.csh
# qsub -l nodes=nq1+nq2 -q p-long start.csh
# qsub -l ncpus=4,mem=1gb,walltime=0:05:00 -q parallel start.csh
# qsub -l ncpus=16,mem=1gb,walltime=0:05:00 -q parallel start.csh
#
if ($?PBS_O_WORKDIR) then
  cd $PBS_O_WORKDIR
endif

# Determine whether this is MPI, how many CPUS etc.
source getconf.csh

#
#  If we don't have a tmp subdirectory: stop here (it is too easy to
#  continue with an NFS directory until you fill everything up).
#
if (! -d tmp) then
  echo ""
  echo ">>  STOPPING: need tmp directory"
  echo ">>  Recommended: create tmp as link to directory on a fast scratch disk"
  echo ">>  Not recommended: you can generate tmp with 'mkdir tmp', but that"
  echo ">>  will most likely end up on your NFS file system and be slow"
  echo
  exit 0
endif

# Create list of subdirectories
set subdirs = `printf "%s%s%s\n" "for(i=0;i<$ncpus;i++){" '"tmp/proc";' 'i; }' | bc`
foreach dir ($subdirs)
  # Make sure a sufficient number of subdirectories exist
  if (! -e $dir) then
    mkdir $dir
  else
    # Clean up
    rm -f $dir/VAR* >& /dev/null
    rm -f $dir/vid* >& /dev/null
    rm -f $dir/*.dat >& /dev/null
    rm -f $dir/*.xy $dir/*.xz >& /dev/null
  endif
end
if (-e tmp/n.dat && ! -z tmp/n.dat) mv tmp/n.dat tmp/n.`timestr`
rm -f tmp/*.dat tmp/*.nml tmp/param*.pro tmp/index*.pro >& /dev/null

# Run start.x
date
#
echo "$mpirun $mpirunops $npops src/start.x"
time $mpirun $mpirunops $npops src/start.x
#
echo ""
date

# cut & paste for job submission on the mhd machine
# bsub -n  4 -q 4cpu12h mpijob dmpirun src/start.x
# bsub -n  8 -q 8cpu12h mpijob dmpirun src/start.x
# bsub -n 16 -q 16cpu8h mpijob dmpirun src/start.x
