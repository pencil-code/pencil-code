#!/bin/csh
###                       start.csh
###                      -----------
### Run src/start.x (initialising f for src/run.x) with certain
### parameters.
#
# run.csh -- driver for time stepping
#
##PBS -S /bin/csh -W group_list=UK06005 -q UK06005
##PBS -l ncpus=16,mem=8gb,walltime=0:10:00
#PBS -l ncpus=1
#PBS -q p-long
#PBS -l nodes=nq1+nq2
#
if ($?PBS_O_WORKDIR) then
  cd $PBS_O_WORKDIR
endif

# Determine whether this is MPI, how many CPUS etc.
source getconf.csh

#
#  If we don't have a tmp subdirectory: give warning and make a local one
#
if (! -e tmp) then
  echo ""
  echo ">> WARNING: need tmp directory; make local directory underneath"
  echo ">> IN FUTURE: you may want to make a link to some fast scratch disk"
  echo ""
  mkdir tmp
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
# bsub -n  4 -q 4cpu12h -o start.`timestr` -e start.`timestr` start.csh
# bsub -n  8 -q 8cpu12h -o start.`timestr` -e start.`timestr` start.csh
# bsub -n 16 -q 16cpu8h -o start.`timestr` -e start.`timestr` start.csh
