#!/bin/csh
###                       start.csh
###                      -----------
### Run src/start.x (initialising f for src/run.x) with certain
### parameters.

# Determine whether this is MPI, how many CPUS etc.
source getconf.csh

# Create list of subdirectories
set subdirs = `printf "%s%s%s\n" "for(i=0;i<$ncpus;i++){" '"tmp/proc";' 'i; }' | bc`
foreach dir ($subdirs)
  # Make sure a sufficient number of subdirectories exist
  if (! -e $dir) then
    mkdir $dir
  else
    # Clean up
    rm -f tsnap.dat tvid.dat >& /dev/null
    rm -f $dir/VAR* >& /dev/null
    rm -f $dir/vid* >& /dev/null
  endif
end
rm -f tmp/n.dat >& /dev/null

# Run start.x
date
#
echo "$mpirun $mpirunops $npops src/start.x"
time $mpirun $mpirunops $npops src/start.x
#
date

#Save initial state
foreach dir ($subdirs)
  cp $dir/var.dat $dir/VAR0
end

