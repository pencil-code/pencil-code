#!/bin/csh
#                       start_run.csh
#                      ---------------
#   Run first src/start.x and then src/run.x.
#   On machines with local scratch disks this saves the overhead
#   associated with copying to and from that disk.
#
# Run this script with csh:
#PBS -S /bin/csh
#$ -S /bin/csh
#@$-s /bin/csh
#
# Join stderr and stout:
#$ -j y -o run.log
#@$-eo
#
# Work in submit directory (SGE):
#$ -cwd

# Work in submit directory (PBS):
if ($?PBS_O_WORKDIR) then
  cd $PBS_O_WORKDIR
endif

# Work in submit directory (SUPER-UX's nqs):
if ($?QSUB_WORKDIR) then
  cd $QSUB_WORKDIR
endif

# ---------------------------------------------------------------------- #

# Common setup for start.csh, run.csh, start_run.csh:
# Determine whether this is MPI, how many CPUS etc.
source getconf.csh

# Prevent code from running twice (and removing files by accident)
if (! -e "NEVERLOCK") touch LOCK

#
#  If we don't have a data subdirectory: stop here (it is too easy to
#  continue with an NFS directory until you fill everything up).
#
if (! -d "$datadir") then
  echo ""
  echo ">>  STOPPING: need $datadir directory"
  echo ">>  Recommended: create $datadir as link to directory on a fast scratch"
  echo ">>  Not recommended: you can generate $datadir with 'mkdir $datadir', "
  echo ">>  but that will most likely end up on your NFS file system and be"
  echo ">>  slow"
  echo
  rm -f LOCK
  exit 0
endif

# Create list of subdirectories
set subdirs = `printf "%s%s%s\n" "for(i=0;i<$ncpus;i++){" '"data/proc";' 'i; }' | bc`
foreach dir ($subdirs)
  # Make sure a sufficient number of subdirectories exist
  if (! -e $dir) then
    mkdir $dir
  else
    # Clean up
    # when used with lnowrite=T, for example, we don't want to remove var.dat:
    set list=`/bin/ls $dir/VAR* $dir/*.dat $dir/*.info $dir/slice*`
    #if ($list != "") then
      foreach rmfile ($list)
        if ($rmfile != $dir/var.dat) rm -f $rmfile >& /dev/null
      end
    #endif
  endif
end
if (-e $datadir/time_series.dat && ! -z $datadir/time_series.dat) mv $datadir/time_series.dat $datadir/time_series.`timestr`
rm -f $datadir/*.dat $datadir/*.nml $datadir/param*.pro $datadir/index*.pro >& /dev/null

# If local disk is used, copy executable to $SCRATCH_DIR of master node
if ($local_binary) then
  cp src/start.x $SCRATCH_DIR
  remote-top >& remote-top.log &
endif

# Run start.x
date
echo "$mpirun $mpirunops $npops $start_x $x_ops"
time $mpirun $mpirunops $npops $start_x $x_ops
echo ""
date

# Clean up control and data files
rm -f STOP RELOAD fort.20

# On machines with local scratch directory, initialize automatic
# background copying of snapshots back to the data directory.
# Also, if necessary copy executable to $SCRATCH_DIR of master node
# and start top command on all procs
if ($local_disc) then
  echo "Use local scratch disk"
  copy-snapshots -v >& copy-snapshots.log &
  remote-top >& remote-top.log &
endif
if ($local_binary) then
  echo "ls src/run.x $SCRATCH_DIR before copying:"
  ls -lt src/run.x $SCRATCH_DIR
  cp src/run.x $SCRATCH_DIR
  echo "ls src/run.x $SCRATCH_DIR after copying:"
  ls -lt src/run.x $SCRATCH_DIR
endif

# Run run.x
date
echo "$mpirun $mpirunops $npops $run_x $x_ops"
echo $mpirun $mpirunops $npops $run_x $x_ops >! run_command.log
time $mpirun $mpirunops $npops $run_x $x_ops
date

# On Horseshoe cluster, copy var.dat back to the data directory
if ($local_disc) then
  echo "Use local scratch disk"
  copy-snapshots -v var.dat
  echo "done, will now killall copy-snapshots"
  killall copy-snapshots
endif

# Shut down lam if we have started it
if ($booted_lam) lamhalt

# remove LOCK file
if (-e "LOCK") rm -f LOCK

# cut & paste for job submission on the mhd machine
# bsub -n  4 -q 4cpu12h mpijob dmpirun src/start.x
# bsub -n  8 -q 8cpu12h mpijob dmpirun src/start.x
# bsub -n 16 -q 16cpu8h mpijob dmpirun src/start.x

# cut & paste for job submission for PBS
# qsub -l ncpus=64,mem=32gb,walltime=1:00:00 -W group_list=UK07001 -q UK07001 start.csh
# qsub -l nodes=4:ppn=1,mem=500mb,cput=24:00:00 -q p-long start.csh
# qsub -l ncpus=4,mem=1gb,cput=100:00:00 -q parallel start.csh
# qsub -l nodes=128,mem=64gb,walltime=1:00:00 -q workq start.csh
