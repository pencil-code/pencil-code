#!/bin/csh
# CVS: $Id: start_run.csh,v 1.41 2006-12-06 19:53:17 wlyra Exp $

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
# If the file NOERASE exists, the old directories are not erased
#   (start.x also knows then that var.dat is not created)
foreach dir ($procdirs $subdirs)
  # Make sure a sufficient number of subdirectories exist
  set ddir = "$datadir/$dir"
  if (! -e $ddir) then
    mkdir $ddir
  else
    # Clean up
    # when used with lnowrite=T, for example, we don't want to remove var.dat:
    set list = \
        `/bin/ls $ddir/VAR* $ddir/TAVG* $ddir/*.dat $ddir/*.info $ddir/slice* $ddir/PVAR*`
    if (! -e NOERASE) then
      foreach rmfile ($list)
        if ($rmfile != $ddir/var.dat) rm -f $rmfile >& /dev/null
      end
    endif
  endif
end

# Clean up previous runs
if (! -e NOERASE) then
  if (-e $datadir/time_series.dat && ! -z $datadir/time_series.dat) \
      mv $datadir/time_series.dat $datadir/time_series.`timestr`
  rm -f $datadir/*.dat $datadir/*.nml $datadir/param*.pro $datadir/index*.pro \
        $datadir/averages/* >& /dev/null
  rm -f $datadir/move-me.list $datadir/moved-files.list >& /dev/null
endif

# If local disk is used, copy executable to $SCRATCH_DIR of master node
if ($local_binary) then
  cp src/start.x $SCRATCH_DIR
endif
# Copy output from `top' on run host to a file we can read from login server
if ($remote_top) then
  remote-top >& remote-top.log &
endif

# Run start.x
date
echo "$mpirun $mpirunops $npops $mpirunops2 $start_x $x_ops"
time $mpirun $mpirunops $npops $mpirunops2 $start_x $x_ops
set start_status=$status	# save for exit
if ($start_status) exit $start_status	# something went wrong
echo ""
date

# Clean up control and data files
rm -f STOP RELOAD fort.20

if ($local_disc) then
  # We still need to copy (at least one of) the var.dat files back, so
  # the background process copy-snapshots will know how large the snapshots
  # ought to be. Certainly far from elegant..
  $copysnapshots -v var.dat     >&  copy-snapshots.log
  if ($lparticles) $copysnapshots -v pvar.dat >>& copy-snapshots.log
  $copysnapshots -v timeavg.dat >>& copy-snapshots.log

  # On machines with local scratch directory, initialize automatic
  # background copying of snapshots back to the data directory.
  # Also, if necessary copy executable to $SCRATCH_DIR of master node
  # and start top command on all procs.
  echo "Use local scratch disk"
  $copysnapshots -v >>& copy-snapshots.log &
endif
# Copy output from `top' on run host to a file we can read from login server
if ($remote_top) then
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
echo "$mpirun $mpirunops $npops $mpirunops2 $run_x $x_ops"
echo $mpirun $mpirunops $npops $mpirunops2 $run_x $x_ops >! run_command.log
time $mpirun $mpirunops $npops $mpirunops2 $run_x $x_ops
set run_status=$status		# save for exit
date

# On machines with local scratch disc, copy var.dat back to the data
# directory
if ($local_disc) then
  echo "Copying all var.dat, VAR*, TIMEAVG*, dxyz.dat, timeavg.dat and crash.dat back from local scratch disks"
  $copysnapshots -v var.dat     >&! copy-snapshots2.log
  if ($lparticles) $copysnapshots -v pvar.dat >>& copy-snapshots2.log
  $copysnapshots -v -1          >>& copy-snapshots2.log
  $copysnapshots -v dxyz.dat    >>& copy-snapshots2.log
  $copysnapshots -v timeavg.dat >>& copy-snapshots2.log
  $copysnapshots -v crash.dat   >>& copy-snapshots2.log
  echo "done, will now killall copy-snapshots"
  # killall copy-snapshots   # Linux-specific
  set pids=`ps -U $USER -o pid,command | grep -E 'remote-top|copy-snapshots' | sed 's/^ *//' | cut -d ' ' -f 1`
  echo "Shutting down processes ${pids}:"
  foreach p ($pids)  # need to do in a loop, and check for existence, since
                     # some systems (Hitachi) abort this script when trying
                     # to kill non-existent processes
    echo "  pid $p"
    if ( `ps -p $p | fgrep -c $p` ) kill -KILL $p
  end

  if ($remove_scratch_root) then
    rm -rf $SCRATCH_DIR/*
    rm -rf $SCRATCH_DIR
  endif
endif
echo "Done"

# Shut down lam if we have started it
if ($booted_lam) lamhalt

# remove LOCK file
if (-e "LOCK") rm -f LOCK

# look for RERUN file 
if (-e "RERUN") then 
  if (-s "RERUN") then
    cd `cat RERUN`
  endif

  rm -f RERUN

  echo "Rerunning; current run status: $run_status"
  #  goto rerun
  ./run.csh
endif  

exit $run_status		# propagate status of mpirun

# cut & paste for job submission on the mhd machine
# bsub -n  4 -q 4cpu12h mpijob dmpirun src/start_run.x
# bsub -n  8 -q 8cpu12h mpijob dmpirun src/start_run.x
# bsub -n 16 -q 16cpu8h mpijob dmpirun src/start_run.x

# cut & paste for job submission for PBS
# qsub -l ncpus=64,mem=32gb,walltime=1:00:00 -W group_list=UK07001 -q UK07001 start_run.csh
# qsub -l nodes=4:ppn=1,mem=500mb,cput=24:00:00 -q p-long start_run.csh
# qsub -l ncpus=4,mem=1gb,cput=100:00:00 -q parallel start_run.csh
# qsub -l nodes=128,mem=64gb,walltime=1:00:00 -q workq start_run.csh
