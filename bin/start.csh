#!/bin/csh
# CVS: $Id$

#                       start.csh
#                      -----------
#   Run src/start.x (initialising f for src/run.x).
#   Start parameters are set in start.in.
#
# Run this script with csh:
#PBS -S /bin/csh
#$ -S /bin/csh
#@$-s /bin/csh
#
# Work in submit directory (SGE):
#$ -cwd -V
# Work in submit directory (PBS):
if ($?PBS_O_WORKDIR) then
  cd $PBS_O_WORKDIR
endif

# Work in submit directory (SUPER-UX's nqs):
if ($?QSUB_WORKDIR) then
  cd $QSUB_WORKDIR
endif

# Work in submit directory (IBM Loadleveler):
if ($?LOADL_STEP_INITDIR) then
  cd $LOADL_STEP_INITDIR
endif

# Common setup for start.csh, run.csh, start_run.csh:
# Determine whether this is MPI, how many CPUS etc.
pwd
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
  if (-e "LOCK") rm -f LOCK
  if (-e "data/LOCK") rm -f data/LOCK
  exit 0
endif

#
#  Execute some script or command specified by the user.
#  E.g. run src-to-data to save disk space in /home by moving src/ to
#  data/ and linking
#
if ($?PENCIL_START1_CMD) then
  echo "Running $PENCIL_START1_CMD"
  $PENCIL_START1_CMD
endif

# ---------------------------------------------------------------------- #

#  For testing backwards compatibility, do not excecute start.x, an do not
#  delete existing var.dat (etc.) files.
#  Rename time_series.dat, so run.x can write a fresh one to compare
#  against.
#
if (-e NOSTART) then
  echo "Found NOSTART file. Won't run start.x"
  if (-e $datadir/time_series.dat) \
      mv $datadir/time_series.dat $datadir/time_series.`timestr`
  exit
endif

# ---------------------------------------------------------------------- #

# Create list of subdirectories
# If the file NOERASE exists, the old directories are not erased
#   (start.x also knows then that var.dat is not created)
foreach dir ($procdirs $subdirs)
  # Make sure a sufficient number of subdirectories exist
  set ddir = "$datadir/$dir"
  if (! -e $ddir) then
    mkdir $ddir
  else if (! -e NOERASE) then
    # Clean up
    rm -f $datadir/../driver/pts_[1-9].dat $datadir/../driver/seed_[1-9].dat >& /dev/null
    if ($dir == "allprocs") then
      rm -f $ddir/VAR[0-9]* $ddir/PERS[0-9]* $ddir/grid.dat $ddir/dim.dat $ddir/varN.list >& /dev/null
    else
      rm -f $ddir/VAR[0-9]* $ddir/PERS[0-9]* $ddir/TAVG[0-9]* $ddir/*.info $ddir/slice* $ddir/PVAR[0-9]* $ddir/SPVAR[0-9]* $ddir/varN.list >& /dev/null
      # in some cases var.dat needs to be conserved (eg. lnowrite=T)
      set list = `/bin/ls $ddir/*.dat`
      foreach rmfile ($list)
        if ($rmfile != $ddir/var.dat && $rmfile != $ddir/pers.dat) rm -f $rmfile >& /dev/null
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
  if ($lcopysnapshots_exp) rm -f $datadir/move-me.list $datadir/moved-files.list >& /dev/null
  rm -f ioerrors.log >& /dev/null
endif

# If local_binary is used, copy executable to $SCRATCH_DIR of master node
if ($local_binary) then
  echo "Copying start.x to $SCRATCH_DIR"
  cp src/start.x $SCRATCH_DIR
  remote-top >& remote-top.log &
endif

# Run start.x
rm -f 'ERROR'
date
echo "$mpirun $mpirunops $npops $mpirunops2 $start_x $x_ops"
time $mpirun $mpirunops $npops $mpirunops2 $start_x $x_ops
set start_status=$status        # save for exit
echo ""
date

# Not sure it makes any sense to continue after mpirun had an error:
if ($start_status) then
  echo "Error status $start_status found -- aborting"
  exit $start_status
endif


# If local disk is used, copy var.dat back to the data directory
if ($local_disc) then
  echo "Copying var.dat etc. back to data directory"
  $copysnapshots -v var.dat     >&! copy-snapshots.log
  if (-e $SCRATCH_DIR/proc0/VAR0) $copysnapshots -v VAR0 >&! copy-snapshots.log
  if ($lparticles) $copysnapshots -v pvar.dat >>& copy-snapshots.log
  if ($lparticles_nbody) $copysnapshots -v spvar.dat >>& copy-snapshots.log
  $copysnapshots -v global.dat  >>& copy-snapshots.log
  $copysnapshots -v timeavg.dat >>& copy-snapshots.log
  $copysnapshots -v dxyz.dat    >>& copy-snapshots.log

  if ($remove_scratch_root) then
    rm -rf $SCRATCH_DIR
  endif
endif

# Shut down mpd if we have started it
if ($?booted_mpd) then
  echo "Shutting down mpd .."
  mpdallexit
  echo "..done"
else
endif

# remove LOCK file
if (-e "LOCK") rm -f LOCK
if (-e "data/LOCK") rm -f data/LOCK

# Detect error status flagged by code (for cases where this does not get
# propagated to the mpirun status):
if (-e 'ERROR') then
  echo "Found ERROR file from start.x"
  set start_status2 = 16
else
  set start_status2 = 0
endif

rm -f resubmit.log
rm -f rs

exit ( $start_status | $start_status2 )        # propagate status of mpirun


# cut & paste for job submission on the mhd machine
# bsub -n  4 -q 4cpu12h mpijob dmpirun src/start.x
# bsub -n  8 -q 8cpu12h mpijob dmpirun src/start.x
# bsub -n 16 -q 16cpu8h mpijob dmpirun src/start.x

# cut & paste for job submission for PBS
# qsub -l ncpus=64,mem=32gb,walltime=1:00:00 -W group_list=UK07001 -q UK07001 start.csh
# qsub -l nodes=4:ppn=1,mem=500mb,cput=24:00:00 -q p-long start.csh
# qsub -l ncpus=4,mem=1gb,cput=100:00:00 -q parallel start.csh
# qsub -l nodes=128,mem=64gb,walltime=1:00:00 -q workq start.csh
# For Andromeda
# qsub $QSUB_OPTIONS $qstring   -pe $PE $nodes $FILE
