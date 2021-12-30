#!/bin/csh
# CVS: $Id$

#                       start_run.csh
#                      ---------------
#   Run first src/start.x (unless data/proc0/var.dat already exists) and
#   then src/run.x.
#   On machines with local scratch disks this saves the overhead
#   associated with copying to and from that disk.
#
# Run this script with csh:
#PBS -S /bin/csh
#PBS -r n
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

# Work in submit directory (IBM Loadleveler):
if ($?LOADL_STEP_INITDIR) then
  cd $LOADL_STEP_INITDIR
endif

# ====================================================================== #
# Starting points when rerunning in new directory.
# This needs to come *before* "source getconf.csh", because
#  (i) it checks for LOCK file, and
# (ii) it sets $datadir/directory_snap
newdir:

# Common setup for start.csh, run.csh, start_run.csh:
# Determine whether this is MPI, how many CPUS etc.
source getconf.csh

# ====================================================================== #
#  Assume that if data/proc0/var.dat exists, we are not supposed
#  to overwrite it, and should try to restart from the existing
#  snapshot instead. In that case, if $local_disc was set in getconf.csh,
#  we need to copy these data first to the locak scratch disc, and
#  then restart, i.e. jump to rerun.
#
if (-e "$datadir"/proc0/var.dat) then
#
#  If necessary, distribute var.dat from the server to the various nodes
#  Don't indent these lines so that it is easier to vimdiff against run.csh
#
  if ($local_disc) then
    if ($one_local_disc) then     # one common local disc
      foreach node ($nodelist)
        foreach d (`cd $datadir; ls -d proc* allprocs`)
          if (-e $datadir/$d/var.dat) $SCP $datadir/$d/var.dat ${node}:$SCRATCH_DIR/$d/
          if (-e $datadir/$d/global.dat) $SCP $datadir/$d/global.dat ${node}:$SCRATCH_DIR/$d/
          if ($lparticles) $SCP $datadir/$d/pvar.dat ${node}:$SCRATCH_DIR/$d/
          if ($lpointmasses) $SCP $datadir/$d/qvar.dat ${node}:$SCRATCH_DIR/$d/
          $SCP $datadir/$d/timeavg.dat ${node}:$SCRATCH_DIR/$d/
        end
        if (-e $datadir/allprocs/dxyz.dat) $SCP $datadir/allprocs/dxyz.dat ${node}:$SCRATCH_DIR/allprocs
      end
    else # one local disc per MPI process (Horseshoe, etc);
         # still doesn't cover Copson
      set i = -1
      foreach node ($nodelist)
        set i=`expr $i + 1`
        echo "i = $i"
        set j = 0
        while ($j != $nprocpernode)
          set k = `expr $nprocpernode \* $i + $j`
          if ($?notserial_procN) set k = `expr $i + $nnodes \* $j`
          $SCP $datadir/proc$k/var.dat ${node}:$SCRATCH_DIR/proc$k/
          if (-e $datadir/proc$k/global.dat) then
            $SCP $datadir/proc$k/global.dat ${node}:$SCRATCH_DIR/proc$k/
          endif
          if ($lparticles) then
            $SCP $datadir/proc$k/pvar.dat ${node}:$SCRATCH_DIR/proc$k/
          endif
          if ($lpointmasses) then
            $SCP $datadir/proc$k/qvar.dat ${node}:$SCRATCH_DIR/proc$k/
          endif
          echo "$SCP $datadir/proc$k/var.dat ${node}:$SCRATCH_DIR/proc$k/"
          if (-e $datadir/proc$k/timeavg.dat) then
            $SCP $datadir/proc$k/timeavg.dat ${node}:$SCRATCH_DIR/proc$k/
          endif
          set j=`expr $j + 1`
        end
        if (-e $datadir/allprocs/dxyz.dat) then
          $SCP $datadir/allprocs/dxyz.dat ${node}:$SCRATCH_DIR/allprocs/
        endif
      end
    endif
  endif
#
# goto rerun instead of producing a new initial condition
#
  goto rerun
endif

#AB: the following 4 lines prevent me from jobtransferring into a fresh directory
# First, check if NEWDIR is set and if yes, go there.
#if (-e "NEWDIR") then
#  goto checknewdir
#endif

# If we don't have a data subdirectory: create it, if set to default.
if ((! -d "$datadir") && ("$datadir" == "data")) then
  echo "Creating default datadir: '$datadir'"
  mkdir "$datadir"
endif

#  If we still don't have a data subdirectory: stop here (it is too easy to
#  continue with an NFS directory until you fill everything up).
if (! -d "$datadir") then
  echo ""
  echo ">>  STOPPING: need $datadir directory"
  echo ">>  Recommended: create $datadir as link to directory on a fast scratch"
  echo ">>  Not recommended: you can generate $datadir with 'mkdir $datadir', "
  echo ">>  but that will most likely end up on your NFS file system and be"
  echo ">>  slow"
  echo
  rm -f LOCK data/LOCK IO_LOCK
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
    if ($dir == "allprocs" || $dir == "reduced") then
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
  rm -f $datadir/../driver/pts_[1-9].dat $datadir/../driver/seed_[1-9].dat >& /dev/null
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
  cp $start_x $SCRATCH_DIR
  remote-top >& remote-top.log &
endif

# Run start.x
rm -f ERROR COMPLETED IO_LOCK
${PENCIL_HOME}/utils/pc_print_revision_file $start_x
date
echo "$mpirun $mpirunops $npops $mpirunops2 $start_x $x_ops"
echo "" >> pc_commands.log
date +'# %Y-%m-%d %H:%M:%S' >> pc_commands.log
echo "$mpirun $mpirunops $npops $mpirunops2 $start_x $x_ops" >> pc_commands.log
time $mpirun $mpirunops $npops $mpirunops2 $start_x $x_ops
set start_status=$status        # save for exit
echo ""
date

# Not sure it makes any sense to continue after mpirun had an error:
if ($start_status) then
  echo "Error status $start_status found -- aborting"
  exit $start_status
endif

# Detect error status flagged by code (for cases where this does not get
# propagated to the mpirun status):
if (-e 'ERROR') then
  echo "Found ERROR file from start.x"
  set start_status2 = 16
else
  set start_status2 = 0
endif

if ($start_status2) exit $start_status2 # propagate error status

# ---------------------------------------------------------------------- #
rerun:

# if the file ADDITIONAL_RUN_COMMAND.csh exist, execute it.
# common file content could be:
# copy-proc-to-proc var.dat . ../another_run_directory
# if one wants to use the var.dat files from another_run_directory
# after having said start.x
if (-e ADDITIONAL_RUN_COMMAND.csh) then
  source ADDITIONAL_RUN_COMMAND.csh
endif


# Clean up control and data files
# NB. Don't remove NEWDIR it may have been put there on purpose so as
#     to catch a crash and run something else instead.
rm -f STOP RELOAD RERUN fort.20

# On machines with local scratch directory, initialize automatic
# background copying of snapshots back to the data directory.
# Also, if necessary copy executable to $SCRATCH_DIR of master node
# and start top command on all procs.
if ($local_disc) then
  echo "Use local scratch disk"
  $copysnapshots -v >&! copy-snapshots.log &
endif
# Copy output from `top' on run host to a file we can read from login server
if ($remote_top) then
  remote-top >&! remote-top.log &
endif
if ($local_binary) then
  echo "ls $run_x $SCRATCH_DIR before copying:"
  ls -lt $run_x $SCRATCH_DIR
  cp $run_x $SCRATCH_DIR
  echo "ls $run_x $SCRATCH_DIR after copying:"
  ls -lt $run_x $SCRATCH_DIR
endif

# Write $PBS_JOBID or $LOADL_STEP_ID to file
# (important when run is migrated within the same job)
#if ($?PBS_JOBID) then
#  echo $PBS_JOBID "  # RUN STARTED on "$PBS_O_QUEUE `date` \
#    >> $datadir/jobid.dat
#endif
if ($?LOADL_STEP_ID) then
  echo $LOADL_STEP_ID " # RUN STARTED on "$LOADL_STEP_CLASS `date` \
    >> $datadir/jobid.dat
endif
if ($?SLURM_JOB_ID) then
  echo $SLURM_JOB_ID " # RUN STARTED on "$SLURMD_NODENAME `date` \
    >> $datadir/jobid.dat
endif
# EASY job (PDC):
if ($?SP_JID) then
  echo $SP_JID " # RUN STARTED on " `date` >> $datadir/jobid.dat
endif

# Write time and current working directory into log file
(date; echo $cwd; echo "")>> $PENCIL_HOME/.run_directories.log

# Write time string into log file
timestr>> $datadir/runtime.dat

# Run run.x
rm -f ERROR COMPLETED
date
${PENCIL_HOME}/utils/pc_print_revision_file $run_x
echo "$mpirun $mpirunops $npops $mpirunops2 $run_x $x_ops"
echo "" >> pc_commands.log
date +'# %Y-%m-%d %H:%M:%S' >> pc_commands.log
echo "$mpirun $mpirunops $npops $mpirunops2 $run_x $x_ops" >> pc_commands.log
time $mpirun $mpirunops $npops $mpirunops2 $run_x $x_ops
set run_status=$status          # save for exit
date

# Create symlinks for deprecated slices
pc_deprecated_slice_links

# Write $PBS_JOBID to file (important when run is migrated within the same job)
#if ($?PBS_JOBID) then
#  echo $PBS_JOBID " # RUN FINISHED on "$PBS_O_QUEUE `date` >> $datadir/jobid.dat
#endif
if ($?SLURM_JOB_ID) then
  echo $SLURM_JOB_ID " # RUN FINISHED on "$SLURMD_NODENAME `date` >> $datadir/jobid.dat
endif
# EASY job (PDC):
if ($?SP_JID) then
  echo $SP_JID " # RUN FINISHED on " `date` >> $datadir/jobid.dat
endif

# look for RERUN file
# With this method one can only reload a new executable.
# One cannot change directory, nor are the var.dat files returned to server.
# See the NEWDIR method below for more options.
if (-e "RERUN") then
  rm -f RERUN
  echo
  echo "======================================================================="
  echo "Rerunning in the *same* directory; current run status: $run_status"
  echo "We are *still* in: " `pwd`
  echo "======================================================================="
  echo
  goto rerun
endif
# ---------------------------------------------------------------------- #

# On machines with local scratch disc, copy var.dat back to the data
# directory
if ($local_disc) then
  echo "Copying all var.dat, VAR*, TIMEAVG*, dxyz.dat, timeavg.dat and crash.dat back from local scratch disks"
  $copysnapshots -v var.dat     >&! copy-snapshots2.log
  if ($lparticles) $copysnapshots -v pvar.dat >>& copy-snapshots2.log
  if ($lpointmasses) $copysnapshots -v qvar.dat >>& copy-snapshots2.log
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
    rm -rf $SCRATCH_DIR
  endif
endif

echo "Done"

checknewdir:

# look for NEWDIR file
# if NEWDIR contains a directory name, then continue run in that directory
if (-e "NEWDIR") then
  if (-s "NEWDIR") then
    # Remove LOCK files before going to other directory, because we are done in the old one.
    rm -f LOCK data/LOCK
    set olddir=$cwd
    cd `cat NEWDIR`
    rm $olddir/NEWDIR
    (echo "stopped run:"; date; echo "new run directory:"; echo $cwd; echo "")\
       >> $olddir/$datadir/directory_change.log
    (date; echo "original run script is in:"; echo $olddir; echo "")\
       >> $datadir/directory_change.log
    echo
    echo "====================================================================="
    echo "Rerunning in new directory; current run status: $run_status"
    echo "We are now in: " `pwd`
    echo "====================================================================="
    echo
    goto newdir
  else
    rm -f NEWDIR LOCK data/LOCK
    echo
    echo "====================================================================="
    echo "Rerunning in the *same* directory; current run status: $run_status"
    echo "We are *still* in: " `pwd`
    echo "====================================================================="
    echo
    echo "Rerunning; current run status: $run_status"
    goto newdir
  endif
endif
# ====================================================================== #

# Shut down lam if we have started it
if ($booted_lam) lamhalt

# remove LOCK files
rm -f LOCK data/LOCK

# Detect error status flagged by code (for cases where this does not get
# propagated to the mpirun status):
if (-e 'ERROR') then
  echo "Found ERROR file from run.x"
  set run_status2 = 16
else
  set run_status2 = 0
endif

exit ( $run_status | $run_status2 ) # propagate status of mpirun

# cut & paste for job submission on the mhd machine
# bsub -n  4 -q 4cpu12h mpijob dmpirun $start_x
# bsub -n  4 -q 4cpu12h mpijob dmpirun $run_x
# bsub -n  8 -q 8cpu12h mpijob dmpirun $start_x
# bsub -n  8 -q 8cpu12h mpijob dmpirun $run_x
# bsub -n 16 -q 16cpu8h mpijob dmpirun $start_x
# bsub -n 16 -q 16cpu8h mpijob dmpirun $run_x

# cut & paste for job submission for PBS
# qsub -l ncpus=64,mem=32gb,walltime=1:00:00 -W group_list=UK07001 -q UK07001 start_run.csh
# qsub -l nodes=4:ppn=1,mem=500mb,cput=24:00:00 -q p-long start_run.csh
# qsub -l ncpus=4,mem=1gb,cput=100:00:00 -q parallel start_run.csh
# qsub -l nodes=128,mem=64gb,walltime=1:00:00 -q workq start_run.csh
# eval `env-setup lam`; qsub -v PATH -pe lam 8 -j y -o run.log run.csh
# cut & paste for job submission for SGE with parallel environment openmpi (e.g. andromeda)
qsub -pe openmpi 1 ./start_run.csh
qsub -pe openmpi 8 ./start_run.csh
