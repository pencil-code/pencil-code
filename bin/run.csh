#!/bin/csh
# CVS: $Id$

#                       run.csh
#                      ---------
#   Run src/run.x (timestepping for src/run.x).
#   Run parameters are set in run.in.
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

#
#  If necessary, distribute var.dat from the server to the various nodes
#
if ($local_disc) then
  if ($one_local_disc) then     # one common local disc
    foreach node ($nodelist)
      foreach d (`cd $datadir; \ls -d proc* allprocs`)
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

# ---------------------------------------------------------------------- #
rerun:

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
if ($?PBS_JOBID) then
  echo $PBS_JOBID "  RUN STARTED on "$PBS_O_QUEUE `date` \
    >> $datadir/jobid.dat
endif
if ($?LOADL_STEP_ID) then
  echo $LOADL_STEP_ID "  RUN STARTED on "$LOADL_STEP_CLASS `date` \
    >> $datadir/jobid.dat
endif
if ($?SLURM_JOBID) then
  #echo $SLURM_JOBID "  RUN STARTED on "$SLURMD_NODENAME `date` \
  echo $SLURM_JOBID "  RUN STARTED on " `date` \
    >> $datadir/jobid.dat
endif
# EASY job (PDC):
if ($?SP_JID) then
  echo $SP_JID "  RUN STARTED on " `date` >> $datadir/jobid.dat
endif

# Write time and current working directory into log file
(date; echo $cwd; echo "")>> $PENCIL_HOME/.run_directories.log

# Write time string into log file
timestr>> $datadir/runtime.dat

# Run run.x
rm -f ERROR COMPLETED
${PENCIL_HOME}/utils/pc_print_revision_file $run_x
date
echo "$mpirun $mpirunops $npops $mpirunops2 $run_x $x_ops"
echo "#" >> pc_commands.log
date +'# %Y-%m-%d %H:%M:%S' >> pc_commands.log
echo "# $mpirun $mpirunops $npops $mpirunops2 $run_x $x_ops" >> pc_commands.log
time $mpirun $mpirunops $npops $mpirunops2 $run_x $x_ops
#gdb $mpirun $mpirunops $npops $mpirunops2 $run_x $x_ops
set run_status=$status          # save for exit
date

# Create symlinks for deprecated slices
pc_deprecated_slice_links

# Write $PBS_JOBID to file (important when run is migrated within the same job)
if ($?PBS_JOBID) then
  echo $PBS_JOBID " RUN FINISHED on "$PBS_O_QUEUE `date` >> $datadir/jobid.dat
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

# look for NEWDIR file
# if NEWDIR contains a directory name, then continue run in that directory
if (-e "NEWDIR") then
  if (-s "NEWDIR") then
    # Remove LOCK files before going to other directory
    rm -f LOCK data/LOCK
    set olddir=$cwd
    set newdir=`cat NEWDIR`
    touch "$datadir/directory_change.log"
    (echo "stopped run:"; date; echo "new run directory:"; echo $newdir; echo "")\
       >> "$datadir/directory_change.log"
    cd "$newdir"
    touch "$datadir/directory_change.log"
    (date; echo "original run script is in:"; echo $olddir; echo "")\
       >> "$datadir/directory_change.log"
    echo
    echo "====================================================================="
    echo "Rerunning in new directory:"
    pwd
    echo "Current status: $run_status"
    echo "====================================================================="
    echo
    rm "$olddir/NEWDIR"
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

# Shut down mpd if we have started it
if ($?booted_mpd) then
  echo "Shutting down mpd .."
  mpdallexit
  echo "..done"
else
endif

# remove LOCK files
rm -f LOCK data/LOCK

# Detect error status flagged by code (for cases where this does not get
# propagated to the mpirun status):
if (-e 'ERROR') then
  echo "Found ERROR file from run.x"
  set run_status2 = 16
else
  set run_status2 = 0
  if(-e "RESUBMIT") then
    rm -f rs
    echo "Resubmitting from node:$hn to masternode $masterhost"
    echo "$resub $resubop" > rs
    chmod +x rs
    if(-e "resubmit.log") then
      $run_resub >>resubmit.log
      date >>resubmit.log
      echo "=====================" >> resubmit.log
    else
      $run_resub >resubmit.log
      date >>resubmit.log
      echo "=====================" >> resubmit.log
    endif
  else
  endif
endif

exit ( $run_status | $run_status2 ) # propagate status of mpirun


# cut & paste for job submission on the mhd machine
# bsub -n  4 -q 4cpu12h -o run.`timestr` -e run.`timestr` run.csh
# bsub -n  8 -q 8cpu12h mpijob dmpirun $run_x
# bsub -n 16 -q 16cpu8h mpijob dmpirun $run_x
# bsub -n  8 -q 8cpu12h -o run.log -w 'exit(123456)' mpijob dmpirun $run_x

# qsub -l ncpus=64,mem=32gb,walltime=500:00:00 -W group_list=UK07001 -q UK07001 run.csh
# qsub -l nodes=4:ppn=1,mem=500mb,cput=24:00:00 -q p-long run.csh
# qsub -l ncpus=16,mem=1gb,cput=400:00:00 -q parallel run.csh
# qsub -l nodes=128,walltime=10:00:00 -q workq run.csh
# qsub -l nodes=1:ppn=4 run.csh
# eval `env-setup lam`; qsub -v PATH -pe lam 8 -j y -o run.log run.csh
