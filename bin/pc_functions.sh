#!/bin/sh

match() {
  if perl -e "if ( \"$1\" =~ /$2/ ) { exit(0) } else { exit (1) }"
  then
    return 0;
  else
    return 1;
  fi
}


if match "$0" "pc_functions.sh" ; then
  echo "----------------------- pc_functions.sh ------------------------------ "
  echo " Now I don't believe you wanted to do that now did you?                "
  echo "                                                                       "
  echo " This script is designed to be sourced by other scripts to provide     "
  echo " shell functions and definitions useful to Pencil-Code scripts.        "
  echo "                                                                       "
  echo " Since I'm not intended to be sourced interactively and don't do       "
  echo " anything when executed I'm going to die now.  Good bye.....           "
  echo "---------------------------------------------------------------------- "
  unset match
  kill $$  
fi

#------------------------------------------------------------------------------
#-- Information Routines                        
#------------------------------------------------------------------------------
show_config()
{
  # Display all relevant configuration information
  echo '$mpi            = ' "<$mpi>"
  echo '$ncpus          = ' "<$ncpus>"
  echo '$npops          = ' "<$npops>"
  echo '$local_disc     = ' "<$local_disc>"
  echo '$one_local_disc = ' "<$one_local_disc>"
  echo '$remote_top   	= ' "<$remote_top>"
  echo '$local_binary 	= ' "<$local_binary>"
  echo '$datadir      	= ' "<$datadir>"
  echo '$SCRATCH_DIR  	= ' "<$SCRATCH_DIR>"
  echo '$hn           	= ' "<$hn>"
  echo '$start_x       	= ' "<$start_x>"
  echo '$run_x       	= ' "<$run_x>"
  echo '$mpirun       	= ' "<$mpirun>"
  echo '$mpirunops    	= ' "<$mpirunops>"
  echo '$x_ops        	= ' "<$x_ops>"
  echo '$NODELIST     	= ' "<$NODELIST>"
  echo '$SSH          	= ' "<$SSH>"
  echo '$SCP          	= ' "<$SCP>"
  echo '$PARENT_PID   	= ' "<$PARENT_PID>"
}

#------------------------------------------------------------------------------
#-- Scratch Directory Routines                        
#------------------------------------------------------------------------------
prepare_scratch_disk()
{
  # If local disc is used, write name into $datadir/directory_snap.
  # This will be read by the code, if the file exists.
  # Remove file, if not needed, to avoid confusion.
  [ -f $datadir/directory_snap ] && rm $datadir/directory_snap
  [ "$local_disc" == "yes" ] &&  echo $SCRATCH_DIR >$datadir/directory_snap


  if [ "$local_binary" == "yes" ]; then
    start_x=$SCRATCH_DIR/start.x
    run_x=$SCRATCH_DIR/run.x
  fi

  # Created subdirectories on local scratch disc (start.csh will also create
  # them under $datadir/)
  if [ "$local_disc" == "yes" ]; then
    echo "Creating directory structure on scratch disc(s)"
    for host in $nodelist 
    do
      $SSH $host "if (! -e $SCRATCH_DIR ) mkdir -p $SCRATCH_DIR; cd $SCRATCH_DIR; mkdir -p $procdirs $subdirs" 
    done
  fi
}

tidy_scratch_disk()
{
  # Tidy scratch directories
  if [ "$local_disc" == "yes" ]; then
    if [ $remove_scratch_root ]; then
      echo -n "Tidying scratch directory ($SCRATCH_DIR)...   "
      rm -rf $SCRATCH_DIR/*
      rm -rf $SCRATCH_DIR
      echo "Done."
    fi
  fi
}


#------------------------------------------------------------------------------
#-- Data Directory Routines                        
#------------------------------------------------------------------------------
determine_datadir()
{
  # Determine data directory (defaults to `data')
  if [ -r datadir.in ]; then
    datadir=`cat datadir.in | sed 's/ *\([^ ]*\).*/\1/'`
  else
    datadir=data
  fi
  echo "datadir = $datadir"
  subdirs=("allprocs" "averages" "idl")
  procdirs=`printf "%s%s%s\n" "for(i=0;i<$ncpus;i++){" '"proc";' 'i; }' | bc`
}

prepare_datadir()
{
  # Create list of subdirectories
  # If the file NOERASE exists, the old directories are not erased
  #   (start.x also knows then that var.dat is not created)
  for dir in $procdirs $subdirs
  do
    # Make sure a sufficient number of subdirectories exist
    ddir="$datadir/$dir"
    if [ ! -e $ddir ]; then
      mkdir $ddir
    else
      # Clean up
      # when used with lnowrite=T, for example, we don't want to remove var.dat:
      list=`/bin/ls $ddir/VAR* $ddir/TAVG* $ddir/*.dat $ddir/*.info $ddir/slice* 2>/dev/null`
      if [ ! -e NOERASE ]; then
        for rmfile in $list
        do
          if [ $rmfile != $ddir/var.dat ]; then rm -f $rmfile >& /dev/null; fi
        done
      fi
    fi
  done

  # Clean up previous runs
  if [ ! -e NOERASE ]; then
    if [ -e "$datadir/time_series.dat" ] && [ ! -s "$datadir/time_series.dat" ] ; then
        mv $datadir/time_series.dat $datadir/time_series.`timestr`
    fi
    rm -f $datadir/*.dat $datadir/*.nml $datadir/param*.pro $datadir/index*.pro \
          $datadir/averages/* >& /dev/null
  fi
}

#------------------------------------------------------------------------------
#-- Execution Control and Monitoring Routines                        
#------------------------------------------------------------------------------
background_remote_top()
{
  # Copy output from `top' on run host to a file we can read from login server
  if [ "$remote_top" == "yes" ]; then
    remote-top >& remote-top.log &
  fi
}

tidy_rundir()
{
  # Clean up control and data files
  # NB. Don't remove NEWDIR it may have been put there on purpose so as
  #     to catch a crash and run something else instead.
  rm -f STOP RELOAD RERUN fort.20
}

save_jobid()
{
  # Write $PBS_JOBID to file (important when run is migrated within the same job)
  if [ -n "$PBS_JOBID" ]; then
    echo $PBS_JOBID " $1  " `date` >> $datadir/jobid.dat
  elif [ -n "$JOBID" ]; then
    echo $JOBID " $1  " `date` >> $datadir/jobid.dat
  fi
}

pencil_code_start()
{
  # Run start.x timestamping and timing appropriately
  date
  echo "$mpirun $mpirunops $npops $start_x $x_ops"
  echo "$mpirun $mpirunops $npops $start_x $x_ops" > run_command.log
  time $mpirun $mpirunops $npops $start_x $x_ops
  set start_status=$?		# save for exit
  date
}

pencil_code_run()
{
  # Run run.x timestamping and timing appropriately
  date
  echo "$mpirun $mpirunops $npops $run_x $x_ops" 
  echo "$mpirun $mpirunops $npops $run_x $x_ops" > run_command.log
  time $mpirun $mpirunops $npops $run_x $x_ops
  set run_status=$?		# save for exit
  date
}


check_RERUN()
{
  # look for RERUN file 
  # With this method one can only reload a new executable.
  # One cannot change directory, nor are the var.dat files returned to server.
  # See the NEWDIR method below for more options.
  if [ -e RERUN ]; then 
    rm -f RERUN
    echo
    echo "=============================================================================="
    echo "Rerunning in the *same* directory; current run status: $run_status"
    echo "We are *still* in: " `pwd`
    echo "=============================================================================="
    echo
  else
    run_finished=yes
  fi  
}

check_NEWDIR()
{
  # look for NEWDIR file 
  # if NEWDIR contains a directory name, then continue run in that directory
  if [ -e NEWDIR ]; then 
    if [ -s NEWDIR ]; then
      remove_lock
      olddir=$cwd
      cd `< NEWDIR`
      rm $olddir/NEWDIR
      (echo "stopped run:"; date; echo "new run directory:"; echo $cwd; echo "")\
         >> $olddir/$datadir/directory_change.log
      (date; echo "original run script is in:"; echo $olddir; echo "")\
         >> $datadir/directory_change.log
      echo
      echo "=============================================================================="
      echo "Rerunning in new directory; current run status: $run_status"
      echo "We are now in: " `pwd`
      echo "=============================================================================="
      echo
    else
      rm -f NEWDIR
      echo
      echo "=============================================================================="
      echo "Rerunning in the *same* directory; current run status: $run_status"
        echo "We are *still* in: " `pwd`
    echo "=============================================================================="
      echo
      echo "Rerunning; current run status: $run_status"
    fi
  else
    job_finished=yes
  fi  
}


#------------------------------------------------------------------------------
#-- Data/Binary Distribution and Collection Routines                        
#------------------------------------------------------------------------------
distribute_data_to_nodes()
{
  #
  #  If necessary, distribute var.dat from the server to the various nodes
  #
  if [ "$local_disc" == "yes" ]; then
    if [ "$one_local_disc" == "yes" ]; then     # one common local disc 
      for node in $nodelist
      do
        for d in `cd $datadir; ls -d proc* allprocs`
        do
          $SCP $datadir/$d/var.dat ${node}:$SCRATCH_DIR/$d/
          $SCP $datadir/$d/timeavg.dat ${node}:$SCRATCH_DIR/$d/
        done
        $SCP $datadir/allprocs/dxyz.dat ${node}:$SCRATCH_DIR/allprocs
      done
    else # one local disc per MPI process (Horseshoe, etc);
      i=0
      for node in $nodelist
      do
        echo "i = $i"
        j=$nprocpernode
        while [ $j != 0 ]
        do
          $SCP $datadir/proc$i/var.dat ${node}:$SCRATCH_DIR/proc$i/
          $SCP $datadir/proc$i/timeavg.dat ${node}:$SCRATCH_DIR/proc$i/
          i=$(( $i + 1 ))
          j=$(( $j - 1 ))
        done
        $SCP $datadir/allprocs/dxyz.dat ${node}:$SCRATCH_DIR/allprocs/      
      done
    fi
  fi
}

distribute_binary()
{ 
  # If necessary copy executable to $SCRATCH_DIR of master node
  if [ "$local_binary" == "yes" ]; then
    echo "ls $1 $SCRATCH_DIR before copying:"
    ls -lt $1 $SCRATCH_DIR
    cp $1 $SCRATCH_DIR
    echo "ls $1 $SCRATCH_DIR after copying:"
    ls -lt $1 $SCRATCH_DIR
  fi
}


background_copy_snapshots()
{
  # On machines with local scratch directory, initialize automatic
  # background copying of snapshots back to the data directory.
  if [ "$local_disc" == "yes" ]; then
    echo "Use local scratch disk"
    copy-snapshots -v >& copy-snapshots.log &
  fi
}

final_copy_snapshots()
{
  # On machines with local scratch disc, copy var.dat back to the data
  # directory
  if [ "$local_disc" == "yes" ]; then
    echo -n "Copying final var.dat back from local scratch disk...   "
    copy-snapshots -v var.dat 2>&1 > copy-snapshots2.log
    copy-snapshots -v dxyz.dat 2>&1 >> copy-snapshots2.log
    copy-snapshots -v timeavg.dat 2>&1 >> copy-snapshots2.log
    copy-snapshots -v crash.dat 2>&1 >> copy-snapshots2.log
    echo "Done."
  fi
}

unbackground_copy_snapshots()
{
  # Kill all backgrounded copy-snapshots
  if [ "$local_disc" == "yes" ]; then
    echo -n "Killing backgrounded copy_snapshots...   "
    pids=`ps -U $USER -o pid,command | grep -E 'remote-top|copy-snapshots' | sed 's/^ *//' | cut -d ' ' -f 1`
    echo "Shutting down processes ${pids}:"
    for p in $pids     # need to do in a loop, and check for existence, since
    do                 # some systems (Hitachi) abort this script when trying
                       # to kill non-existent processes
      echo "  pid $p"
      # should not need existence check when scripting with sh 
      #[ `ps -p $p | fgrep -c $p` ] 
      kill -SIGKILL $p
    done
    echo "Done."
  fi
}

#------------------------------------------------------------------------------
#-- Run Lock Manipulation 
#------------------------------------------------------------------------------
check_not_locked()
{
  # Prevent code from running twice (and removing files by accident)
  if [ -e LOCK ]; then
    echo ""
    echo "getconf.csh: found LOCK file"
    echo "This may indicate that the code is currently running in this directory"
    echo "If this is a mistake (eg after a crash), remove the LOCK file by hand:"
    echo "rm LOCK"
    # exit                        # won't work in a sourced file
    kill -SIGKILL $$			# full-featured suicide
  fi
}

create_lock()
{
  # Prevent code from running twice (and removing files by accident)
  [ ! -e NEVERLOCK ] && touch LOCK
}

remove_lock()
{
  # Prevent code from running twice (and removing files by accident)
  [ -e LOCK ] && rm -f LOCK 
}


#------------------------------------------------------------------------------
#-- System / Run Enquiry Routines                        
#------------------------------------------------------------------------------
check_sgi_fix()
{
  local os=`uname -s`
  if match "$os" "IRIX" ; then
    touch SGIFIX
  else
    rm -f SGIFIX
  fi
}

determine_mpi()
{
  # Are we running the MPI version?
  if [ -e src/Makefile.local ]; then 
    egrep -c '^[ 	]*MPICOMM[ 	]*=[ 	]*mpicomm' src/Makefile.local >/dev/null
    if [ $? -eq 0 ]; then
      mpi=yes
      echo "Running under MPI"
    else
      mpi=no
    fi
  else
    mpi=no
    echo "ERROR: Cannot find src/Makefile.local - Unable to determine MPI usage - assume mpi=$mpi"
  fi
}

determine_ncpus()
{
  # Determine number of CPUS
  if [ -e src/cparam.local ]; then 
    ncpus=`perl -ne '$_ =~ /^\s*integer\b[^\\!]*ncpus\s*=\s*([0-9]*)/i && print $1' src/cparam.local`
    nprocpernode=1
    echo $ncpus CPUs
  else
    ncpus=1
    echo "ERROR: Cannot find src/cparam.local - Unable to determine MPI usage - assume ncpus=1"
  fi
}

determine_nodelist()
{
  # Get list of nodes; filters lines such that it would also handle
  # machines.XXX files or lam-bhost.der, although this is hardly necessary.
  if [ -e "$PBS_NODEFILE" ]; then
    [ "$debug" == "yes" ] && echo "PBS job"
    nodelist=`cat $PBS_NODEFILE | grep -v '^#' | sed 's/:[0-9]*//'`
  elif [ -e "$PE_HOSTFILE" ]; then
    [ "$debug" == "yes" ] && echo "SGE Parallel Environment job - $PE"
    nodelist = `cat $PE_HOSTFILE | grep -v '^#' | sed 's/\ .*//'`
  elif [ -n "$JOB_ID" ]; then
    if [ -e "$HOME/.score/ndfile.$JOB_ID" ]; then
      [ "$debug" == "yes" ] && echo "Scout job"
      nodelist=`cat $HOME/.score/ndfile.$JOB_ID | grep -v '^#' | sed 's/:[0-9]*//'`
    else
      # E.g. for wd running SGE on Kabul cluster
      echo "Apparently not a scout job"
    fi
  else
    [ "$debug" == "yes" ] && echo "Setting nodelist to ($hn)"
    nodelist="$hn"
  fi
}

determine_mpi_processor_options()
{
  ## MPI specific setup
  if [ "$mpi" == "yes" ]; then
    # Some mpiruns need special options
    if match "$mpirun" "mpirun" ; then
      npops="-np $ncpus"
    elif match "$mpirun" "mpiexec" ; then
      npops="-n $ncpus"
    elif match "$mpirun" "scout" ; then
      nnode=`expr $NSLOTS - 1`
      nprocpernode=`expr $ncpus / $nnode` 
      npops="-nodes=${nnode}x${nprocpernode}"
    elif match "$mpirun" "poe" ; then
      nprocpernode=1
      npops="-procs $ncpus"
    else
      echo "getconf.csh: No clue how to tell $mpirun to use $ncpus nodes"
    fi

  else # no MPI
    echo "Non-MPI version"
    mpirun=''
    mpirunops=''
    npops=''
    ncpus=1
    nprocpernode=1
  fi
}

check_reference_data()
{
  if [ -e reference.out ] && [ -e data/time_series.dat ]; then
    diff reference.out data/time_series.dat > /dev/null 
    if [ $? -ne 0 ]; then 
      echo "Found reference data: reference.out"
      echo "WARNING: Time series output does not match reference data."
    else
      echo "Found reference data - reference.out, which matches time_series.dat"
    fi
  fi
}

check_datadir()
{
 #
 #  If we don't have a data subdirectory: stop here (it is too easy to
 #  continue with an NFS directory until you fill everything up).
 #
 if [ ! -d "$datadir" ]; then
   echo ""
   echo ">>  STOPPING: need $datadir directory"
   echo ">>  Recommended: create $datadir as link to directory on a fast scratch"
   echo ">>  Not recommended: you can generate $datadir with 'mkdir $datadir', "
   echo ">>  but that will most likely end up on your NFS file system and be"
   echo ">>  slow"
   echo
   exit 0
 fi
}


ishost() {
  if perl -e "if ( \"$hn\" =~ /$1/ ) { exit(0) } else { exit (1) }"
  then
    return 0;
  else
    return 1;
  fi
}


#------------------------------------------------------------------------------
#-- Variable Defaults                        
#------------------------------------------------------------------------------
TIMEFORMAT="real %3lR, system %3S, user %3U, cpu usage %P%%"

