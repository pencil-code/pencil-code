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
#  kill $$
  exit 1
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
  echo '$remote_top     = ' "<$remote_top>"
  echo '$local_binary   = ' "<$local_binary>"
  echo '$datadir        = ' "<$datadir>"
  echo '$SCRATCH_DIR    = ' "<$SCRATCH_DIR>"
  echo '$hn             = ' "<$hn>"
  echo '$start_x        = ' "<$start_x>"
  echo '$run_x          = ' "<$run_x>"
  echo '$mpirun         = ' "<$mpirun>"
  echo '$mpirunops      = ' "<$mpirunops>"
  echo '$x_ops          = ' "<$x_ops>"
  echo '$NODELIST       = ' "<$NODELIST>"
  echo '$SSH            = ' "<$SSH>"
  echo '$SCP            = ' "<$SCP>"
  echo '$PARENT_PID     = ' "<$PARENT_PID>"
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
  [ "$local_disc" = "yes" ] &&  echo $SCRATCH_DIR >$datadir/directory_snap


  if [ "$local_binary" = "yes" ]; then
    start_x=$SCRATCH_DIR/start.x
    run_x=$SCRATCH_DIR/run.x
  fi

  # Created subdirectories on local scratch disc (start.csh will also create
  # them under $datadir/)
  if [ "$local_disc" = "yes" ]; then
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
  if [ "$local_disc" = "yes" ]; then
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
  subdirs="allprocs averages idl"
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
    if [ ! -d $ddir ]; then
      mkdir $ddir
    else
      # Clean up
      # when used with lnowrite=T, for example, we don't want to remove var.dat:
      list=`/bin/ls $ddir/VAR* $ddir/TAVG* $ddir/*.dat $ddir/*.info $ddir/slice* 2>/dev/null`
      if [ ! -f NOERASE ]; then
        for rmfile in $list
        do
          if [ $rmfile != $ddir/var.dat ]; then rm -f $rmfile &>/dev/null; fi
        done
      fi
    fi
  done

  # Clean up previous runs
  if [ ! -f NOERASE ]; then
    if [ -f "$datadir/time_series.dat" ] && [ ! -s "$datadir/time_series.dat" ] ; then
        mv $datadir/time_series.dat $datadir/time_series.`timestr`
    fi
    rm -f $datadir/*.dat $datadir/*.nml $datadir/param*.pro $datadir/index*.pro \
          $datadir/averages/* &>/dev/null
  fi
}

#------------------------------------------------------------------------------
#-- Execution Control and Monitoring Routines
#------------------------------------------------------------------------------
background_remote_top()
{
  # Copy output from `top' on run host to a file we can read from login server
  if [ "$remote_top" = "yes" ]; then
    remote-top &>remote-top.log &
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
  date
  echo "#" >> pc_commands.log
  date +'# %Y-%m-%d %H:%M:%S' >> pc_commands.log
  echo "# $mpirun $mpirunops $npops $start_x $x_ops" >> pc_commands.log
  time $mpirun $mpirunops $npops $start_x $x_ops
  start_status=$?               # save for exit
  date
}

pencil_code_run()
{
  # Run run.x timestamping and timing appropriately
  date
  echo "$mpirun $mpirunops $npops $run_x $x_ops"
  date
  echo "#" >> pc_commands.log
  date +'# %Y-%m-%d %H:%M:%S' >> pc_commands.log
  echo "# $mpirun $mpirunops $npops $run_x $x_ops" >> pc_commands.log
  time $mpirun $mpirunops $npops $run_x $x_ops
  run_status=$?         # save for exit
  date
}

queue_submit()
{
  # Submit a pencil-code run to the queue system
  # This function should be over ridden in pc_config.sh

  echo "------------------------------------------------------------------"
  echo " ERROR: NO QUEUE SUBMISSION DEFINED           "
  echo "                                              "
  echo " You have not defined how to submit to the queue on this machine. "
  echo " To do so add a definition of the queue_submit shell function     "
  echo " to the appropriate machine dependent section of pc_config.sh.    "
  echo "                                              "
  echo " eg.                                          "
  echo "  queue_submit()                              "
  echo "  {                                           "
  echo "    qsub -pe mpi \$ncpus \$1                  "
  echo "  }                                           "
  echo "                                              "
  echo "  You may use and variables like \$ncpus defined by running the    "
  echo "  pc_config.sh script. \$1 will contain the name of the executable "
  echo "  script suggested for submission. \$2 will contain the suggested  "
  echo "  for the job. \$3, ..., etc. will contain any other parameters    "
  echo "  passed on the pc_qsub command line.                              "
  echo "    NB. one may use 'shift 2' to get rid of \$1, \$2 then \$@ to   "
  echo "        refer to the rest of the parameters.                       "
  echo "-------------------------------------------------------------------"
}

check_RERUN()
{
  # look for RERUN file
  # With this method one can only reload a new executable.
  # One cannot change directory, nor are the var.dat files returned to server.
  # See the NEWDIR method below for more options.
  if [ -f RERUN ]; then
    rm -f RERUN
    echo
    echo "------------------------------------------------------------------------------"
    echo "Rerunning in the *same* directory; current run status: $run_status"
    echo "We are *still* in: " `pwd`
    echo "------------------------------------------------------------------------------"
    echo
  else
    run_finished=yes
  fi
}

check_NEWDIR()
{
  # look for NEWDIR file
  # if NEWDIR contains a directory name, then continue run in that directory
  if [ -f NEWDIR ]; then
    if [ -s NEWDIR ]; then
      remove_lock
      olddir=$cwd
      cd `< NEWDIR`
      rm $olddir/NEWDIR
      (echo "stopped run:"; date; echo "new run directory:"; echo $cwd; echo "")\
         2>&1 >> $olddir/$datadir/directory_change.log
      (date; echo "original run script is in:"; echo $olddir; echo "")\
         2>&1 >> $datadir/directory_change.log
      echo
      echo "------------------------------------------------------------------------------"
      echo "Rerunning in new directory; current run status: $run_status"
      echo "We are now in: " `pwd`
      echo "------------------------------------------------------------------------------"
      echo
    else
      rm -f NEWDIR
      echo
      echo "------------------------------------------------------------------------------"
      echo "Rerunning in the *same* directory; current run status: $run_status"
        echo "We are *still* in: " `pwd`
    echo "------------------------------------------------------------------------------"
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
  if [ "$local_disc" = "yes" ]; then
    if [ "$one_local_disc" = "yes" ]; then     # one common local disc
      echo -n "Distributing data to nodes...  "
      for node in $nodelist
      do
        for d in `cd $datadir; ls -d proc* allprocs`
        do
          $SCP ${remote_data_path}${datadir}/$d/var.dat ${node}:$SCRATCH_DIR/$d/ 2>/dev/null
          $SCP ${remote_data_path}$datadir/$d/timeavg.dat ${node}:$SCRATCH_DIR/$d/ 2>/dev/null
        done
        $SCP ${remote_data_path}$datadir/allprocs/dxyz.dat ${node}:$SCRATCH_DIR/allprocs 2>/dev/null
      done
      echo " Done."
    else # one local disc per MPI process (Horseshoe, etc);
      echo -n "Distributing data to nodes...  "
      i=0
      for node in $nodelist
      do
        j=$nprocpernode
        while [ $j -gt 0 ]
        do
          echo -n "$i "
          $SCP ${remote_data_path}$datadir/proc$i/var.dat ${node}:$SCRATCH_DIR/proc$i/ 2>/dev/null
          $SCP ${remote_data_path}$datadir/proc$i/timeavg.dat ${node}:$SCRATCH_DIR/proc$i/ 2>/dev/null
          i=`expr $i + 1`
          j=`expr $j - 1`
        done
        $SCP $datadir/allprocs/dxyz.dat ${node}:$SCRATCH_DIR/allprocs/ 2>/dev/null
      done
      echo " Done."
    fi
  fi
}

distribute_binary()
{
  # If necessary copy executable to $SCRATCH_DIR of master node
  if [ "$local_binary" = "yes" ]; then
    echo "ls $1 $SCRATCH_DIR before copying:"
    ls -lt $1 $SCRATCH_DIR
    cp $1 $SCRATCH_DIR
    echo "ls $1 $SCRATCH_DIR after copying:"
    ls -lt $1 $SCRATCH_DIR
  fi
}


copy_snapshots ()
{
  # Copy snapshots VAR# and TAVG# from /scratch to $PBS_O_WORKDIR (from
  # where the code started).
  # Relies on PBS (and on tcsh) and is needed on Horseshoe (the Odense
  # cluster).
  debug=no

  if [ "$1" = "-v" ] || [ "$1" = "--verbose" ]; then
    debug=yes
    shift
  fi

  #[ "$debug" = "yes" ] && (set verbose;  set echo)

  varfile=$1
  [ "$debug" = "yes" ] && echo "varfile = <$varfile>"

  pwd=`pwd`
  targetdir=$pwd/data
  nodelist=`echo $NODELIST | sed 's/:/ /g'` # unpack NODELIST

  if [ "$debug" = "yes" ]; then
    echo "SCRATCH_DIR = <$SCRATCH_DIR>"
    echo "targetdir   = <$targetdir>"
    echo "nodelist    = <$nodelist>"
  fi

  if [ -n "$varfile" ]; then            # explicit filename given
    for node in $nodelist
    do
      echo "---------------------- $node ---------------------------"
      if [ "$debug" = "yes" ]; then
        echo "node=$node"
        printf "\n$SSH $node ls -ltd $SCRATCH_DIR $SCRATCH_DIR/proc*/$varfile $SCRATCH_DIR $SCRATCH_DIR/allprocs/$varfile :\n"
        $SSH $node "ls -ltd $SCRATCH_DIR $SCRATCH_DIR/proc*/$varfile $SCRATCH_DIR $SCRATCH_DIR/allprocs/$varfile"
        printf "\n$SSH $node ls -ltd $targetdir/proc*/$varfile $targetdir/allprocs/ :\n"
        $SSH $node "ls -ltd $targetdir/proc*/ $targetdir/allprocs/"
        echo
      fi
      # Copy all files you find (not efficient on dual-CPU systems, unless
      # we would delete files after copying..)
      # Extremely awkward due to quoting rules:
      cmd1="cd $SCRATCH_DIR; "
      cmd2='for f in `'
      cmd3="ls proc*/$varfile allprocs/$varfile"
      cmd4='`; do cp '
      cmd5="$SCRATCH_DIR/"
      cmd6='$f '
      cmd7="$targetdir/"
      cmd8='$f; done'
      remcmd="$cmd1$cmd2$cmd3$cmd4$cmd5$cmd6$cmd7$cmd8"
      [ $debug = "yes" ]  && echo "Now running <$SSH $node sh -c $remcmd>"
      $SSH $node sh -c "'$remcmd'"
      printf "\n$SSH $node ls -ltd $targetdir/proc*/$varfile $targetdir/allprocs/$varfile :\n"
      $SSH $node "ls -ltd $targetdir/proc*/$varfile $targetdir/allprocs/$varfile"
      echo "--------------------------------------------------------"
    done
  else                          # no explicit file given -- copy VARN, TAVGN
    if [ -f COPY_IN_PROGRESS ]; then
      echo "ERROR: Copy already in progress! Exiting..."
      return 1
    fi

    touch COPY_IN_PROGRESS
    while [ -f COPY_IN_PROGRESS ];   # loop until killed
    do
      sleep 60 &                        # only check every minute
      save_sleep_pid=$!
      trap "kill $save_sleep_pid; rm -f COPY_IN_PROGRESS" SIGQUIT
      wait $save_sleep_pid
      [ -f COPY_IN_PROGRESS ] && date >> COPY_IN_PROGRESS

      ## Is there something to copy? (It is sufficient to copy once the
      ## files show up for the master process of this job -- might depend on
      ## the queuing system etc.)
      if ( (ls $SCRATCH_DIR/proc0/VAR[0-9]*     || \
              ls $SCRATCH_DIR/proc0/TIMEAVG[0-9]* || \
              ls $SCRATCH_DIR/allprocs/VAR[0-9]*  || \
              ls $SCRATCH_DIR/allprocs/TIMEAVG[0-9]* ) >& /dev/null ) ; then
        ## Decide whether to check for file size (old scheme; won't work with
        ## TAVGN unless they have the same size as var.dat), or to use list of
        ## snapshot files (new scheme, will not work with old binaries).
        if [ -f $SCRATCH_DIR/proc0/varN.list ]     || \
           [ -f $SCRATCH_DIR/proc0/tavgN.list ]    || \
           [ -f $SCRATCH_DIR/allprocs/tavgN.list ] || \
           [ -f $SCRATCH_DIR/allprocs/tavgN.list ];  then
          ## New scheme
          echo "New scheme for copying (based on varN.list)"
          for node in $nodelist
          do
            echo "---------------------- $node ---------------------------"
            echo '$targetdir' $targetdir
            for d in `cd $SCRATCH_DIR; ls -d allprocs proc*`
            do
              fdir="$SCRATCH_DIR/$d"
              echo "cd $SCRATCH_DIR; ls allprocs proc* -->"
              cd $SCRATCH_DIR; ls allprocs proc*
              echo "d = <$d>"
              echo "fdir = <$fdir>"
              for f in `$SSH $node "cat $fdir/*.list"`
              do
                file="$fdir/$f"
                echo "f     = <$f>"
                echo "file = <$file>"
                [ "$debug" = "yes" ] && echo "$SSH $node 'if (-e $file) mv -f $file $targetdir/$d/'"
                $SSH $node "if (-e $file) mv -f $file $targetdir/$d/"
              done
            done
            echo "--------------------------------------------------------"
          done
        else
          ## Old scheme
          echo "Old scheme for copying (size-based)-- will be removed at some point"
        fi

        if [ `ps -p $PARENT_PID | fgrep -c $PARENT_PID` -le 0 ]; then
          echo "No parent process (pid $PARENT_PID) -- exiting"
          exit 1
        fi
      fi
    done
  fi

  [ "$debug" = "yes" ] && echo "copy_snapshots: done"
}

#copy_snapshots ()
#{
#  # Copy snapshots VAR# and TAVG# from /scratch to $PBS_O_WORKDIR (from
#  # where the code started).
#  # Relies on PBS (and on tcsh) and is needed on Horseshoe (the Odense
#  # cluster).
#  local debug=no
#
#  if [ "$1" = "-v" ] || [ "$1" = "--verbose" ]; then
#    debug=yes
#    shift
#  fi
#
#  #[ "$debug" = "yes" ] && (set verbose;  set echo)
#
#  local varfile=$1
#  [ "$debug" = "yes" ] && echo "varfile = <$varfile>"
#
#  local pwd=`pwd`
#  local targetdir=$pwd/data
##  local nodelist=`echo $NODELIST | sed 's/:/ /g'` # unpack NODELIST
#
#  if [ "$debug" = "yes" ]; then
#    echo "SCRATCH_DIR = <$SCRATCH_DIR>"
#    echo "targetdir   = <$targetdir>"
#    echo "nodelist    = <$nodelist>"
#  fi
#
#  if [ -n "$varfile" ]; then           # explicit filename given
#    for node in $nodelist
#    do
#      echo "---------------------- $node ---------------------------"
#      if [ "$debug" = "yes" ]; then
#        echo "node=$node"
#        printf "\n$SSH $node ls -ltd $SCRATCH_DIR $SCRATCH_DIR/proc*/$varfile $SCRATCH_DIR $SCRATCH_DIR/allprocs/$varfile :\n"
#        $SSH $node "ls -ltd $SCRATCH_DIR $SCRATCH_DIR/proc*/$varfile $SCRATCH_DIR $SCRATCH_DIR/allprocs/$varfile" 2>/dev/null
#        printf "\n$SSH $node ls -ltd $targetdir/proc*/$varfile $targetdir/allprocs/ :\n"
#        $SSH $node "ls -ltd $targetdir/proc*/ $targetdir/allprocs/" 2>/dev/null
#        echo
#      fi
#      # Copy all files you find (not efficient on dual-CPU systems, unless
#      # we would delete files after copying..)
#      # Extremely awkward due to quoting rules:
#      cmd1="cd $SCRATCH_DIR; "
#      cmd2='for f in `'
#      cmd3="ls proc*/$varfile allprocs/$varfile 2>/dev/null"
#      cmd4='`; do cp '
#      cmd5="$SCRATCH_DIR/"
#      cmd6='$f '
#      cmd7="$targetdir/"
#      cmd8='$f; done'
#      remcmd="$cmd1$cmd2$cmd3$cmd4$cmd5$cmd6$cmd7$cmd8"
#      [ $debug = "yes" ]  && echo "Now running <$SSH $node sh -c $remcmd>"
#      $SSH $node sh -c "'$remcmd'"
#      printf "\n$SSH $node ls -ltd $targetdir/proc*/$varfile $targetdir/allprocs/$varfile :\n"
#      $SSH $node "ls -ltd $targetdir/proc*/$varfile $targetdir/allprocs/$varfile" 2>/dev/null
#      echo "--------------------------------------------------------"
#    done
#  else                         # no explicit file given -- copy VARN, TAVGN
#    if [ -f COPY_IN_PROGRESS ]; then
#      echo "ERROR: Copy already in progress! Exiting..."
#      return 1
#    fi
#
#    touch COPY_IN_PROGRESS
#    while [ -f COPY_IN_PROGRESS ];   # loop until killed
#    do
#      sleep 60 &                    # only check every minute
#      local save_sleep_pid=$!
#      trap "kill $save_sleep_pid; rm -f COPY_IN_PROGRESS" SIGQUIT
#      wait $save_sleep_pid
#      [ -f COPY_IN_PROGRESS ] && date >> COPY_IN_PROGRESS
#
#      ## Is there something to copy? (It is sufficient to copy once the
#      ## files show up for the master process of this job -- might depend on
#      ## the queuing system etc.)
#      if ( (ls $SCRATCH_DIR/proc0/VAR[0-9]*     || \
#              ls $SCRATCH_DIR/proc0/TIMEAVG[0-9]* || \
#              ls $SCRATCH_DIR/allprocs/VAR[0-9]*  || \
#              ls $SCRATCH_DIR/allprocs/TIMEAVG[0-9]* ) &> /dev/null ) ; then
#        ## Decide whether to check for file size (old scheme; won't work with
#        ## TAVGN unless they have the same size as var.dat), or to use list of
#        ## snapshot files (new scheme, will not work with old binaries).
#        if [ -f $SCRATCH_DIR/proc0/varN.list ]     || \
#           [ -f $SCRATCH_DIR/proc0/tavgN.list ]    || \
#           [ -f $SCRATCH_DIR/allprocs/tavgN.list ] || \
#           [ -f $SCRATCH_DIR/allprocs/tavgN.list ];  then
#          ## New scheme
#          echo "New scheme for copying (based on varN.list)"
#          for node in $nodelist
#          do
#            echo "---------------------- $node ---------------------------"
#            echo '$targetdir' $targetdir
#            for d in `cd $SCRATCH_DIR; ls -d allprocs proc* 2>/dev/null`
#            do
#              fdir="$SCRATCH_DIR/$d"
#              echo "cd $SCRATCH_DIR; ls allprocs proc* -->"
#              cd $SCRATCH_DIR; ls allprocs proc*
#              echo "d = <$d>"
#              echo "fdir = <$fdir>"
#              for f in `$SSH $node "cat $fdir/*.list"`
#              do
#                file="$fdir/$f"
#                echo "f     = <$f>"
#                echo "file = <$file>"
#                [ "$debug" = "yes" ] && echo "$SSH $node 'if (-e $file) mv -f $file $targetdir/$d/'"
#                $SSH $node "if (-e $file) mv -f $file $targetdir/$d/"
#              done
#            done
#            echo "--------------------------------------------------------"
#          done
#        else
#          ## Old scheme
#          echo "Old scheme for copying (size-based)-- will be removed at some point"
#        fi
#
#        if [ `ps -p $PARENT_PID | fgrep -c $PARENT_PID` -le 0 ]; then
#          rm -f COPY_IN_PROGRESS
#          echo "No parent process (pid $PARENT_PID) -- exiting"
#          exit 1
#        fi
#      fi
#    done
#  fi
#
#  rm -f COPY_IN_PROGRESS
#
#  [ "$debug" = "yes" ] && echo "copy_snapshots: done"
#}

background_copy_snapshots()
{
  # On machines with local scratch directory, initialize automatic
  # background copying of snapshots back to the data directory.
  if [ "$local_disc" = "yes" ]; then
    echo -n "Starting background copy snapshots...  "
    copy_snapshots -v &> copy-snapshots.log &
    bkgnd_copy_snapshots_pid=$!
    echo "Started."
  fi
}

final_copy_snapshots()
{
  # On machines with local scratch disc, copy var.dat back to the data
  # directory
  if [ "$local_disc" = "yes" ]; then
    echo -n "Copying final var.dat back from local scratch disk...   "
    copy_snapshots -v var.dat &> copy-snapshots2.log
    copy_snapshots -v dxyz.dat 2>&1 >> copy-snapshots2.log
    copy_snapshots -v timeavg.dat 2>&1 >> copy-snapshots2.log
    copy_snapshots -v crash.dat 2>&1 >> copy-snapshots2.log
    echo "Done."
  fi
}

unbackground_copy_snapshots()
{
  # Kill all backgrounded copy-snapshots
  if [ "$local_disc" = "yes" ]; then
    echo -n "Stopping backgrounded copy_snapshots...   Waiting...  "
#    kill $bkgnd_copy_snapshots_pid
    kill -SIGQUIT $bkgnd_copy_snapshots_pid 2>&1 >/dev/null
    wait $bkgnd_copy_snapshots_pid
    echo "Done."
  fi
}

#------------------------------------------------------------------------------
#-- Run Lock Manipulation
#------------------------------------------------------------------------------
check_is_run_directory()
{
  isrundir=yes

  if [ ! -f start.in ]; then
    echo "ERROR: Cannot find start.in"
    isrundir=no
  fi

  if [ ! -f run.in ]; then
    echo "ERROR: Cannot find run.in"
    isrundir=no
  fi

  if [ ! "$isrundir" = "yes" ]; then
    echo "This does not appear to be a run directory! Exiting... "
    exit 1
#    kill -SIGKILL $$                   # full-featured suicide
  fi
}

check_not_locked()
{
  # Prevent code from running twice (and removing files by accident)
  if [ -f LOCK ]; then
    echo ""
    echo "getconf.csh: found LOCK file"
    echo "This may indicate that the code is currently running in this directory"
    echo "If this is a mistake (eg after a crash), remove the LOCK file by hand:"
    echo "rm LOCK"
    exit 1                        # won't work in a sourced file
#    kill -SIGKILL $$                   # full-featured suicide
  fi
}

create_lock()
{
  # Prevent code from running twice (and removing files by accident)
  [ ! -f NEVERLOCK ] && touch LOCK
}

remove_lock()
{
  # Prevent code from running twice (and removing files by accident)
  [ -f LOCK ] && rm -f LOCK
}


#------------------------------------------------------------------------------
#-- System / Run Enquiry Routines
#------------------------------------------------------------------------------
check_sgi_fix()
{
  os=`uname -s`
  if match "$os" "IRIX" ; then
    touch SGIFIX
  else
    rm -f SGIFIX
  fi
}

determine_mpi()
{
  # Are we running the MPI version?
  if [ -f src/Makefile.local ]; then
    egrep -c '^[        ]*MPICOMM[      ]*=[    ]*mpicomm' src/Makefile.local >/dev/null
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
  if [ -f src/cparam.local ]; then
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
echo -n "Determine nodelist... "
  if [ -f "$PBS_NODEFILE" ]; then
    echo "PBS job"
    nodelist=`cat $PBS_NODEFILE | grep -v '^#' | sed 's/:[0-9]*//'`
  elif [ -f "$PE_HOSTFILE" ]; then
    echo "SGE Parallel Environment job - $PE"
    nodelist=`cat $PE_HOSTFILE | grep -v '^#' | sed 's/\ .*//'`
  elif [ -n "$JOB_ID" ]; then
    if [ -f "$HOME/.score/ndfile.$JOB_ID" ]; then
      echo "Scout job"
      nodelist=`cat $HOME/.score/ndfile.$JOB_ID | grep -v '^#' | sed 's/:[0-9]*//'`
    else
      # E.g. for wd running SGE on Kabul cluster
      echo "Apparently not a scout job"
    fi
  else
    echo "Setting nodelist to ($hn)"
    nodelist="$hn"
  fi
}

determine_mpi_processor_options()
{
  ## MPI specific setup
  if [ "$mpi" = "yes" ]; then
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
  if [ -f reference.out ] && [ -f data/time_series.dat ]; then
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

