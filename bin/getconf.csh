#!/bin/csh

# Name:   getconf.csh
# Author: wd (Wolfgang.Dobler@ncl.ac.uk)
# Date:   16-Dec-2001
# $Id: getconf.csh,v 1.120 2004-08-10 14:34:58 nilshau Exp $
#
# Description:
#  Initiate some variables related to MPI and the calling sequence, and do
#  other setup needed by both start.csh and run.csh.

set debug = 1
# set verbose
# set echo

# Just as a keepsake 
set dollar = '$'

# Set up PATH for people who don't include $PENCIL_HOME/bin by default
if ($?PENCIL_HOME) setenv PATH ${PATH}:${PENCIL_HOME}/bin

# Save working directory for other scripts we call
setenv PENCIL_WORKDIR `pwd`

newdir:
# Prevent code from running twice (and removing files by accident)
if (-e "LOCK") then
  echo ""
  echo "getconf.csh: found LOCK file"
  echo "This may indicate that the code is currently running in this directory"
  echo "If this is a mistake (eg after a crash), remove the LOCK file by hand: rm LOCK"
  echo ""
  echo "THERE IS STILL ONE CHANCE; if a NEWDIR file exist we change directory:"
  if (-e "NEWDIR") then 
    if (-s "NEWDIR") then
      set olddir=$cwd
      cd `cat NEWDIR`
      # Up here, datadir is not yet defined, so do it locally:
      # Determine data directory (defaults to `data')
      if (-r datadir.in) then
        set datadir = `cat datadir.in | sed 's/ *\([^ ]*\).*/\1/'`
      else
        set datadir = "data"
      endif
      # Write some info into log files of current and new directories
      (echo "redirected run:"; date; echo "redirected run directory:"; echo $cwd; echo "")\
         >> $olddir/$datadir/directory_change.log
      (date; echo "redirected FROM original run script in:"; echo $olddir; echo "")\
         >> $datadir/directory_change.log
      echo "WELL DONE: redirected to new directory:"
      echo `pwd`
      echo "... wait 20 sec, to allow manual escape if we have produced a loop!"
      echo
      # In new run directory, need to check for yet another possible LOCK file
      # Sleep for 20 sec, to allow manual escape if we have produced a loop!
      sleep 20
      goto newdir
    endif
  else
    echo "BAD LUCK: there is no NEWDIR file"
    echo
    # exit                        # won't work in a sourced file
    kill -KILL $$			# full-featured suicide
  endif  
endif

# Are we running the MPI version?
set mpi = `egrep -c '^[ 	]*MPICOMM[ 	]*=[ 	]*mpicomm' src/Makefile.local`
# Determine number of CPUS
set ncpus = `perl -ne '$_ =~ /^\s*integer\b[^\\!]*ncpus\s*=\s*([0-9]*)/i && print $1' src/cparam.local`
echo "$ncpus CPUs"

# Location of executables and other default settings; can be overwritten
# below
set start_x = "src/start.x"
set run_x   = "src/run.x"
set x_ops = ""         # arguments to both start.x and run.x
set mpirun = 'mpirun'

# Settings for machines with local data disks
# local_disc     = 1 means processes uses their own local scratch disc(s) for
#                    large files like var.dat, VARN, etc.
# one_local_disc = 1 means one common scratch disc for all processes
# local_binary   = 1 means copy executable to local scratch area of master node
# remote_top     = 1 means get `top' output in regular intervals
set local_disc     = 0
set one_local_disc = 1		# probably the more common case
set remote_top     = 0
set local_binary   = 0
set remove_scratch_root = 0
setenv SCRATCH_DIR /scratch
setenv SSH ssh
setenv SCP scp

# Any lam daemons to shutdown at the end of run.csh?
set booted_lam = 0

echo `uname -a`
set hn = `uname -n`
if ($mpi) echo "Running under MPI"
set mpirunops = ''

# Get list of nodes; filters lines such that it would also handle
# machines.XXX files or lam-bhost.der, although this is hardly necessary.
if ($?PBS_NODEFILE) then
  if ($debug) echo "PBS job"
  set nodelist = `cat $PBS_NODEFILE | grep -v '^#' | sed 's/:[0-9]*//'`
else if ($?PE_HOSTFILE) then
  if ($debug) echo "SGE Parallel Environment job - $PE"
  set nodelist = `cat $PE_HOSTFILE | grep -v '^#' | sed 's/\ .*//'`
else if ($?JOB_ID) then
  if (-e $HOME/.score/ndfile.$JOB_ID) then
    if ($debug) echo "Scout job"
    set nodelist = `cat $HOME/.score/ndfile.$JOB_ID | grep -v '^#' | sed 's/:[0-9]*//'`
  else
    # E.g. for wd running SGE on Kabul cluster
    echo "Apparently not a scout job"
  endif
else
  if ($debug) echo "Setting nodelist to ($hn)"
  set nodelist = ("$hn")
endif
# Output information about number of cpus per node
set nnodes = $#nodelist
set nprocpernode = `expr $ncpus / $nnodes`
echo "$nnodes nodes, $nprocpernode CPU(s) per node"

## ------------------------------
## Choose machine specific settings

if ($hn =~ mhd*.st-and.ac.uk) then
  echo "St Andrews machine"
  set mpirun = "dmpirun"

else if ($hn =~ *.kis.uni-freiburg.de) then
  set mpirun = /opt/local/mpich/bin/mpirun

else if ($hn =~ giga[0-9][0-9].ncl.ac.uk) then
  echo "Newcastle e-Science Cluster"
  if ($?PE) then
    echo "SGE job"
    set local_disc = 1
    set one_local_disc = 0
    set local_binary = 0

  #  cat $PE_HOSTFILE | sed 's/\([[:alnum:].-]*\)\ \([0-9]*\).*/for ( i=0 \; i < 2 \; i++ ){print "\1\\n"};/' | bc > hostfile
  #  set nodelist = `cat hostfile`

    if ($PE =~ mpi2) then
      set nprocpernode = 2
    else if ($PE =~ mpi) then
      set nprocpernode = 4
    endif
  else
    echo "Non-SGE, running on `hostname`"
  endif
  set mpirun = /addon/shared/lam/bin/mpirun
  set mpirunops = "-O -ssi rpi tcp -s n0 -x LD_ASSUME_KERNEL=2.4.1" #Fix bug in Redhat 9
  if ($local_disc) then
    #setenv SCRATCH_DIR `cat $TMPDIR/scratch` 
    setenv SCRATCH_DIR /work/$JOB_ID 
    set remove_scratch_root = 1
  endif

else if ($hn =~ giga[0-9][0-9]) then
  echo "Nordita cluster"
  if ($?PBS_NODEFILE) then
    echo "PBS job"
    cat $PBS_NODEFILE > lamhosts
    set local_disc = 1
    set one_local_disc = 0
    set local_binary = 0
  else
    echo "Non-PBS, running on `hostname`"
    echo `hostname` > lamhosts
  endif
  # lamboot with logging to file 
  echo "lambooting.."
  echo 'pidof lamd | xargs ps -lfwww :' > lamboot.log
  pidof lamd | xargs ps -lfwww >> lamboot.log
  echo '------------------------------------------------' >> lamboot.log
  lamboot -v lamhosts >>& lamboot.log # discard output, or auto-test
                                      # mysteriously hangs on Nq0
  echo '------------------------------------------------' >> lamboot.log
  echo 'pidof lamd | xargs ps -lfwww :' >> lamboot.log
  pidof lamd | xargs ps -lfwww >> lamboot.log
  #
  set booted_lam = 1
  echo "lamnodes:"
  lamnodes
  set mpirun = /opt/lam/bin/mpirun
  set mpirun = /usr/bin/mpirun
  set mpirunops = "-O -c2c -s n0 -x LD_ASSUME_KERNEL=2.4.1"#Fix bug in Redhat 9
  if ($local_disc) then
    setenv SCRATCH_DIR "/var/tmp"
  endif

else if (($hn =~ sleipner) || ($hn =~ fenris) || ($hn =~ hugin) || ($hn =~ munin)) then
  set mpirun = /usr/bin/mpiexec
  set local_disc = 1
  set one_local_disc = 1
  set local_binary = 0
  if ($?SCRDIR) then
    if (-d $SCRDIR) then
      setenv SCRATCH_DIR "$SCRDIR"
    else
      echo 'NO SUCH DIRECTORY: $SCRDIR'="<$SCRDIR> -- ABORTING"
      kill $$                     # full-featured suicide
    endif
  else
    echo 'NO SUCH VARIABLE: $SCRDIR'="<$SCRDIR> -- ABORTING"
    kill $$                     # full-featured suicide
  endif
  setenv SSH rsh
  setenv SCP rcp
  setenv LANG en_US


else if ($hn =~ node*.clusters.com) then
  echo "fire in Bergen"
  set mpirunops = ''
  set mpirun = 'mpirunpbs'
  set npops = "-np $ncpus"
#  set mpi = 1


else if ( ($hn =~ cincinnatus*) || ($hn =~ owen*) \
          || ($hn =~ master) || ($hn =~ node*) \
	  || ($hn =~ ns0*) ) then
  if ($mpi) then
    # Choose appropriate mpirun version (LAM vs. MPICH)
    if (`fgrep -c lam_mpi src/start.x` > 0) then # lam
      if (-x /usr/lib/lam/bin/mpirun) set mpirun=/usr/lib/lam/bin/mpirun
      if (-x /opt/lam/bin/mpirun)     set mpirun=/opt/lam/bin/mpirun
      set mpirunops = "-c2c -O"
    else if (`egrep -c 'MPICHX|MPICH_DEBUG_ERRS' src/start.x` > 0) then # mpich
      if (-x /usr/lib/mpich/bin/mpirun) set mpirun=/usr/lib/mpich/bin/mpirun
      if (-x /opt/mpich/bin/mpirun)     set mpirun=/opt/mpich/bin/mpirun
      if ($?SGE_O_WORKDIR) then	# sge job
	set mpirunops = "-nolocal -machinefile $SGE_O_WORKDIR/machines-$JOB_NAME-$JOB_ID"
      else			# interactive run
	set mpirunops = "-nolocal" # or we get one CPU less
      endif
      else
      set mpirun 'Cannot_find_out_which_mpirun_to_use'
    endif
  endif

else if ($hn =~ nq* || $hn =~ ns*) then
  echo "Nordita cluster"
  if ($?PBS_NODEFILE) then
    echo "PBS job"
    cat $PBS_NODEFILE > lamhosts
    set local_disc = 1
    set one_local_disc = 0
    set local_binary = 0
    set notserial_procN = 1
  else
    echo "Non-PBS, running on `hostname`"
    echo `hostname` > lamhosts
  endif
  # lamboot with logging to file 
  echo "lambooting.."
  echo 'pidof lamd | xargs ps -lfwww :' > lamboot.log
  pidof lamd | xargs ps -lfwww >> lamboot.log
  echo '------------------------------------------------' >> lamboot.log
  lamboot -v lamhosts >>& lamboot.log # discard output, or auto-test
                                      # mysteriously hangs on Nq0
  echo '------------------------------------------------' >> lamboot.log
  echo 'pidof lamd | xargs ps -lfwww :' >> lamboot.log
  pidof lamd | xargs ps -lfwww >> lamboot.log
  #
  set booted_lam = 1
  echo "lamnodes:"
  lamnodes
  set mpirun = /opt/lam/bin/mpirun
  set mpirun = /usr/bin/mpirun
  set mpirunops = "-O -c2c -s n0 -x LD_ASSUME_KERNEL=2.4.1"#Fix bug in Redhat 9
  if ($local_disc) then
    setenv SCRATCH_DIR "/var/tmp"
  endif

else if (($hn =~ s[0-9]*p[0-9]*) || ($hn =~ 10_[0-9]*_[0-9]*_[0-9]*)) then
  echo "Horseshoe cluster"
  if ($mpi) then
    if ($?RUNNINGMPICH) then
      set mpirunops = "-machinefile $PBS_NODEFILE"
      set mpirun = /usr/local/lib/MPICH/bin/mpirun
      set start_x=$PBS_O_WORKDIR/src/start.x
      set run_x=$PBS_O_WORKDIR/src/run.x
    else
      echo "Using LAM-MPI"
      cat $PBS_NODEFILE > lamhosts
      lamboot -v lamhosts
      set booted_lam = 1
      echo "lamnodes:"
      lamnodes
      # set mpirunops = "-O -c2c -s n0 N -v" #(direct; should be faster)
      # set mpirunops = "-O -s n0 N -lamd" #(with debug options on; is slower)
      set mpirunops = "-O -s n0 N -ger" #(with debug options on; is slower)
      set mpirun = mpirun
    endif
    setenv SCRATCH_DIR /scratch
    set local_disc     = 1
    set one_local_disc = 0
    set remote_top     = 1
    set local_binary   = 1
    setenv SSH rsh 
    setenv SCP rcp
  else # (no MPI)
    echo "Batch job: non-MPI single processor run"
    cat $PBS_NODEFILE > lamhosts
    lamboot -v lamhosts 
    set booted_lam = 1
    echo "lamnodes:" 
    lamnodes
    set mpirunops = ''
    set mpirun = ''
  endif

else if ($hn =~ giga*) then
  echo "Giga Cluster"
  setenv SCRATCH_DIR /work
  if ($?JOB_ID) then
# Need to find out/configure where the SGE is writing the node list
    if (-e $HOME/.score/ndfile.$JOB_ID) then
      set local_disc=1
    else
      echo "WARNING: Cannot find ~/.score/ndfile.$JOB_ID, continuing without local disk access"
      set local_disc=0
    endif
  else
    set local_disc=0
  endif

else if (($hn =~ copson*.st-and.ac.uk) || ($hn =~ comp*.st-and.ac.uk)) then
  echo "Copson Cluster - St. Andrews"
  if ($?PE) then                            # Are we running under SGE?   
    if ($PE =~ gm) then                    # Using Myrinet?
      setenv SSH /usr/bin/rsh
      setenv SCP /usr/bin/rcp
      cat $PE_HOSTFILE | sed 's/\([[:alnum:].-]*\)\ \([0-9]*\).*/for ( i=0 \; i < 2 \; i++ ){print "\1\\n"};/' | bc > hostfile
      set mpirun = /usr/local/mpich-gm_INTEL/bin/mpirun 
      set mpirunops = "-local -machinefile hostfile"
      setenv SCRATCH_DIR `cat $TMPDIR/scratch`
      set local_disc=1     
      set one_local_disc=0
      set nprocpernode=2
      # Hack to give common scratch space path on each node
      #foreach host ($nodelist)
      #   $SSH $host "rm -rf $SCRATCH_DIR; ln -s ${dollar}TMPDIR $SCRATCH_DIR; ls -lR /tmp/pencil*" 
      #end
      echo '--------------- MPI_HOSTFILE ----------------'
      cat hostfile 
      echo '----------- MPI_HOSTFILE - END --------------'
    else if ($PE =~ score) then             # Using SCore?
      #set mpirunops = "-wait -F $HOME/.score/ndfile.$JOB_ID -e /tmp/scrun.$JOB_ID"
      #echo '--------------- PE_HOSTFILE ----------------'
      #cat $PE_HOSTFILE
      #echo '----------- PE_HOSTFILE - END --------------'
      set mpirunops = "-wait -F $PE_HOSTFILE -e $TMPDIR/scrun.$JOB_ID"
      set mpirun = /opt/score/bin/scout 
      #setenv SCRATCH_DIR /scratch/$JOB_ID
      setenv SCRATCH_DIR $TMPDIR
      set local_disc=1
      set one_local_disc=0
    endif
  else
      echo $hn >! hostfile
      set mpirun = /usr/local/mpich-gm_INTEL/bin/mpirun 
      set mpirunops = "-local -machinefile hostfile"
     set local_disc=0
  endif

else if ($hn =~ obelix) then
  echo "Obelix Cluster - Calgary"
  set local_disc = 0
  set one_local_disc = 0
  if ($?PBS_NODEFILE) then
    cp $PBS_NODEFILE lamhosts
    lamboot -v lamhosts
    set booted_lam = 1
    set mpirun = /opt/intel/compiler70/lam/bin/mpirun
  endif

else if ($hn == rasmussen) then
  echo "Rasmussen"
  limit stacksize unlimited
  set mpirun = 'mpirun'
  set mpirunops = ''
  set npops = ''
  set ncpus = 1

else if ($hn == hwwsr8k) then
  echo "Hitachi in Stuttgart"
  set nnodes = `echo "precision=0; ($ncpus-1)/8+1" | bc` # Hitachi-specific
  set mpirun = mpiexec
  set mpirunops = "-p multi -N $nnodes"
  set x_ops = '-Fport(nmlist(2))'
  set local_disc = 1
  set one_local_disc = 1	# (the default anyway)
  set local_binary = 0
  if (($?SCRDIR) && (-d $SCRDIR)) then
    setenv SCRATCH_DIR "$SCRDIR"
  else
    echo 'NO SUCH DIRECTORY: $SCRDIR'="<$SCRDIR> -- ABORTING"
    kill $$			# full-featured suicide
  endif
  setenv SSH rsh 
  setenv SCP rcp

else if ($hn =~ hwwsx5*) then
  echo "NEC-SX5 in Stuttgart"
  set mpirun = mpiexec
  set mpirunops = ''
  set local_disc = 1
  set one_local_disc = 1        # (the default anyway)
  set local_binary = 0 
  if (($?SCRDIR) && (-d $SCRDIR)) then
    setenv SCRATCH_DIR "$SCRDIR"
  else
    echo 'NO SUCH DIRECTORY: $SCRDIR'="<$SCRDIR> -- ABORTING"
    kill $$                     # full-featured suicide
  endif
  setenv SSH rsh
  setenv SCP rcp

else if ($hn =~ morvern || $hn =~ renton || $hn =~ lanark) then
  echo "LAM MPI on Newcastle desktops)" 
  #NB: need set mpi here: egrep on Makefile.local fails if using ~/.adapt-mkfile.inc 
  set mpi = 1
  echo "lamnodes:" 
  lamnodes
  set mpirun = 'mpirun'
  set mpirunops = ''
  set npops = ''
  set local_disc = 0
  set one_local_disc = 0 
  if ($local_disc) then
    setenv SCRATCH_DIR /var/tmp
  endif

else if ($hn =~ cosmo) then
  set mpirunops = '-x NLSPATH'

else if ($hn =~ mhd) then
  echo "mhd node at Newcastle (alpha linux running LAM MPI)" 
  #echo `hostname` "cpu=4" >! lamhosts
  #lamboot -v lamhosts
  #set booted_lam = 1
  echo "lamnodes:" 
  lamnodes
  set nprocpernode = 4
  set mpirun = mpirun
  #set mpirunops = "-ssi rpi tcp -s n0 N -v"
  set mpirunops = "-v"
  set npops = ''
  set local_disc = 0
  set one_local_disc = 0 
  if ($local_disc) then
    setenv SCRATCH_DIR /var/tmp
  endif

else
  echo "Generic setup; hostname is <$hn>"
  if ($mpi) echo "Use mpirun"
  set mpirun = mpirun
endif

## MPI specific setup
if ($mpi) then
  # Some mpiruns need special options
  if ($mpirun =~ *mpirun*) then
    set npops = "-np $ncpus"
  else if ($mpirun =~ *mpiexec*) then
    set npops = "-n $ncpus"
  else if ($mpirun =~ *scout*) then
    set nnode = `expr $NSLOTS - 1`
    set nprocpernode = `expr $ncpus / $nnode` 
    set npops = "-nodes=${nnode}x${nprocpernode}"
  else if ($mpirun =~ *poe*) then
    set nprocpernode = 1
    set npops = "-procs $ncpus"
  else
    echo "getconf.csh: No clue how to tell $mpirun to use $ncpus nodes"
  endif

else # no MPI

  echo "Non-MPI version"
  set mpirun = ''
  set mpirunops = ''
  set npops = ''
  set ncpus = 1
  set nprocpernode = 1

endif
## End of machine specific settings
## ------------------------------

# Determine data directory (defaults to `data')
if (-r datadir.in) then
  set datadir = `cat datadir.in | sed 's/ *\([^ ]*\).*/\1/'`
else
  set datadir = "data"
endif
echo "datadir = $datadir"

# If local disc is used, write name into $datadir/directory_snap.
# This will be read by the code, if the file exists.
# Remove file, if not needed, to avoid confusion.
if (-f $datadir/directory_snap) rm $datadir/directory_snap
if ($local_disc) then
  echo $SCRATCH_DIR >$datadir/directory_snap
endif

if ($local_binary) then
  set start_x = $SCRATCH_DIR/start.x
  set run_x   = $SCRATCH_DIR/run.x
endif


# Created subdirectories on local scratch disc (start.csh will also create
# them under $datadir/)
set subdirs = ("allprocs" "averages" "idl")
set procdirs = \
    `printf "%s%s%s\n" "for(i=0;i<$ncpus;i++){" '"proc";' 'i; }' | bc`
if ($local_disc) then
  echo "Creating directory structure on scratch disc(s)"
  foreach host ($nodelist)
    $SSH $host "if (! -e $SCRATCH_DIR ) mkdir -p $SCRATCH_DIR; cd $SCRATCH_DIR; mkdir -p $procdirs $subdirs" 
  end
endif

# Apply the SGI namelist read fix if running IRIX
set os = `uname -s`
if ($os =~ IRIX*) then
  touch SGIFIX
else
  rm -f SGIFIX
endif

# Propagate current pid to copy-snapshots:
setenv PARENT_PID $$

# Wrap up nodelist as (scalar, colon-separateds) environment variable
# NODELIST for transport to sub-processes.
setenv NODELIST `echo $nodelist | perl -ne 'print join(":",split(/\s/,$_)),"\n"'`

if ($debug) then
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
  echo '$mpirun       	= ' "<$mpirun>"
  echo '$mpirunops    	= ' "<$mpirunops>"
  echo '$x_ops        	= ' "<$x_ops>"
  echo '$NODELIST     	= ' "<$NODELIST>"
  echo '$SSH          	= ' "<$SSH>"
  echo '$SCP          	= ' "<$SCP>"
  echo '$PARENT_PID   	= ' "<$PARENT_PID>"
endif

exit

# End of file getconf.csh
