#!/bin/csh

# Name:   getconf.csh
# Author: wd (Wolfgang.Dobler@ncl.ac.uk)
# Date:   16-Dec-2001
# $Id: getconf.csh,v 1.174 2006-09-13 16:22:39 ajohan Exp $
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
  echo "(if it is left over from a crash, remove it by hand: rm LOCK)"
  echo ""
  echo "Will *not* start in this directory, as code may already be running."
  echo "Checking for NEWDIR file to tell us to run somewhere else:"
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
      echo "Found NEWDIR file: redirecting to new directory:"
      echo `pwd`
      echo "... wait 20 sec, to allow manual escape if we have produced a loop!"
      echo
      # In new run directory, need to check for yet another possible LOCK file
      # Sleep for 20 sec, to allow manual escape if we have produced a loop!
      sleep 20
      goto newdir
    endif
  else
    echo "No NEWDIR file -- exiting"
    echo
    # exit                        # won't work in a sourced file
    (sleep 1; kill -KILL $$ >& /dev/null) &	  # schedule full-featured suicide
    kill -TERM $$	  	  # .. but try exiting in civilized manner
  endif  
endif
# The LOCK is ours:
if (! -e "NEVERLOCK") touch LOCK

# Are we running the MPI version?
#set mpi = `egrep -c '^[ 	]*MPICOMM[ 	]*=[ 	]*mpicomm' src/Makefile`
set mpi = `fgrep -c 'mpicomm_init: MPICOMM neighbors' src/start.x`
# Determine number of CPUS
set ncpus = `perl -ne '$_ =~ /^\s*integer\b[^\\!]*ncpus\s*=\s*([0-9]*)/i && print $1' src/cparam.local`
echo "$ncpus CPUs"

# Location of executables and other default settings; can be overwritten
# below
set start_x = "src/start.x"
set run_x   = "src/run.x"
set x_ops = ""         # arguments to both start.x and run.x
set mpirun = 'mpirun'


# Check if experimental copy-snapshots is used
set copysnapshots="copy-snapshots"
set lcopysnapshots_exp=0
if { ( grep lcopysnapshots_exp=T start.in >& /dev/null ) } then
  set copysnapshots="copy-snapshots_exp"
  set lcopysnapshots_exp=1
endif

# Check for particles
set lparticles=0
if { ( grep particles_init_pars start.in >& /dev/null ) } set lparticles=1

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
set mpirunops  = ''  # options before -np $ncpus
set mpirunops2 = ''  # options after -np $ncpus

# Try to reconstruct submit host (to distinguish between the different
# clusters that have nodes called `node7', etc). Possibly, $masterhost and
# $masternode could be merged, but in fact these might be different host
# names
set masterhost = ''
if ($?PBS_O_HOST) then
  if ($PBS_O_HOST =~ obelix*) set masterhost = 'obelix'
endif
if ($?PBS_JOBID) then
  if ($PBS_JOBID =~ *.obelix*) set masterhost = 'obelix'
endif

# Get list of nodes; filters lines such that it would also handle
# machines.XXX files or lam-bhost.der, although this is hardly necessary.
if ($?PBS_NODEFILE) then
  if ($debug) echo "PBS job"
  set nodelist = `cat $PBS_NODEFILE | grep -v '^#' | sed 's/:[0-9]*//'`
else if ($?LOADL_PROCESSOR_LIST) then
  set nodelist = ""
  set lastproc = ""
  foreach proc ($LOADL_PROCESSOR_LIST)
    if ("$proc" != "$lastproc") then
      set nodelist = "$nodelist $proc"
    endif
    set lastproc = "$proc"
  end
  unset proc
  unset lastproc
else if ($?PE_HOSTFILE) then
  if ($debug) echo "SGE Parallel Environment job - $PE"
  set nodelist = `cat $PE_HOSTFILE | grep -v '^#' | sed 's/\ .*//'`
else if ($?JOB_ID) then
  if (-e $HOME/.score/ndfile.$JOB_ID) then
    if ($debug) echo "Scout job"
    set nodelist = `cat $HOME/.score/ndfile.$JOB_ID | grep -v '^#' | sed 's/:[0-9]*//'`
  else
    # E.g. for wd running SGE on Kabul cluster
    echo "Apparently not parallel environment job"
    if ($debug) echo "Setting nodelist to ($hn)"
    set nodelist = ("$hn")
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
  echo "MHD machine - St Andrews "
  set mpirun = "dmpirun"

else if ($hn =~ *.kis.uni-freiburg.de) then
  echo "KIS machines - Freiburg"
  set mpirun = /opt/local/mpich/bin/mpirun

else if ($hn =~ hamlet) then
  echo "hamlet pc, Heidelberg (formerly in Copenhagen)"
  set mpirun = /home/ajohan/mpich/bin/mpirun

else if ($hn =~ giga[0-9][0-9].ncl.ac.uk) then
  echo "Newcastle e-Science Cluster"
  if ($?PBS_JOBID) then
  set mpirun = mpiexec
  set mpirunops = 
  #set mpirun = mpirun
  #set mpirunops = "-machinefile $PBS_NODEFILE -O -ssi rpi tcp -s n0 -x LD_ASSUME_KERNEL=2.4.1" #Fix bug in Redhat 9
  #set mpirunops = "-machinefile $PBS_NODEFILE" 
    set one_local_disc = 0
    set local_binary = 0
    #setenv SCRATCH_DIR /work/$PBS_JOBID 
#    echo "SGE job"
#    set local_disc = 1
#    set one_local_disc = 0
#    set local_binary = 0
#
#  #  cat $PE_HOSTFILE | sed 's/\([[:alnum:].-]*\)\ \([0-9]*\).*/for ( i=0 \; i < 2 \; i++ ){print "\1\\n"};/' | bc > hostfile
#  #  set nodelist = `cat hostfile`
#
#    if ($PE =~ mpi2) then
#      set nprocpernode = 2
#    else if ($PE =~ mpi) then
#      set nprocpernode = 4
#    endif
  else
    echo "Non-SGE, running on `hostname`"
    if ( -e MANUALMPI ) then
      set nprocpernode = 2
      set local_disc = 0
      set one_local_disc = 0
      set local_binary = 0
      setenv JOB_ID MANUALMPI
    endif
  endif
#  set mpirun = mpirun
#  set mpirunops = "-O -ssi rpi tcp -s n0 -x LD_ASSUME_KERNEL=2.4.1" #Fix bug in Redhat 9
#  set local_disc = 1
  if ($local_disc) then
    #setenv SCRATCH_DIR `cat $TMPDIR/scratch` 
    set remove_scratch_root = 1
  endif

else if ($hn =~ giga[0-9][0-9]) then
  echo "Nordita cluster - Copenhagen"
  if ($?PBS_NODEFILE) then
    echo "PBS job"
    cat $PBS_NODEFILE >! lamhosts
    set local_disc = 1
    set one_local_disc = 0
    set local_binary = 0
  else
    echo "Non-PBS, running on `hostname`"
    echo `hostname` >! lamhosts
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
  echo "IBM machine - Aarhus"
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

else if ($hn =~ p690*) then
  echo "SP1600 - CSC, Espoo, Finland (IBM with AIX UNIX)"
  set mpirunops = "-euilib us -shared_memory yes"
  set mpirunops = 
  set mpirun = poe
  set local_disc = 0
  set one_local_disc = 0
  set local_binary = 0
#  setenv SSH ssh
#  setenv SCP scp
#  setenv LANG en_US
#  setenv SCRATCH_DIR /scratch/${USER}
#  set masternode=sp.sp4
#  echo "Setting master node to sp.sp4, the only node that is accesible by rsh"

else if ($hn =~ sp[0-9]*) then
  echo "SP4 - CINECA, Bologna (IBM with AIX UNIX)"
  set mpirun = poe
  set local_disc = 0
  set one_local_disc = 1
  set local_binary = 0
#  setenv SSH ssh
#  setenv SCP scp
  setenv LANG en_US
  setenv SCRATCH_DIR /scratch/${USER}
#  set masternode=sp.sp4
  echo "Setting master node to sp.sp4, the only node that is accesible by rsh"

else if ($hn =~ node*.steno.dcsc.ku.dk) then
  echo "DCSC cluster - Copenhagen (Denmark)"
  setenv SCRATCH_DIR /scratch/$LOADL_STEP_ID
  set hostfile = "$LOADL_STEP_INITDIR/hostfile"
  echo $LOADL_PROCESSOR_LIST | tr ' ' '\n' > $hostfile
  set mpirunops      = "-hostfile $hostfile"
  set mpirun         = "/usr/local/topspin/mpi/mpich/bin/mpirun_ssh"
  set local_disc     = 1
  set one_local_disc = 0
  set remote_top     = 1
  set local_binary   = 0
  set nprocpernode   = 4

else if (($hn =~ columbia*) || ($hn =~ cfe*)) then
  echo "Columbia cluster"
  set mpirun = mpirun
  set start_x=$PBS_O_WORKDIR/src/start.x
  set run_x=$PBS_O_WORKDIR/src/run.x
  set local_disc = 0

else if ($hn =~ node[0-9]*) then
  echo "CLX - CINECA, Bologna (IBM Linux Cluster)"
  if ( $?PBS_JOBID ) then 
    echo "Running job: $PBS_JOBID"
  endif
  set mpirun = mpiexec
  set mpirunops =
  set local_disc = 0
  set one_local_disc = 1
  set local_binary = 0

else if ($hn =~ psi*) then
  echo "RZG - Garching (IBM pSeries Regatta with AIX UNIX)"
  set mpirun = mpiexec
  set local_disc = 1
  set one_local_disc = 1
  set local_binary = 0
  setenv SSH rsh
  setenv SCP rcp
  setenv LANG en_US
  setenv SCRATCH_DIR /ptmp/$USER
  set masternode=psi24
  echo "Setting master node to psi24, the only node that is accesible by rsh"

else if ($hn =~ sepeli.csc.fi) then
  echo "Sepeli - CSC, Espoo, Finland"
  set mpirunops = ''
  set mpirun = 'mpirun'
  set npops = "-np $ncpus"
  set local_disc = 0
  set one_local_disc = 0
  set local_binary = 0
else if ($hn =~ compute-*.local) then
  echo "Sepeli Nodes - CSC, Espoo, Finland"
  set mpirunops = "-machinefile `ls -d /tmp/${JOB_ID}*`/machines"
  set mpirun = 'mpirun'
  set npops = "-np $ncpus"
  set local_disc = 0
  set one_local_disc = 0
  set local_binary = 0
#  setenv SCRATCH_DIR $TMPDIR

else if ($hn =~ corona*) then
  echo "Corona SunFire - CSC, Espoo, Finland"
  # generic fgrep -c fails
  set mpi = 1
  set tmpi = `fgrep -c 'nompicomm' src/Makefile.local`
  if ($tmpi == 1) set mpi = 0 
  set mpirun = 'mprun'
  set npops = "-np $ncpus"
  set mpirunops = ""

else if ( ($hn =~ node*.clusters.com) || ($hn =~ fire) ) then
  echo "Fire - Bergen"
  set mpirunops = ''
  set mpirun = 'mpirunpbs'
  set npops = "-np $ncpus"
#  set mpi = 1


else if ( ($hn =~ cincinnatus*) || ($hn =~ owen*) \
          || ($hn =~ master) || ($hn =~ node* && $masterhost == '') \
	  || ($hn =~ ns0*) ) then
  echo "KIS cluster or Ns0"
  if ($mpi) then
    # Choose appropriate mpirun version (LAM vs. MPICH)
    if (`fgrep -c lam_mpi src/start.x` > 0) then # lam
      if (-x /usr/lib/lam/bin/mpirun) set mpirun=/usr/lib/lam/bin/mpirun
      if (-x /opt/lam/bin/mpirun)     set mpirun=/opt/lam/bin/mpirun
      set mpirunops = "-boot -c2c -O"
    else if (`egrep -c 'MPICHX|MPICH_DEBUG_ERRS' src/start.x` > 0) then # mpich
      if (-x /usr/lib/mpich/bin/mpirun)   set mpirun=/usr/lib/mpich/bin/mpirun
      if (-x /opt/mpich/bin/mpirun)       set mpirun=/opt/mpich/bin/mpirun
      if (-x /opt/mpich/ch-p4/bin/mpirun) set mpirun=/opt/mpich/ch-p4/bin/mpirun
      if ($?SGE_O_WORKDIR) then	# sge job
	set mpirunops = "-nolocal -machinefile $SGE_O_WORKDIR/machines-$JOB_NAME-$JOB_ID"
      else			# interactive run
	set mpirunops = "-nolocal" # or we get one CPU less with older mpich
      endif
      if ($mpirun == '/opt/mpich/ch-p4/bin/mpirun') then
	# Don't use -nolocal with newer mpich versions, otherwise we are
	# not able to run on single-CPU test machines
	set mpirunops = "`echo $mpirunops | sed 's/ \?-nolocal \?/ /'`"
      endif
      else
      set mpirun 'Cannot_find_out_which_mpirun_to_use'
    endif
  endif

else if ($hn == "frontend") then
  echo "KIS Bagdad cluster - Freiburg"
  if ($mpi) then
    # Choose appropriate mpirun version (LAM vs. MPICH)
    if (`fgrep -c lam_mpi src/start.x` > 0) then # lam
      if (-x /usr/local/share/lam-mpi/bin/mpirun) set mpirun=/usr/local/share/lam-mpi/bin/mpirun
      if (-x /opt/lam/bin/mpirun)     set mpirun=/opt/lam/bin/mpirun
      set mpirunops = "-O"
    else if (`egrep -c 'MPICHX|MPICH_DEBUG_ERRS' src/start.x` > 0) then # mpich
      if (-x /usr/local/share/mpich-g95/bin/mpirun)   set mpirun=/usr/local/share/mpich-g95/bin/mpirun
    endif
  endif

else if ($hn =~ nq* || $hn =~ ns*) then
  echo "Nordita cluster - Copenhagen"
  if ($?PBS_NODEFILE) then
    echo "PBS job"
    cat $PBS_NODEFILE >! lamhosts
    set local_disc = 1
    set one_local_disc = 0
    set local_binary = 0
    set notserial_procN = 1
  else
    echo "Non-PBS, running on `hostname`"
    echo `hostname` >! lamhosts
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
  echo "Horseshoe cluster - Odense (Denmark)"
  if ($mpi) then
    if ($?RUNNINGMPICH) then
      set mpirunops = "-machinefile $PBS_NODEFILE"
      set mpirun = /usr/local/lib/MPICH/bin/mpirun
      set start_x=$PBS_O_WORKDIR/src/start.x
      set run_x=$PBS_O_WORKDIR/src/run.x
    else
      echo "Using LAM-MPI"
      cat $PBS_NODEFILE >! lamhosts
      lamboot -v lamhosts
      set booted_lam = 1
      echo "lamnodes:"
      lamnodes
      # set mpirunops = "-O -c2c -s n0 N -v" #(direct; should be faster)
      # set mpirunops = "-O -s n0 N -lamd" #(with debug options on; is slower)
      set mpirunops = "-O -s n0 N -ger" #(with debug options on; is slower)
      set mpirun = mpirun
    endif
    # use unique scratch directory name, just in case it wasn't cleaned up
    if ($?PBS_JOBID) then
      setenv SCRATCH_DIR /scratch/$PBS_JOBID
    else
      setenv SCRATCH_DIR /scratch
    endif
    set local_disc     = 1
    set one_local_disc = 0
    set remote_top     = 1
    set local_binary   = 1
    setenv SSH rsh 
    setenv SCP rcp
  else # (no MPI)
    echo "Batch job: non-MPI single processor run"
    cat $PBS_NODEFILE >! lamhosts
    lamboot -v lamhosts 
    set booted_lam = 1
    echo "lamnodes:" 
    lamnodes
    set mpirunops = ''
    set mpirun = ''
  endif

else if (($hn =~ copson*.st-and.ac.uk) || ($hn =~ comp*.st-and.ac.uk)) then
  echo "Copson Cluster - St. Andrews"
  if ($?PE) then                            # Are we running under SGE?   
    if ($PE =~ gm-test) then                    # Using Myrinet?
      setenv SSH /usr/bin/rsh
      setenv SCP /usr/bin/rcp
      cat $PE_HOSTFILE | sed 's/\([[:alnum:].-]*\)\ \([0-9]*\).*/for ( i=0 \; i < 2 \; i++ ){print "\1\\n"};/' | bc > hostfile
      set mpirun = /usr/local/mpich-gm_INTEL/bin/mpirun 
      echo "Setting mpirun... $mpirun"
      set mpirunops = "-local -machinefile $TMPDIR/machines"
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
    else if ($PE =~ gm) then                    # Using Myrinet?
      #setenv SSH ssh
      #setenv SCP scp
      setenv SSH /usr/bin/rsh
      setenv SCP /usr/bin/rcp
#      setenv MPIHOME /usr/local/mpich-gm-1.2.6..14/pgi-intel-7.1/bin
#      setenv PATH ${SGE_O_PATH}:${PATH} 
      set mpirun = /usr/local/mpi_wrappers/mpirun
#${MPIHOME}/mpirun
      set mpirunops = "-local -machinefile $TMPDIR/machines"

      setenv SCRATCH_DIR `cat $TMPDIR/scratch`
      set local_disc=1     
      set one_local_disc=0
      set nprocpernode=2
      echo '--------------- MPI_HOSTFILE ----------------'
      cat $TMPDIR/machines
      cat $TMPDIR/machines >! hostfile
      echo '----------- MPI_HOSTFILE - END --------------'
#      env | sort >! env.run
    else if ($PE =~ score) then             # Using SCore?
      #set mpirunops = "-wait -F $HOME/.score/ndfile.$JOB_ID -e /tmp/scrun.$JOB_ID"
      #echo '--------------- PE_HOSTFILE ----------------'
      #cat $PE_HOSTFILE
      #echo '----------- PE_HOSTFILE - END --------------'
      set mpirunops = "-wait -F $PE_HOSTFILE -e $TMPDIR/scrun.$JOB_ID"
      set mpirun = /opt/score/bin/scout 
      echo "Setting mpirun... $mpirun"
      #setenv SCRATCH_DIR /scratch/$JOB_ID
      setenv SCRATCH_DIR $TMPDIR
      set local_disc=1
      set one_local_disc=0
    endif
  else
      echo $hn >! hostfile
      set mpirun = /usr/local/mpich-gm_INTEL/bin/mpirun 
      echo "Setting mpirun... $mpirun"
      set mpirunops = "-local -machinefile hostfile"
     set local_disc=0
  endif

else if ($hn =~ obelix || \
          ($hn =~ node[0-9][0-9] && $masterhost =~ obelix)) then
  echo "Obelix Cluster - Calgary"
  set local_disc = 0
  set one_local_disc = 0
  if ($?PBS_NODEFILE) then
    cp $PBS_NODEFILE lamhosts
    lamboot -v lamhosts
    set booted_lam = 1
    # set mpirun = /opt/intel/compiler70/lam/bin/mpirun
    set mpirun = /opt/gnu/lam-7.0/bin/mpirun
    set mpirun = /opt/intel/lam-7.1.1/bin/mpiexec
    set mpirunops = '-boot'
  endif

else if ($hn == pencil) then
  echo "Pencil [sic] - Calgary"
  set mpirun = mpiexec
  set mpirunops = '-boot'

else if ($hn == rasmussen) then
  echo "Rasmussen - Copenhagen (SGI machine)"
  limit stacksize unlimited
  set mpirun = 'mpirun'
  set mpirunops = ''
  set npops = ''
  set ncpus = 1

else if ($hn == hwwsr8k) then
  echo "Hitachi - Stuttgart"
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
  echo "NEC-SX5 - Stuttgart"
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
  echo "LAM MPI - Newcastle desktops)" 
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

else if ($hn =~ hades) then
  echo "Hades - Wlad's laptop"
  set mpirun = 'mpiexec'
  set mpirunops = '-boot'

else if ($hn =~ kolmogorov) then
  echo "Kolmogorov on Wolfgang's desk"
  set mpirun = 'mpiexec'
  set mpirunops = '-boot'

else if ($hn =~ mhd) then
  echo "mhd node - Newcastle (alpha linux running LAM MPI)"
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

else if ($hn =~ opto[1-4]) then
  echo "opto[1-4] 4xAMD opteron procs@2190.948 MHz at MPIA"
  set nprocpernode = 4
  set mpirun = mpirun_rsh
  if ($ncpus <= 4) then
    echo `repeat $ncpus echo $hn` > hosts.list
  else
    set ncpus2=`echo $ncpus-4 | bc`
    if ($hn =~ opto4) set hn2=opto4
    if ($hn =~ opto3) set hn2=opto3
    echo `repeat 4 echo $hn; repeat $ncpus2 echo $hn2` > hosts.list
  endif
  set mpirun = '~/mpich/bin/mpirun'
  set mpirunops = '-machinefile hosts.list'
  setenv SSH rsh
  setenv SCP rcp
  setenv SCRATCH_DIR /var/tmp/$USER

else if ($hn =~ swift) then
  echo "swift laptop (Anders)"
  set mpirun = mpiexec
  setenv SSH ssh
  setenv SCP scp
  setenv SCRATCH_DIR /var/tmp/$USER

else if ($hn =~ lfc*) then
  echo "opteron cluster at MPIK with SGE queue."
  if ($#nodelist == 1) then
    echo "Apparently an interactive run."
    set nodelist = `repeat $ncpus echo $nodelist`
  endif
  echo "Writing host list to file hosts.list."
  echo $nodelist > hosts.list
  set nprocpernode = 1
  if ($mpi) then
    set local_disc = 1
    set one_local_disc = 0
  else
    set local_disc = 0
    set one_local_disc = 1
  endif
  set mpirun = mpirun
  set mpirunops = '-nolocal'
  set mpirunops2 = '-machinefile hosts.list'
  setenv SSH ssh_csh_aj
  setenv SCP scp
  setenv SCRATCH_DIR /var/tmp/$USER
#  setenv SGE_CWD_PATH `pwd`
#  setenv SGE_O_WORKDIR `pwd`
#  setenv JOB_SCRIPT `pwd`/start.csh

else if ($hn =~ rio* || $hn =~ pia*) then
  echo "opteron cluster at RZG with SGE queue."
  if ($#nodelist == 1) then
    echo "Apparently an interactive run."
    set nodelist = `repeat $ncpus echo $nodelist`
  else
    set nprocpernode = 1
    if ($mpi) then
      set local_disc = 1
      set one_local_disc = 0
    else
      set local_disc = 0
      set one_local_disc = 1
    endif
    set mpirun = /afs/ipp/amd64_sles9/soft/mvapich-0.9.5/gcc/bin/mpirun
    set mpirunops2 = '-hostfile $TMPDIR/machines'
    set nodelist=`cat $TMPDIR/machines`
    echo $nodelist
    setenv SSH 'ssh -x'
    setenv SCP scp
    setenv SCRATCH_DIR /var/tmp/$USER
    foreach node ($nodelist)
      echo $SSH $node 'killall start.x run.x'
      $SSH $node 'killall start.x run.x'
      echo $SSH $node '\rm -rf /var/tmp/'$USER'/*'
      $SSH $node '\rm -rf /var/tmp/$USER/*'
    end
  endif

else if ($hn =~ nut.uppmax.uu.se) then
  echo "NUT - primary node of RA, Uppsala"
  #http://www.uppmax.uu.se/ComputerSystems/Ra/Ra.html  
  set mpirun = /opt/scali/bin/mpimon
  set mpi_suffix = "$nodelist $nprocpernode"
  setenv SCRATCH_DIR data/scratch

else if ($hn =~ ra*) then
  #copied these lines from Ander's rio & pia, 
  #which also use the portland compilers.
  echo "RA Linux cluster at Uppmax with SGE queue - Uppsala."
  #http://www.uppmax.uu.se/ComputerSystems/Ra/Ra.html  
  set mpirun = /opt/scali/bin/mpimon
  setenv SCRATCH_DIR data/scratch
  if ($#nodelist == 1) then
    echo "Apparently an interactive run."
    set nodelist = `repeat $ncpus echo $nodelist`
    set mpi_suffix = "$nodelist"
  else
    set nodelist = `cat $TMPDIR/machines`
    echo $nodelist
    set nprocpernode = 1
    set mpi_suffix="$nodelist"
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
  else if ($mpirun =~ *mpimon*) then
    set npops = "-stdin all -inherit_limits"
    set x_ops = "-- $mpi_suffix"
  else if ($mpirun =~ *scout*) then
    set nnode = `expr $NSLOTS - 1`
    set nprocpernode = `expr $ncpus / $nnode` 
    set npops = "-nodes=${nnode}x${nprocpernode}"
  else if ($mpirun =~ *poe*) then
    set nprocpernode = 1
    set x_ops = "$mpirunops -procs $ncpus"
    set mpirunops =
    set npops = ""
  else
    echo "getconf.csh: No clue how to tell $mpirun to use $ncpus nodes"
  endif

else # no MPI

  echo "Non-MPI version"
  set mpirun = ''
  set mpirunops = ''
  set mpirunops2 = ''
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

# Propagate current pid to copy-snapshots:
setenv PARENT_PID $$

# Make SCRATCH_DIR unique, so different
# jobs running simultaneously will not interfere with each other. 
if ($?JOB_ID) then
  setenv SCRATCH_DIR ${SCRATCH_DIR}/pencil-$JOB_ID
else
  setenv SCRATCH_DIR ${SCRATCH_DIR}/pencil-$PARENT_PID
endif

# If local disc is used, write name into $datadir/directory_snap.
# This will be read by the code, if the file exists.
# Remove file, if not needed, to avoid confusion.
if (-f $datadir/directory_snap) rm $datadir/directory_snap
if ($local_disc) then
  echo $SCRATCH_DIR >! $datadir/directory_snap
endif

if ($local_binary) then
  set start_x = $SCRATCH_DIR/start.x
  set run_x   = $SCRATCH_DIR/run.x
endif

# Set nodelist to first node (or master node, if specified in the machine
# dependent part), if common scratch disc.
if ($one_local_disc) then
  if ($?masternode) then
    set nodelist=$masternode
    echo "Setting nodelist to master node $masternode"
  else
    set nodelist=$nodelist[1]
    echo "Setting nodelist to node $nodelist"
  endif
endif

# Create subdirectories on local scratch disc (start.csh will also create
# them under $datadir/)
set subdirs = ("allprocs" "averages" "idl")
set procdirs = \
    `perl -e 'for $i (0..'"$ncpus"'-1) { print "proc$i\n"}'`
if ($local_disc) then
  if ($one_local_disc) then
    echo "Creating directory structure on common scratch disc"
  else
    echo "Creating directory structure on local scratch disc(s)"
  endif
  #set command = \'"if (! -e $SCRATCH_DIR ) mkdir -p $SCRATCH_DIR; cd $SCRATCH_DIR; mkdir -p $procdirs $subdirs"\'
  foreach host ($nodelist)
    $SSH $host "if (! -e $SCRATCH_DIR ) mkdir -p $SCRATCH_DIR; cd $SCRATCH_DIR; mkdir -p $procdirs $subdirs" 
    #TH: Not everybody uses csh!
    #$SSH $host /bin/csh -c $command
  end
endif

# Apply the SGI namelist read fix if running IRIX
set os = `uname -s`
if ($os =~ IRIX*) then
  touch SGIFIX
else
  rm -f SGIFIX
endif

# Wrap up nodelist as (scalar, colon-separated) environment variable
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
  echo '$masterhost	= ' "<$masterhost>"
  echo '$mpirun       	= ' "<$mpirun>"
  echo '$mpirunops    	= ' "<$mpirunops>"
  echo '$mpirunops2   	= ' "<$mpirunops2>"
  echo '$x_ops        	= ' "<$x_ops>"
  echo '$NODELIST     	= ' "<$NODELIST>"
  echo '$SSH          	= ' "<$SSH>"
  echo '$SCP          	= ' "<$SCP>"
  echo '$PARENT_PID   	= ' "<$PARENT_PID>"
  echo '$copysnapshots  = ' "<$copysnapshots>"
  echo '$particles      = ' "<$lparticles>"
endif

exit

# End of file getconf.csh
