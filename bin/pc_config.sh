#!/bin/sh

# Name:   pc_config
# Author: Antony Mee (A.J.Mee@ncl.ac.uk)
# Date:   05-Apr-2004
# $Id$
#
# Description:
#  Initiate some variables related to MPI and the calling sequence, and do
#  other setup needed by both pc_start and pc_run execution scripts
#
# History:
#  Reworked csh version by WD (Wolfgang.Dobler@kis.uni-freiburg.de)
#   (see tagged cvs revision: cvs up -r sh_version getconf.csh)
#

#source $PENCIL_HOME/bin/pc_functions.sh
. pc_functions.sh

# Check we are actually in a run directory
check_is_run_directory
# Set up PATH for people who don't include $PENCIL_HOME/bin by default
[ -n "$PENCIL_HOME" ] && export PATH=${PATH}:${PENCIL_HOME}/bin
# Save working directory for other scripts we call
export PENCIL_WORKDIR=`pwd`
# Are we running the MPI version?
determine_mpi
# Determine number of CPUS
determine_ncpus

##--------------------------- DEFAULT SETTINGS ----------------------------
##---------------- can be overwritten per machine below -------------------
# Location of executables
start_x="src/start.x"
run_x="src/run.x"
x_ops=''         # arguments to both start.x and run.x

# MPI executable
mpirun='mpirun'
mpirunops=''

# Settings for machines with local data disks
# local_disc     = yes means processes uses their own local scratch disc(s) for
#                      large files like var.dat, VARN, etc.
# one_local_disc = yes means one common scratch disc for all processes
# local_binary   = yes  means copy executable to local scratch area of master node
# remote_top     = yes means get `top' output in regular intervals
local_disc=no
one_local_disc=yes              # probably the more common case
local_binary=no
remove_scratch_root=no
export SCRATCH_DIR=/scratch

# Does one require top snapshots to be grabbed from all the nodes periodically?
remote_top=no

# Which tools should be used for remote copy and shell access?
export SSH=ssh
export SCP=scp

# Any lam daemons to shutdown at the end of run.csh?
booted_lam=no

##----------------------- END OF DEFAULT SETTINGS -------------------------
##-------------------------------------------------------------------------

echo `uname -a`
hn=`uname -n`

# Get list of nodes; filters lines such that it would also handle
# machines.XXX files or lam-bhost.der, although this is hardly necessary.
determine_nodelist

##---------------------- MACHINE SPECIFIC SETTINGS ------------------------
##---------------- settings based upon the host name ----------------------
##-- ishost uses it's parameter as a regular expression to match against --
##-- the current machines hostname                                       --

if ishost "mhd*.st-and.ac.uk" ; then
  echo "St Andrews machine"
  mpirun="dmpirun"

elif ishost "kis.uni-freiburg.de" ; then
  mpirun=/opt/local/mpich/bin/mpirun

elif ishost "giga[0-9][0-9].ncl.ac.uk" ; then
  echo "Newcastle e-Science Cluster"
  queue_submit() # Override queue submission shell function
  {
    qsub -pe mpi2 `expr \(\( $ncpus + 1 \) / 2 \) * 2`  -N $2 $1
  }
  export LD_LIBRARY_PATH=/addon/shared/intel/compiler70/ia32/bin:$LD_LIBRARY_PATH
  if [ -n "$PE" ]; then
    echo "SGE job"
    local_disc=yes
    one_local_disc=no
    local_binary=no

  #  cat $PE_HOSTFILE | sed 's/\([[:alnum:].-]*\)\ \([0-9]*\).*/for ( i=0 \; i < 2 \; i++ ){print "\1\\n"};/' | bc >hostfile
  #  set nodelist = `cat hostfile`

    if [ "$PE" = "mpi2" ]; then
      nprocpernode=2
    elif [ "$PE" = "mpi" ]; then
      nprocpernode=4
    fi
  else
    echo "Non-SGE, running on $hn"
  fi
  mpirun=/addon/shared/lam/bin/mpirun
  mpirunops="-O -ssi rpi tcp -s n0 -x LD_ASSUME_KERNEL=2.4.1" #Fix bug in Redhat 9
  if [ "$local_disc" = "yes" ]; then
    export SCRATCH_DIR=/work/$JOB_ID
    remove_scratch_root=no
  fi

elif ishost "giga[0-9][0-9]" ; then
  echo "Nordita cluster"
  if [ [ -n $PBS_NODEFILE ] && [ -e $PBS_NODEFILE ] ]; then
    echo "PBS job"
    cat $PBS_NODEFILE >lamhosts
    local_disc=yes
    one_local_disc=no
    local_binary=no
  else
    echo "Non-PBS, running on `hostname`"
    echo `hostname` >lamhosts
  fi
  # lamboot with logging to file
  echo "lambooting.."
  echo 'pidof lamd | xargs ps -lfwww :' >lamboot.log
  pidof lamd | xargs ps -lfwww 2>&1 >>lamboot.log
  echo '------------------------------------------------' >>lamboot.log
  lamboot -v lamhosts 2>&1 >>lamboot.log # discard output, or auto-test
                                      # mysteriously hangs on Nq0
  echo '------------------------------------------------' >>lamboot.log
  echo 'pidof lamd | xargs ps -lfwww :' >>lamboot.log
  pidof lamd | xargs ps -lfwww 2>&1 >>lamboot.log
  #
  booted_lam=yes
  echo "lamnodes:"
  lamnodes
  mpirun=/opt/lam/bin/mpirun
  mpirun=/usr/bin/mpirun
  mpirunops="-O -c2c -s n0 -x LD_ASSUME_KERNEL=2.4.1"   #Fix bug in Redhat 9
  [ "$local_disc" = "yes" ] && export SCRATCH_DIR="/var/tmp"

elif ishost sleipner  || ishost fenris || ishost hugin ishost munin ; then
  mpirun = /usr/bin/mpiexec
  local_disc=yes
  one_local_disc=yes
  local_binary=no
  if [ -d $SCRDIR ]; then
    export SCRATCH_DIR="$SCRDIR"
  else
    echo 'NO SUCH DIRECTORY: $SCRDIR'="<$SCRDIR> -- ABORTING"
    kill $$                     # full-featured suicide
  fi
  export SSH=rsh
  export SCP=rcp
  export LANG=en_US

elif ishost cincinnatus || ishost owen || ishost master || ishost node || ishost ns0 ; then
  if [ "$mpi" = "yes" ]; then
    # Choose appropriate mpirun version (LAM vs. MPICH)
    if [ `fgrep -c lam_mpi src/start.x` -gt 0 ]; then # lam
      mpirun=/usr/lib/lam/bin/mpirun
      mpirunops="-c2c -O"
    elif [ `egrep -c 'MPICHX|MPICH_DEBUG_ERRS' src/start.x` -gt 0 ]; then # mpich
      [ -x /usr/lib/mpich/bin/mpirun ] && mpirun=/usr/lib/mpich/bin/mpirun
      [ -x /opt/mpich/bin/mpirun ] && mpirun=/opt/mpich/bin/mpirun
      if [ -d $SGE_O_WORKDIR ]; then    # sge job
        mpirunops="-nolocal -machinefile $SGE_O_WORKDIR/machines-$JOB_NAME-$JOB_ID"
      else                      # interactive run
            mpirunops="-nolocal" # or we get one CPU less
      fi
      else
      mpirun='Cannot_find_out_which_mpirun_to_use'
    fi
  fi

elif ishost "nq[0-9]*" || ishost "ns[0-9]*" ; then
  echo "Nordita cluster"
  if [ -e $PBS_NODEFILE ]; then
    echo "PBS job"
    cat $PBS_NODEFILE >lamhosts
    local_disc=yes
    one_local_disc=no
    local_binary=no
  else
    echo "Non-PBS, running on `hostname`"
    echo `hostname` >lamhosts
  fi
  # lamboot with logging to file
  echo "lambooting.."
  echo 'pidof lamd | xargs ps -lfwww :' &>lamboot.log
  pidof lamd | xargs ps -lfwww 2>&1 >>lamboot.log
  echo '------------------------------------------------' 2>&1 >>lamboot.log
  lamboot -v lamhosts 2>&1 >>lamboot.log # discard output, or auto-test
                                    # mysteriously hangs on Nq0
  echo '------------------------------------------------' 2>&1 >>lamboot.log
  echo 'pidof lamd | xargs ps -lfwww :' >>lamboot.log
  pidof lamd | xargs ps -lfwww 2>&1 >>lamboot.log
  #
  booted_lam=yes
  echo "lamnodes:"
  lamnodes
  mpirun=/opt/lam/bin/mpirun
  mpirun=/usr/bin/mpirun
  mpirunops="-O -c2c -s n0 -x LD_ASSUME_KERNEL=2.4.1"  #Fix bug in Redhat 9
  [ "$local_disc" = "yes" ] && export SCRATCH_DIR="/var/tmp"

elif ishost "s[0-9]*p[0-9]*" || ishost "10_[0-9]*_[0-9]*_[0-9]*" ; then
  echo "Horseshoe cluster"
  if [ "$mpi" = "yes" ]; then
#ajwm  NOT SURE IF THIS IS THE CORRECT CONDITION?...
    if [ $RUNNINGMPICH -ne 0 ]; then
      mpirunops="-machinefile $PBS_NODEFILE"
      mpirun=/usr/local/lib/MPICH/bin/mpirun
      start_x=$PBS_O_WORKDIR/src/start.x
      run_x=$PBS_O_WORKDIR/src/run.x
    else
      echo "Using LAM-MPI"
      cat $PBS_NODEFILE >lamhosts
      lamboot -v lamhosts
      booted_lam=yes
      echo "lamnodes:"
      lamnodes
      # set mpirunops = "-O -c2c -s n0 N -v" #(direct; should be faster)
      # set mpirunops = "-O -s n0 N -lamd" #(with debug options on; is slower)
      mpirunops="-O -s n0 N -ger" #(with debug options on; is slower)
      mpirun=mpirun
    fi
    export SCRATCH_DIR=/scratch
    local_disc=yes
    one_local_disc=no
    remote_top=yes
    local_binary=yes
    export SSH=rsh
    export SCP=rcp
  else # (no MPI)
    echo "Batch job: non-MPI single processor run"
    cat $PBS_NODEFILE >lamhosts
    lamboot -v lamhosts
    booted_lam=yes
    echo "lamnodes:"
    lamnodes
    mpirunops=''
    mpirun=''
  fi

elif ishost "giga[0-9]*" ; then
  echo "Giga Cluster"
  export SCRATCH_DIR=/work
# Need to find out/configure where the SGE is writing the node list
  if [ [ -n "$JOB_ID" ] && [ -e "$HOME/.score/ndfile.$JOB_ID" ] ]; then
    local_disc=yes
  else
    echo "WARNING: Cannot find ~/.score/ndfile.$JOB_ID, continuing without local disk access"
    local_disc=no
  fi

elif ishost "copson\.st-and\.ac\.uk" || ishost "comp[0-9]*.st-and.ac.uk" ; then
  echo "Copson Cluster - St. Andrews"
  queue_submit() # Override queue submission shell function
  {
    qsub -pe gm $ncpus -N $2 $1
  }
   if [ -n "$PE" ]; then                            # Are we running under SGE?
    if [ "$PE" = "gm" ]; then                    # Using Myrinet?
      export SSH=/usr/bin/rsh
      export SCP=/usr/bin/rcp
      cat $PE_HOSTFILE | sed 's/\([[:alnum:].-]*\)\ \([0-9]*\).*/for ( i=0 \; i < 2 \; i++ ){print "\1\\n"};/' | bc >hostfile
      mpirun=/usr/local/mpich-gm_INTEL/bin/mpirun
      mpirunops="-local -machinefile hostfile"
      export SCRATCH_DIR=`cat $TMPDIR/scratch`
      local_disc=yes
      one_local_disc=no
      nprocpernode=2
      echo '--------------- MPI_HOSTFILE ----------------'
      cat hostfile
      echo '----------- MPI_HOSTFILE - END --------------'
    elif [ "$PE" = "score" ]; then             # Using SCore?
      #set mpirunops = "-wait -F $HOME/.score/ndfile.$JOB_ID -e /tmp/scrun.$JOB_ID"
      #echo '--------------- PE_HOSTFILE ----------------'
      #cat $PE_HOSTFILE
      #echo '----------- PE_HOSTFILE - END --------------'
      mpirunops="-wait -F $PE_HOSTFILE -e $TMPDIR/scrun.$JOB_ID"
      mpirun=/opt/score/bin/scout
      #setenv SCRATCH_DIR /scratch/$JOB_ID
      export SCRATCH_DIR=$TMPDIR
      local_disc=yes
      one_local_disc=no
    fi
  else
      echo $hn >hostfile
      mpirun=/usr/local/mpich-gm_INTEL/bin/mpirun
      mpirunops="-local -machinefile hostfile"
      local_disc=no
  fi

elif ishost rasmussen ; then
  echo "Rasmussen"
  ulimit -s unlimited
  mpirun='mpirun'
  mpirunops=''
  npops=''
  ncpus=1

elif ishost hwwsr8k ; then
  echo "Hitachi in Stuttgart"
  nnodes=`echo "precision=0; ($ncpus-1)/8+1" | bc` # Hitachi-specific
  mpirun=mpiexec
  mpirunops="-p multi -N $nnodes"
  x_ops='-Fport(nmlist(2))'
  local_disc=yes
  one_local_disc=yes    # (the default anyway)
  local_binary=no
  if [ -d $SCRDIR ]; then
    export SCRATCH_DIR="$SCRDIR"
  else
    echo 'NO SUCH DIRECTORY: $SCRDIR'="<$SCRDIR> -- ABORTING"
    kill $$                     # full-featured suicide
  fi
  export SSH=rsh
  export SCP=rcp

elif ishost hwwsx5 ; then
  echo "NEC-SX5 in Stuttgart"
  mpirun=mpiexec
  mpirunops=''
  local_disc=yes
  one_local_disc=yes        # (the default anyway)
  local_binary=no
  if [ -d $SCRDIR ]; then
    export SCRATCH_DIR="$SCRDIR"
  else
    echo 'NO SUCH DIRECTORY: $SCRDIR'="<$SCRDIR> -- ABORTING"
    kill $$                     # full-featured suicide
  fi
  export SSH=rsh
  export SCP=rcp

elif ishost cosmo ; then
  mpirunops='-x NLSPATH'

else
  echo "Generic setup; hostname is <$hn>"
  [ "$mpi" = "yes" ] && echo "Use mpirun"
  mpirun=mpirun
fi

##------------------ END OF MACHINE SPECIFIC SETTINGS ---------------------
##-------------------------------------------------------------------------



# Determine appropriate options for the MPI execution environment
#   eg. number of processors
determine_mpi_processor_options

# Determine data directory (defaults to `data')
# plus a list of subdirs and procdirs
determine_datadir


# Apply the SGI namelist read fix if running IRIX
check_sgi_fix

# Propagate current pid to copy-snapshots:
export PARENT_PID=$$

# Wrap up nodelist as (scalar, colon-separateds) environment variable
# NODELIST for transport to sub-processes.
export NODELIST=`echo $nodelist | perl -ne 'print join(":",split(/\s/,$_)),"\n"'`

if [ "$debug" = "yes" ]; then show_config; fi

if match "$0" "pc_config" ; then
  if [ ! "$debug" = "yes" ]; then show_config; fi
  echo "----------------------------- pc_config ------------------------------ "
  echo " You have just run the Pencil-Code pre-execution configuration script. "
  echo "                                                                       "
  echo " This script creates some file directory structures and set a          "
  echo " collection of shell and environment variables in preperation for      "
  echo " the execution of one of the Pencil-Code binaries.                     "
  echo "                                                                       "
  echo " You will find a list of the variables that have been set above this   "
  echo " message.                                                              "
  echo "---------------------------------------------------------------------- "
fi
