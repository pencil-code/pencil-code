#!/bin/csh

# Name:   getconf.csh
# Author: wd (Wolfgang.Dobler@ncl.ac.uk)
# Date:   16-Dec-2001
# Description:
#  Initiate some variables related to MPI and the calling sequence. This
# is used by both start.csh and run.csh

# Are we running the MPI version?
set mpi = `egrep -c '^[ 	]*MPICOMM[ 	]*=[ 	]*mpicomm' src/Makefile.local`

# location of executables; can be overwritten below
set start_x = "src/start.x"
set run_x = "src/run.x"

# settings for machines with local data disks
set local_disc = 0  #use global file system by default
setenv SCRATCH_DIR /scratch
setenv SSH ssh
setenv SCP scp

# choose machine specific settings
echo `uname -a`
set hn = `hostname`
if ($mpi) then
  echo "Running under MPI"
  set mpirunops = ''
  if ($hn =~ mhd*.st-and.ac.uk) then
    echo "St Andrews machine"
    set mpirun = "dmpirun"

  else if ($hn =~ *.kis.uni-freiburg.de) then
    set mpirun = /opt/local/mpich/bin/mpirun

  else if (($hn =~ cincinnatus*) || ($hn =~ owen*) || ($hn =~ master)) then
    set mpirun = /usr/lib/lam/bin/mpirun
    set mpirunops = "-c2c"
    set mpirunops = "-c2c -O"
#    set mpirunops = " c0-7"
#    set mpirunops = "-c2c c8-13"

  else if ($hn =~ nq*) then
    echo "Use options for the Nordita cluster"
    if ($?PBS_NODEFILE ) then
      set nodelist = `cat $PBS_NODEFILE`
      cat $PBS_NODEFILE > lamhosts
      set local_disc = 1
    endif
    lamboot -v lamhosts
    echo "lamnodes:"
    lamnodes
    set mpirun = /usr/bin/mpirun
    set mpirunops = "-O -c2c -s n0"
    if ($local_disc) then
       setenv SCRATCH_DIR "/var/tmp"
       set start_x = $SCRATCH_DIR/start.x
       set run_x = $SCRATCH_DIR/run.x
    endif

  else if ($hn =~ s[0-9]*p[0-9]*) then
    if ($?RUNNINGMPICH) then
      echo "Running using MPICH"
      set mpirunops = "-machinefile $PBS_NODEFILE"
      set mpirun = /usr/local/lib/MPICH/bin/mpirun
      set start_x=$PBS_O_WORKDIR/src/start.x
      set run_x=$PBS_O_WORKDIR/src/run.x
    else
      echo "Use LAM-MPI options for the Horseshoe cluster"
      set nodelist = `cat $PBS_NODEFILE`
      cat $PBS_NODEFILE > lamhosts
#echo $nodelist > lamhosts.before
#shift nodelist
#cat $nodelist > lamhosts
      lamboot -v lamhosts
      echo "lamnodes:"
      lamnodes
      set mpirunops = "-O -c2c -s n0 N -v" #(direct; should be faster)
      # set mpirunops = "-O -s n0 N -lamd" #(with debug options on; is slower)
      set mpirun = mpirun
      set start_x = $SCRATCH_DIR/start.x
      set run_x = $SCRATCH_DIR/run.x
    endif

    setenv SCRATCH_DIR /scratch
    set local_disc = 1
    setenv SSH rsh
    setenv SCP rcp

  else
    echo "Use mpirun as the default option"
    set mpirun = mpirun
  endif

  # Some mpiruns need special options
  if (`domainname` == "aegaeis") then
    set mpirunops = '-machinefile ~/mpiconf/mpihosts-martins'
  endif

  # Determine number of CPUS
  set ncpus = `perl -ne '$_ =~ /^\s*integer\b[^\\!]*ncpus\s*=\s*([0-9]*)/i && print $1' src/cparam.local`
  echo $ncpus CPUs
  set npops = "-np $ncpus"
else # no MPI
  echo "Non-MPI version"
  set mpirun = ''
  set mpirunops = ''
  set npops = ''
  set ncpus = 1
endif

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
if ($local_disc) then
  echo $SCRATCH_DIR >$datadir/directory_snap
else
  if (-f $datadir/directory_snap) rm $datadir/directory_snap
endif

exit

# If we are using local disk, the directory name is written into
# data/directory_snap. This is read by the code during run time
if ($local_disc) echo $SCRATCH_DIR >$datadir/directory_snap

# End of file getconf.csh
