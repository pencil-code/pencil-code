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
    echo "lamndodes:"
    lamnodes
    set mpirun = /usr/bin/mpirun
    set mpirunops = "-O -c2c -s n0"
    if ($local_disc) then
       setenv SCRATCH_DIR "/var/tmp"
       set start_x = $SCRATCH_DIR/start.x
       set run_x = $SCRATCH_DIR/run.x
    endif

  else if ($hn =~ s[0-9]*p[0-9]*) then
    echo "Use options for the Horseshoe cluster"
    set nodelist = `cat $PBS_NODEFILE`
    cat $PBS_NODEFILE > lamhosts
    lamboot -v lamhosts
    echo "lamndodes:"
    lamnodes
    set mpirun = mpirun
    #set mpirunops = "-O -s n0 N -lamd"
    set mpirunops = "-O -c2c -s n0 N -v"
    setenv SCRATCH_DIR /scratch
    set start_x = $SCRATCH_DIR/start.x
    set run_x = $SCRATCH_DIR/run.x
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
exit

# End of file getconf.csh
