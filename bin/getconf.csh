#!/bin/csh

# Name:   getconf.csh
# Author: wd (Wolfgang.Dobler@ncl.ac.uk)
# Date:   16-Dec-2001
# Description:
#  Initiate some variables related to MPI and the calling sequence. This
# is used by both start.csh and run.csh

# Are we running the MPI version?
set mpi = `egrep -c '^[ 	]*MPICOMM[ 	]*=[ 	]*mpicomm' src/Makefile.local`

if ($mpi) then
  echo "Running under MPI"
  # Compaq has `dmpirun' instead of `mpirun'; some mpiruns have specail path
  set hn = `hostname`
  if ($hn =~ mhd*.st-and.ac.uk) then
    set mpirun = dmpirun
  else if ($hn =~ *.kis.uni-freiburg.de) then
    set mpirun = /opt/local/mpich/bin/mpirun
  else if (($hn =~ cincinnatus*) || ($hn =~ owen*)) then
    set mpirun = /usr/lib/lam/bin/mpirun
  else if ($hn =~ nq*) then
    set mpirun = /nq1/home/jes/mpich-1.2.1/bin/mpirun
  else
    set mpirun = mpirun
  endif
  # Some mpiruns need special options
  if (`domainname` == "aegaeis") then
    set mpirunops = '-machinefile ~/mpiconf/mpihosts-martins'
  else
    set mpirunops = ''
  endif
  # Determine number of CPUS
  set ncpus = `perl -ne '$_ =~ /^\s*integer\b[^\\!]*ncpus\s*=\s*([0-9]*)/i && print $1' src/cparam.local`
  echo $ncpus CPUs
  # Number of processors
  set npops = "-np $ncpus"
else # no MPI
  echo "Non-MPI version"
  set mpirun = ''
  set mpirunops = ''
  set npops = ''
  set ncpus = 1
endif
  
# End of file getconf.csh
