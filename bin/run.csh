#!/bin/csh
#
# run.csh -- driver for time stepping
#
##PBS -S /bin/csh -W group_list=UK06005 -q UK06005
##PBS -l ncpus=32,mem=16gb,walltime=200:00:00
#PBS -l ncpus=4
#PBS -q p-long
#PBS -l nodes=nq0+nq4+nq2+nq3

#setenv PGHPF_HOST -file=$PBS_NODEFILE
if ($?PBS_NODEFILE) then
  setenv PGHPF_HOST -file=$PBS_NODEFILE
endif
#setenv MPI_HOST nq1,nq2,nq3,nq4
#echo $PBS_NODEFILE

if ($?PBS_O_WORKDIR) then
  cd $PBS_O_WORKDIR
endif

# Determine whether this is MPI, how many CPUS etc.
source getconf.csh

# Clean up control and data files
rm -f STOP RELOAD fort.20

# Run run.x
date
#
echo "$mpirun $mpirunops $npops src/run.x"
echo $mpirun $mpirunops $npops src/run.x >>run.log
time $mpirun $mpirunops $npops src/run.x
#
date

# bsub -n  4 -q 4cpu12h -o run.`timestr` -e run.`timestr` run.csh
# bsub -n  8 -q 8cpu12h -o run.`timestr` -e run.`timestr` run.csh
# bsub -n 16 -q 16cpu8h -o run.`timestr` -e run.`timestr` run.csh
