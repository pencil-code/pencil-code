#!/bin/csh
#
# run.csh -- driver for time stepping
#
##PBS -S /bin/csh -W group_list=UK07001 -q UK07001
##PBS -l ncpus=64,mem=32gb,walltime=300:00:00
#PBS -S /bin/csh
#PBS -q workq
##PBS -l ncpus=4,mem=2mb,walltime=0:00:10
##PBS -l ncpus=4
##PBS -l nodes=4,walltime=0:00:10
#PBS -l nodes=64
##PBS -q p-long
##PBS -l nodes=nq1+nq2+nq3+nq4

#setenv PGHPF_HOST -file=$PBS_NODEFILE
#if ($?PBS_NODEFILE) then
#  setenv PGHPF_HOST -file=$PBS_NODEFILE
#endif
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
echo $mpirun $mpirunops $npops src/run.x >>run_command.log
time $mpirun $mpirunops $npops src/run.x
#
date

# cut & paste for job submission on the mhd machine
# bsub -n  4 -q 4cpu12h -o run.`timestr` run.csh
# bsub -n  8 -q 8cpu12h -o run.`timestr` run.csh
# bsub -n 16 -q 16cpu8h -o run.`timestr` run.csh
# bsub -n  8 -q 8cpu12h -o run.log -w 'exit(123456)' mpijob dmpirun src/run.x
# bsub -n 16 -q 16cpu8h mpijob dmpirun src/run.x
