#!/bin/csh
#
# run.csh -- driver for time stepping
#
#PBS -S /bin/csh

# cut & paste for job submission for PBS
# qsub -l ncpus=64,mem=32gb,walltime=500:00:00 -W group_list=UK07001 -q UK07001 run.csh
# qsub -l nodes=nq1+nq2+nq3+nq4,mem=1gb,cput=24:00:00 -q p-long run.csh
# qsub -l nodes=4 -q p-long run.csh
# qsub -l nodes=nq1+nq2+nq3+nq4 -q p-long run.csh
# qsub -l ncpus=4,nodes=nq1+nq2 -q p-long run.csh
# qsub -l ncpus=4,mem=1gb,cput=100:00:00 -q parallel run.csh
# qsub -l ncpus=16,mem=1gb,cput=400:00:00 -q parallel run.csh
# qsub -l nodes=8,walltime=20:00:00 -q workq run.csh
# qsub -l nodes=128,walltime=10:00:00 -q workq run.csh
# qsub -l nodes=128 -q workq run.csh

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
# bsub -n  4 -q 4cpu12h -o run.`timestr` -e run.`timestr` run.csh
# bsub -n  8 -q 8cpu12h mpijob dmpirun src/run.x
# bsub -n 16 -q 16cpu8h mpijob dmpirun src/run.x
# bsub -n  8 -q 8cpu12h -o run.log -w 'exit(123456)' mpijob dmpirun src/run.x
