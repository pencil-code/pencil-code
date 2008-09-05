#PBS -S /bin/csh
# CVS: $Id$

# (I'm not sure this script is currently being used!?)

# Preliminaries for batch jobs
if ($?PBS_O_WORKDIR) then
  cd $PBS_O_WORKDIR
endif

# copy to your working directory and modify correspondingly

#copy-proc-to-proc var.dat ../rot512_Om0b
#copy-VAR-to-var 2 .

# qsub -l nodes=1,walltime=1:00:00 -q workq restart.csh
