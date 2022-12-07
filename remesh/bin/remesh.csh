#!/bin/csh
#
#  Run remesh.x -- to remesh
#
#PBS -S /bin/csh

# examples of the calling sequence:
# qsub -l nodes=1,walltime=5:00:00 -q workq remesh.csh
# qsub -l nodes=1,walltime=5:00:00 -q giga2 remesh.csh
# qsub -l nodes=1,walltime=5:00:00 -q giga remesh.csh
# qsub -l ncpus=1,cput=1:00:00 -q parallel remesh.csh
# aprun -n 1 remesh/remesh.x

if ($?PBS_O_WORKDIR) then
  cd $PBS_O_WORKDIR
endif

date

# write the current directory to new directory (for information)
set targetdir = `cat remesh.in`
echo $cwd >>! $targetdir/data/remeshed_from.dir
date >>! $targetdir/data/remeshed_from.dir
cp data/param.nml $targetdir/data/

# generate directory tree
set ncpus = `perl -ne '$_ =~ /^\s*integer\b[^\\!]*ncpus\s*=\s*([0-9]*)/i && print $1' $targetdir/src/cparam.local`
echo prepare directory tree for $ncpus processors
cat $targetdir/src/cparam.local

# Create list of subdirectories
set subdirs = (`printf "%s%s%s\n" "for(i=0;i<$ncpus;i++){" '"data/proc";' 'i; }' | bc` 'data/allprocs' 'data/averages')
foreach dir ($subdirs)
  # Make sure a sufficient number of subdirectories exist
  if (! -e $targetdir/$dir) then
    mkdir $targetdir/$dir
  endif
  #cp data/proc0/seed.dat $targetdir/$dir
end

# Run
echo 'Starting remesh/remesh.x'
srun remesh/remesh.x

echo ""
date
