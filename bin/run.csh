#!/bin/csh
#                       run.csh
#                      ---------
#   Run src/run.x (timestepping for src/run.x).
#   Run parameters are set in run.in.
#
# Run this script with csh:
#PBS -S /bin/csh
#$ -S /bin/csh
#@$-s /bin/csh
#
# Join stderr and stout:
#$ -j y -o run.log
#@$-eo
#
# Work in submit directory (SGE):
#$ -cwd

# Work in submit directory (PBS):
if ($?PBS_O_WORKDIR) then
  cd $PBS_O_WORKDIR
endif

# Work in submit directory (SUPER-UX's nqs):
if ($?QSUB_WORKDIR) then
  cd $QSUB_WORKDIR
endif

# ---------------------------------------------------------------------- #

# Common setup for start.csh, run.csh, start_run.csh:
# Determine whether this is MPI, how many CPUS etc.
source getconf.csh

# Prevent code from running twice (and removing files by accident)
if (! -e "NEVERLOCK") touch LOCK

#
#  If necessary, distribute var.dat from the server to the various nodes
#
if ($local_disc) then
  if ($one_local_disc) then	# one common local disc
    foreach node ($nodelist)
      foreach d (`cd $datadir; ls -d proc* allprocs`)
	$SCP $datadir/$d/var.dat ${node}:$SCRATCH_DIR/$d/
      end
    end
  else # one local disc per MPI process (Horseshoe, etc);
       # still doesn't cover Copson
    set i=0
    foreach node ($nodelist)
      echo "i = $i"
      $SCP $datadir/proc$i/var.dat ${node}:$SCRATCH_DIR/proc$i/
      set i=`expr $i + 1`
    end
  endif
endif
# Clean up control and data files
rm -f STOP RELOAD fort.20

# On machines with local scratch directory, initialize automatic
# background copying of snapshots back to the data directory.
# Also, if necessary copy executable to $SCRATCH_DIR of master node
# and start top command on all procs
if ($local_disc) then
  echo "Use local scratch disk"
  copy-snapshots -v >& copy-snapshots.log &
  remote-top >& remote-top.log &
endif
if ($local_binary) then
  echo "ls src/run.x $SCRATCH_DIR before copying:"
  ls -lt src/run.x $SCRATCH_DIR
  cp src/run.x $SCRATCH_DIR
  echo "ls src/run.x $SCRATCH_DIR after copying:"
  ls -lt src/run.x $SCRATCH_DIR
endif

# Run run.x
date
echo "$mpirun $mpirunops $npops $run_x $x_ops"
echo $mpirun $mpirunops $npops $run_x $x_ops >! run_command.log
time $mpirun $mpirunops $npops $run_x $x_ops
set run_status=$status		# save for exit
date

# Copy var.dat back from local scratch to data directory
if ($local_disc) then
  echo "Copying final var.dat back from local scratch disk"
  copy-snapshots -v var.dat
  echo "done, will now killall copy-snapshots"
  killall copy-snapshots
endif

# Shut down lam if we have started it
if ($booted_lam) lamhalt

# remove LOCK file
if (-e "LOCK") rm -f LOCK

exit $run_status		# propagate status of mpirun

# cut & paste for job submission on the mhd machine
# bsub -n  4 -q 4cpu12h -o run.`timestr` -e run.`timestr` run.csh
# bsub -n  8 -q 8cpu12h mpijob dmpirun src/run.x
# bsub -n 16 -q 16cpu8h mpijob dmpirun src/run.x
# bsub -n  8 -q 8cpu12h -o run.log -w 'exit(123456)' mpijob dmpirun src/run.x

# qsub -l ncpus=64,mem=32gb,walltime=500:00:00 -W group_list=UK07001 -q UK07001 run.csh
# qsub -l nodes=4:ppn=1,mem=500mb,cput=24:00:00 -q p-long run.csh
# qsub -l ncpus=16,mem=1gb,cput=400:00:00 -q parallel run.csh
# qsub -l nodes=128,walltime=10:00:00 -q workq run.csh
# eval `env-setup lam`; qsub -v PATH -pe lam 8 -j y -o run.log run.csh
