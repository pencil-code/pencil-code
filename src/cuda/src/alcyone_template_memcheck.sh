#!/bin/bash
#SBATCH -J astar_test 
#SBATCH -N 1 
#SBATCH -n 1
#SBATCH -t 00:15:00
#SBATCH -o output_astaroth_test_%j.out
#SBATCH -e output_astaroth_test_%j.err
#SBATCH --partition="2G_gpu3,2G_gpu6"
#SBATCH --gres=gpu:1

module load cuda/5.5.22

#-------------------------------------------------------
#
# First check that we have enough disk space on /tmp.
# If not the script exits with status 64. The status 
# is also shown in the email message sent on job 
# completion.
#
# The limit (in %) is given in shell variable tmplimit 
# (note: it is an integer).
#
#-------------------------------------------------------

tmplimit=80

space=`df /tmp | awk '!/Filesystem/ && /%/ {sub("%","",$5); print $5+0;}'`
echo
echo "Disk space usage on /tmp on computation node:" $space "%"

if test 0$space -ge $tmplimit
then
    echo
    echo "qsub script check:"
    echo "File space on tmp on node" `hostname` "dangerously low. Aborting job."
    echo "Please ask an administrator to clean up /tmp on that node."
    echo
    df /tmp
    echo
    exit 64
fi



#-------------------------------------------------------
# Create the computation work directory.
#
# Naming pattern is '/tmp/<USERNAME>/__<CWD>_<N>/', 
# where <CWD> is the trailing part of the current 
#             directory and 
#       <ID>  is the process ID of the current shell 
#             (just to get a unique name).
#-------------------------------------------------------

sd=`pwd`
prefix=/tmp/$USER
if [ ! -d $prefix ]; then
   mkdir $prefix
fi
cd=${sd/$HOME/}
cd=${cd//\//__}
cd=$prefix/${cd}_${SLURM_JOB_ID}
mkdir $cd

#-------------------------------------------------------
# Copy stuff to running directory
#-------------------------------------------------------

rsync -avu $sd/ $cd/

echo "Date:                   " `date`
echo "Submission directory:   " $sd
echo "Running directory:      " `pwd`
echo "Host:                   " `hostname`
echo "JOB_ID:                 " $SLURM_JOB_ID

cd $cd

#-------------------------------------------------------
# Run the program and copy the results
#-------------------------------------------------------

cuda-memcheck ./runme 
rsync -avu $cd/ $sd/

#-------------------------------------------------------
# Cleanup, i.e. remove the computation work directory.
# This is important. Otherwise /tmp would gradually 
# fill up.
#-------------------------------------------------------

cd $sd/
rm -rf $cd
