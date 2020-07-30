#!/usr/bin/env python
# coding: utf-8

# In[8]:


import pencil as pc
import numpy as np
import os
from os.path import join, exists
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import cm
import glob
import fileinput
import sys
import subprocess as sub
#TODO rvid_box plots
#TODO sim objects and hdf5
""" Tools to copy old binary format simulations into new hdf5 format for protability and continuation
    Remeshing of old simulations in binary and/or hdf5 format shall follow
"""

#pencil_home = os.getenv("PENCIL_HOME")
#home = os.getenv("SCRATCH")
#os.chdir('/scratch/project_2001062')
wrkdir=os.getcwd()
print(wrkdir)

#os.chdir(pencil_home+'/python/tutorials/mps')
#os.chdir(join(wrkdir,'ism_binary'))
#identify the top directory for sim_path and list dir contents
sim_path = os.getcwd()
sim_list = np.array(glob.glob('*'))
#par_sim = join(sim_path,'restart')
print(sim_path)
print(sim_list)

""" An old simulation which has been executed with fortran binary
    files distributed across a processor tree can be converted to
    hdf5, redundant data removed and the run continued in hdf5.
    This script loops through all items in the sim_path and 
    identifies which objects are existing simulation objects.
    If a simulation directory contains the file 'data/grid.h5'
    it is assumed already to be in hdf5 format so is skipped.
    For some very old simulations errors in loading data may occur
    for obsolete 'data/param.nml' files. If this occurs a dummy
    run with the necessary setup can be compiled started and the
    new param.nml replace the existing file. This can be automated
    see comment (c) below.
"""
#process dir subdirectories if simulation and if not already hdf5
for sim in sim_list:
    print(sim)
    sim_ = join(sim_path,sim)
    h5_already=exists(join(sim_,'data/grid.h5'))
    if pc.sim.is_sim_dir(sim_) and not h5_already:
        #if simulation is to be processed move to directory
        os.chdir(sim_)
        cmd = 'du -csh '+sim_
        os.system(cmd)
        #edit Makefile.local from io_dist to io_hdf5 pre-compiling
        for i, line in enumerate(fileinput.input('src/Makefile.local',
                                                                    inplace=1)):
            sys.stdout.write(line.replace('io_dist','io_hdf5'))
        #print('cleaning '+sim_)

        ##uncomment next two lines if video slices do not exist
        #cmd = 'touch src/read_all_videofiles.x'
        #os.system(cmd)
        #create a bash script to compile pc with src/read_all_videoslices.x
        #may be commented out if executable already present (a) or data/slice_* complete (b)
        #amend/comment the lines as required for your fortran compiler and python envcomment
        #if not exists(join(sim_,'src/read_all_videofiles.x')):
        #    f = open('pc_compiler','w')
        #    f.write('#!/bin/bash\n')
        #    f.write('\n')
        #    f.write('module purge\n')
        #    f.write('module load intel\n module load intel-mpi\n module load hdf5/1.10.4-mpi\n')
        #    f.write('module list\n')
        #    f.write('\n')
        #    #minimal next 4 lines required to add executables depending on compiler
        #    f.write('pc_setupsrc\n')
        #    #f.write('make distclean\n')
        #    f.write('pc_build -f Intel_MPI.conf\n')
        #    f.write('pc_build -f Intel_MPI.conf -t ALL\n')
        #    f.write('\n')
        #    #following lines required to restore python environment depending on your system
        #    f.write('module purge\n')
        #    f.write('module load python-data\n')
        #    #f.write('module load hdf5/1.10.4-mpi\n')
        #    #f.write('export PYTHONPATH=$USERAPPL/Pencilcodepython3/lib/python3.5/site-packages\n')
        #    #f.write('export PYTHONPATH="$PYTHONPATH:$PENCIL_HOME/python"\n')
        #    f.write('module list\n')
        #    f.close()
        #    cmd = 'bash pc_compiler'
        #    os.system(cmd)
        #(a) files above to provide src/read_all_videoslices
        #slices files.
        #may still want to call above cmd if slices need to be updated

        #copy scripts from sim_path which will be edited for each sim
        #typically the data can be converted in the same sim directory
        #if confident, to save disk space redundant files can be removed during the process
        #otherwise they can be retained until the conversion is complete and verified
        #here we create a new sim 'ism_binary2h5' for test purposes

        #script to call fort2h5.py if no mpi and var processing by proc not needed
        #here we shall use only serial python - edit as required for mpi
        cmd = 'cp '+join(sim_path,'loc-fort2h5.py')+' '+join(sim_path,sim_,'loc-fort2h5.py')
        os.system(cmd)
        #script to call fort2h5.py with mpi or var processing by proc needed
        cmd = 'cp '+join(sim_path,'par-fort2h5.py')+' '+join(sim_path,sim_,'par-fort2h5.py')
        #batch script to call fort2h5.py with mpi or var processing by proc needed
        #a python environment compiled for hdf5-parallel must be available to the HPC nodes 
        os.system(cmd)
        cmd = 'cp '+join(sim_path,'puhti_python')+' '+join(sim_path,sim_,'puhti_python')
        os.system(cmd)
        for i, line in enumerate(fileinput.input(
                                               'puhti_python',inplace=1)):
            sys.stdout.write(line.replace('parallelh5py',
                                                          'py'+sim))

        #read cparam.local to determine number ncpus for proc tree of binary files
        for line in open('src/cparam.local','r').readlines():
            if 'ncpus' in line:
                ncpus = int(str.split(str.split(str.split(
                            line,',')[1],'::')[1],'=')[1])
                break
        print(sim_,'ncpus',ncpus)

        #having copied the files to the sim directory we shall edit the call as files are
        #available or in case the conversion is partially complete
        #by default var files, videos slices and averages are included
        #if large 2D averages set also laver2D=True, laver sufficient for smaller 2D datasets
        lvars, lvids, laver = True, True, True
        #(b) lines below to construct data/slices from proc tree binary
        if not exists('data/slice_uu1.xy'):
            lvids = False
            for iproc in range(0, ncpus):
                if exists('data/proc{}/slice_uu1.xy'.format(iproc)):
                    lvids=True
                    cmd = './src/read_all_videofiles.x'
                    os.system(cmd)
                    break
        if not exists('data/proc0/var.dat'):
            lvars = False
            for i, line in enumerate(fileinput.input('loc-fort2h5.py',inplace=1)):
                sys.stdout.write(line.replace('lvars=True','lvars=False'))
        if not exists('data/xyaverages.dat'):
            for i, line in enumerate(fileinput.input('loc-fort2h5.py',inplace=1)):
                sys.stdout.write(line.replace('laver=True','laver=False'))
        if not lvids:
            for i, line in enumerate(fileinput.input('loc-fort2h5.py',inplace=1)):
                sys.stdout.write(line.replace('lvids=True','lvids=False'))

        #edit files to offer new directory path '*2h5' - omit if using same sim directory
        for i, line in enumerate(fileinput.input('loc-fort2h5.py',inplace=1)):
            sys.stdout.write(line.replace('sim2h5(','sim2h5(\n             olddir="'+sim_+'",'))
        for i, line in enumerate(fileinput.input('par-fort2h5.py',inplace=1)):
            sys.stdout.write(line.replace('sim2h5(','sim2h5(\n             olddir="'+sim_+'",'))
        #(c) if the old 'data/param.nml' file leads to reading errors for the old data
        #cmd = 'cp '+join(sim_path,'ism','data','param.nml')+\
        #      join(sim_,'data','param.nml')
        #os.system(cmd)

        if ncpus > 1:
            #use mpi with batch script submission
            if lvars:
                #var files handled by processor so amend the append loc-fort2h5.py
                for i, line in enumerate(fileinput.input(
                                                   'loc-fort2h5.py',inplace=1)):
                    sys.stdout.write(line.replace('lvars=True',
                                                             'lvars=False')) 
                f = open('par-fort2h5.py','a')
                f.write(open('loc-fort2h5.py').read())
                f.close()
                cmd = 'sbatch puhti_python'
                process = sub.Popen(cmd.split(),stdout=sub.PIPE)
                output, error = process.communicate()
                print(cmd,output,error)
            else:
                #apply mpi without snap_by_proc for var files
                for i, line in enumerate(fileinput.input(
                                                       'puhti_python',inplace=1)):
                    sys.stdout.write(line.replace('par-fort2h5.py',
                                                                  'loc-fort2h5.py'))
                cmd = 'sbatch puhti_python'
                process = sub.Popen(cmd.split(),stdout=sub.PIPE)
                output, error = process.communicate()
                print(cmd,output,error)
        else:
            #no mpi command line call
            cmd = 'python loc-fort2h5.py'
            process = sub.Popen(cmd.split(),stdout=sub.PIPE)
            output, error = process.communicate()
            print(cmd,output,error)

