#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pencil as pc
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import cm
import glob
import fileinput
import sys
#TODO rvid_box plots
#TODO sim objects and hdf5
""" Tools to copy old binary format simulations into new hdf5 format for protability and continuation
    Remeshing of old simulations in binary and/or hdf5 format shall follow
"""


# In[ ]:


pencil_home = os.getenv("PENCIL_HOME")


# In[ ]:


os.chdir(pencil_home+'/python/tutorials/mps')


# In[ ]:


#identify the top directory for sim_path and list dir contents
sim_path = os.getcwd()
sim_list = glob.glob('*')
#par_sim = os.path.join(sim_path,'restart')
print(sim_path)
print(sim_list)


# In[ ]:


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
    sim_ = os.path.join(sim_path,sim)
    h5_already=os.path.exists(os.path.join(sim_,'data/grid.h5'))
    if pc.sim.is_sim_dir(sim_) and not h5_already:
        #if simulation is to be processed move to directory
        os.chdir(sim_)
        cmd = 'du -csh '+sim_
        os.system(cmd)
        #edit Makefile.local from io_dist to io_hdf5 pre-compiling
        for i, line in enumerate(fileinput.input('src/Makefile.local',
                                                                    inplace=1)):
            sys.stdout.write(line.replace('io_dist','io_hdf5'))
        print('cleaning '+sim_)

        #create a bash script to compile pc with src/read_all_videoslices.x
        #may be commented out if executable already present (a) or data/slice_* complete (b)
        #amend/comment the lines as required for your fortran compiler and python envcomment 
        f = open('pc_compiler','w')
        f.write('#!/bin/bash\n')
        f.write('\n')
        f.write('module purge\n')
        f.write('module load intel/16.0.0 intelmpi hdf5-par\n')
        f.write('module list\n')
        f.write('\n')
        #minimal next 4 lines required to add executables depending on compiler
        f.write('pc_setupsrc\n')
        f.write('make distclean\n')
        f.write('pc_build -f Intel_MPI.conf\n')
        f.write('pc_build -f Intel_MPI.conf -t ALL\n')
        f.write('\n')
        #following lines required to restore python environment depending on your system
        f.write('module purge\n')
        f.write('module load python-env/3.5.3\n')
        f.write('module load hdf5-par\n')
        f.write('export PYTHONPATH=$USERAPPL/Pencilcodepython3/lib/python3.5/site-packages\n')
        f.write('export PYTHONPATH="$PYTHONPATH:$PENCIL_HOME/python"\n')
        f.write('module list\n')
        f.close()
        cmd = 'bash pc_compiler'
        os.system(cmd)
        #(a) files above to provide src/read_all_videoslices
        #comment cmd call above if executables already present
        if not os.path.exists('data/slice_uu1.xy'):
            lvids = False
            for iproc in range(0, ncpus):
                if os.path.exists('data/proc{}/slice_uu1.xy'.format(iproc)):
                    lvids=True
                    cmd = './src/read_all_videofiles.x'
                    os.system(cmd)
                    break
        #(b) files above to construct data/slices from proc tree binary
        #slices files.
        #may still want to call above cmd if slices need to be updated

        #copy scripts from sim_path which will be edited for each sim
        #typically the data can be converted in the same sim directory
        #if confident, to save disk space redundant files can be removed during the process
        #otherwise they can be retained until the conversion is complete and verified
        #here we create a new sim 'ism_binary2h5' for test purposes

        #script to call fort2h5.py if no mpi and var processing by proc not needed
        #here we shall use only serial python - edit as required for mpi
        cmd = 'cp '+os.path.join(sim_path,'loc-fort2h5.py')              +' '+os.path.join(sim_path,sim_,'loc-fort2h5.py')
        os.system(cmd)
        #script to call fort2h5.py with mpi or var processing by proc needed
        cmd = 'cp '+os.path.join(sim_path,'par-fort2h5.py')              +' '+os.path.join(sim_path,sim_,'par-fort2h5.py')
        #batch script to call fort2h5.py with mpi or var processing by proc needed
        #a python environment compiled for hdf5-parallel must be available to the HPC nodes 
        os.system(cmd)
        cmd = 'cp '+os.path.join(sim_path,'job_python')              +' '+os.path.join(sim_path,sim_,'job_python')
        os.system(cmd)

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
        if not os.path.exists('data/proc0/var.dat'):
            lvars = False
            for i, line in enumerate(fileinput.input('loc-fort2h5.py',inplace=1)):
                sys.stdout.write(line.replace('lvars=True','lvars=False'))
        if not os.path.exists('data/xyaverages.dat'):
            for i, line in enumerate(fileinput.input('loc-fort2h5.py',inplace=1)):
                sys.stdout.write(line.replace('laver=True','laver=False'))
        if not lvids:
            for i, line in enumerate(fileinput.input('loc-fort2h5.py',inplace=1)):
                sys.stdout.write(line.replace('lvids=True','lvids=False'))

        #edit files to offer new directory path '*2h5' - omit if using same sim directory
        for i, line in enumerate(fileinput.input('loc-fort2h5.py',inplace=1)):
            sys.stdout.write(line.replace('sim2h5(','sim2h5(\n             newdir="'+sim_+'2h5",'))
        for i, line in enumerate(fileinput.input('par-fort2h5.py',inplace=1)):
            sys.stdout.write(line.replace('sim2h5(','sim2h5(\n             newdir="'+sim_+'2h5",'))
        #(c) if the old 'data/param.nml' file leads to reading errors for the old data
        #cmd = 'cp '+os.path.join(sim_path,'ism','data','param.nml')+\
        #      os.path.join(sim_,'data','param.nml')
        #os.system(cmd)

        if ncpus > 24:
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
                cmd = 'sbatch job_python'
                os.system(cmd)
            else:
                #apply mpi without snap_by_proc for var files
                for i, line in enumerate(fileinput.input(
                                                       'job_python',inplace=1)):
                    sys.stdout.write(line.replace('par-fort2h5.py',
                                                                  'loc-fort2h5.py'))
                cmd = 'sbatch job_python'
                os.system(cmd)
        else:
            #no mpi command line call
            cmd = 'python loc-fort2h5.py'
            os.system(cmd)


# In[ ]:


""" The following section details remeshing a simulation from an existing.
    This will shortly be revised to use the pencil modules and apply either 
    fortran binary or hdf5 formats. We shall also demonstrate the capacity
    to use the pencil.sim modules to setup and array of simulations, compile
    and run them from a single script. (To follow)
"""


# In[ ]:


#using Plotly likely requires to install plotly in conda
#$ conda install -c plotly plotly-orca psutil requests


# In[ ]:


#identify old sim to remesh
wrkdir=os.path.join(sim_path,'ism_binary')


# In[ ]:


#create new sim from existing parameters
os.chdir(wrkdir)
meshdir = wrkdir+"remesh"
cmd = "pc_newrun "+cmd = "pc_newrun "+meshdir
os.system(cmd)


# In[ ]:


#remesh binary sim using old pencil script TO BE REVISED IN PENCIL
from pencil_old.files import remesh


# In[ ]:


#What are the existing grid dimensions
os.chdir(wrkdir)
old_grid=pc.read.grid(quiet=True)
old_dim=pc.read.dim()
old_param=pc.read.param(quiet=True)
old_index=pc.read.index()


# In[ ]:


#arbitrary regridding permissible except that each direction must be divisible by new processor layout
print(int(old_dim.nxgrid*0.8), int(old_dim.nygrid*0.8), int(old_dim.nzgrid*0.8))
print(old_dim.nprocx,old_dim.nprocy,old_dim.nprocz)
print(old_param.xyz0,old_param.xyz1)


# In[ ]:





# In[ ]:





# In[ ]:


#Before this call edit cparam.local in 
#default var.dat copied and written
arrs=[]
for key in old_index.__dict__.keys():
    if 'keys' not in key:
        arrs.append(key)
print(arrs)
fnew=remesh.interp_var(target_path=newdir, source_path=olddir, arrs=arrs)


# In[ ]:




