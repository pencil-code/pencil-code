#!/usr/bin/env python
# coding: utf-8
""" example of how shock tube analytic solution might be compared to the
simulation results. Assume the directory sod_tests contains one or more
subdirectories of shocktube test simulations with various parameters.

The first call takes the parameters of each simulation and constructs an hdf5
file containing the time and saved states of the numeric and analytic solutions.

The second call produces a set of comparison plots for each variable at each
snapshot and saves the figures in the simulation directory. From the hdf5
files plots could also be made combining multiple parameter compared to the
analytic solution at a given time - not shown here.
"""

import pencil as pc
from pencil.math import natural_sort
from pencil.calc import calc_shocktube
import os
from os.path import join, exists
import sys
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import cm
from matplotlib import artist
import glob
import string

_sod_tests = '/scratch/project_2001062/sod_tests'
if not os.path.isdir(_sod_tests):
    print("This Python file requires the directory {} to exist".format(_sod_tests))
sys.exit(0)

os.chdir(_sod_tests) # root dir of multiple sims
wrkdir = os.getcwd() # variable used for chdir
sims_list = natural_sort(glob.glob('*')) # ordered list of sub dirs
sims = np.array(sims_list) # convert list to np array

indx = [0,1,2] # list of the subset of sims which will be used, can be revised for each action

# create objects for each sim which will be used
for sim in sims:
    print('processing parameters for'+sim)
    os.chdir(join(wrkdir,sim)) # cd to sim dir
    globals()[sim+'par'] = pc.read.param()
    globals()[sim+'gd'] = pc.read.grid(quiet=True,trim=True)
    globals()[sim+'dim'] = pc.read.dim()


for sim in sims[:]:
    print('processing varfiles for'+sim)
    os.chdir(join(wrkdir,sim))
    if exists('data/allprocs/var.h5'):
        os.chdir('data/allprocs')
    else:
        os.chdir('data/proc0')
    varfiles = natural_sort(glob.glob('VAR*'))
    os.chdir(join(wrkdir,sim))
    if par.lmagnetic:
        magic = ['tt', 'pp', 'bb']
    else:
        magic = ['tt', 'pp']
    with h5py.File(join(wrkdir,sim,'data','shocktest.h5'), 'w') as hf:
        for varfile in varfiles:
            var = pc.read.var(varfile,quiet=True,magic=magic,trimall=True)
            # calc pp, uu, rho of analytic solution
            pp, uu, rho = calc_shocktube(
                                         globals()[sim+'gd'].x,
                                         var.t,
                                         par=globals()[sim+'par'],
                                         lreference=False,
                                         DEBUG=False,
                                         lplot=False)
            # calc and append to list ee
            ee = pp/rho/(globals()[sim+'par'].gamma-1)
            grp=hf.create_group(str.strip(varfile,'.h5'))
            grp.create_dataset('pp-sol',data=pp)
            grp.create_dataset('uu-sol',data=uu)
            grp.create_dataset('ee-sol',data=ee)
            grp.create_dataset('rho-sol',data=rho)
            if globals()[sim+'par'].ldensity_nolog:
                grp.create_dataset('rho-num',data=var.rho)
            else:
                rho = np.exp(var.lnrho)
                grp.create_dataset('rho-num',data=rho)
            grp.create_dataset('uu-num',data=var.ux)
            grp.create_dataset('pp-num',data=var.pp)
            ee = var.TT/globals()[sim+'pars'].gamma*globals()[sim+'pars'].cp
            grp.create_dataset('ee-num',data=ee)
            grp.create_dataset('time',data=var.t)
            if globals()[sim+'par'].lmagnetic:
                grp.create_dataset('bb-num',data=var.bb[2])

figsize=[7.5,4.635255]

for sim in sims:
    png='.png'
    os.chdir(join(wrkdir,sim))
    if exists('data/allprocs/var.h5'):
        os.chdir('data/allprocs')
    else:
        os.chdir('data/proc0')
    varfiles = natural_sort(glob.glob('VAR*'))
    os.chdir(join(wrkdir,sim))
    with h5py.File(join(wrkdir,sim,'data','shocktest.h5'),'w') as hf:
        narray=globals()[sim+'gd'].x.size
        n1,n2=0,narray # option to reduce plot range in x
        plt.figure(figsize=figsize)
        for varfile in varfiles[-1]:
            grp = str.strip(varfile,'.h5')
            plt.semilogy(globals()[sim+'gd'].x[n1:n2]-globals()[sim+'gd'].x[0],
                         hf[grp]['rho-num'][n1:n2],'+',label='numeric',
                         alpha=0.5)
            plt.semilogy(globals()[sim+'gd'].x[n1:n2]-globals()[sim+'gd'].x[0],
                         hf[grp]['rho-sol'][n1:n2],'+',
                         ':',lw=4, color='orange',label='analytic')
            plt.text(globals()[sim+'gd'].x[n1+1]-globals()[sim+'gd'].x[0],
                     0.5*hf[grp]['rho-sol'][n1:n2].mean(),
                     r'$t={:.2}$'.format(hf[grp]['time'][()]))
            plt.xlabel(r'$x$ []')
            plt.ylabel(r'$n$ []')
            plt.legend(loc='lower left',framealpha=0.5)
            plt.savefig(join(wrkdir,sim,'fig',grp+'rho'+png))

            plt.figure(figsize=figsize)
            plt.plot(globals()[sim+'gd'].x[n1:n2]-globals()[sim+'gd'].x[0],
                             hf[grp]['uu-num'][n1:n2],'+',label='numeric',
                             alpha=0.5)
            plt.plot(globals()[sim+'gd'].x[n1:n2]-globals()[sim+'gd'].x[0],
                             hf[grp]['uu-sol'][n1:n2],'+',
                             ':',lw=4, color='orange',label='analytic')
            plt.text(globals()[sim+'gd'].x[n1+1]-globals()[sim+'gd'].x[0],
                     0.5*hf[grp]['uu-sol'][n1:n2].mean(),
                     r'$t={:.2}$'.format(hf[grp]['time'][()]))
            plt.xlabel(r'$x$ [kpc]')
            plt.ylabel(r'$u_x$ []')
            plt.legend(loc='upper left',framealpha=0.5)
            plt.savefig(join(wrkdir,sim,'fig',grp+'uu'+png))

            plt.figure(figsize=figsize)
            plt.plot(globals()[sim+'gd'].x[n1:n2]-globals()[sim+'gd'].x[0],
                             hf[grp]['ee-num'][n1:n2],'+',label='numeric',
                             alpha=0.5)
            plt.plot(globals()[sim+'gd'].x[n1:n2]-globals()[sim+'gd'].x[0],
                             hf[grp]['ee-sol'][n1:n2],'+',
                             ':',lw=4, color='orange',label='analytic')
            plt.text(globals()[sim+'gd'].x[n1+1]-globals()[sim+'gd'].x[0],
                     0.5*hf[grp]['ee-sol'][n1:n2].mean(),
                     r'$t={:.2}$'.format(hf[grp]['time'][()]))
            plt.xlabel(r'$x$ []')
            plt.ylabel(r'$e$ []')
            plt.legend(loc='lower left',framealpha=0.5)
            plt.savefig(join(wrkdir,sim,'fig',grp+'ee'+png))

            plt.figure(figsize=figsize)
            plt.plot(globals()[sim+'gd'].x[n1:n2]-globals()[sim+'gd'].x[0],
                             hf[grp]['pp-num'][n1:n2],'+',label='numeric',
                             alpha=0.5)
            plt.plot(globals()[sim+'gd'].x[n1:n2]-globals()[sim+'gd'].x[0],
                             hf[grp]['pp-sol'][n1:n2],'+',
                             ':',lw=4, color='orange',label='analytic')
            plt.text(globals()[sim+'gd'].x[n1+1]-globals()[sim+'gd'].x[0],
                     0.5*hf[grp]['pp-sol'][n1:n2].mean(),
                     r'$t={:.2}$'.format(hf[grp]['time'][()]))
            plt.xlabel(r'$x$ []')
            plt.ylabel(r'$p$ []')
            plt.legend(loc='lower left',framealpha=0.5)
            plt.savefig(join(wrkdir,sim,'fig',grp+'pp'+png))

            plt.figure(figsize=figsize)
            plt.plot(globals()[sim+'gd'].x[n1:n2]-globals()[sim+'gd'].x[0],
                             hf[grp]['uu-num'][n1:n2]/np.sqrt(
                             globals()[sim+'par'].gamma*hf[grp]['pp-num'][-1]/
                             hf[grp]['rho-num'][-1],'+',label='numeric',
                             alpha=0.5))
            plt.text(globals()[sim+'gd'].x[n1+1]-globals()[sim+'gd'].x[0],
                     1,
                     r'$t={:.2}$'.format(hf[grp]['time'][()]))
            plt.xlabel(r'$r$ [kpc]')
            plt.ylabel(r'Ms')
            plt.legend(loc='upper left',framealpha=0.5)
            plt.savefig(join(wrkdir,sim,'fig',grp+'Ms'+png))

            if globals()[sim+'par'].lmagnetic:
                plt.figure(figsize=figsize)
                plt.plot(globals()[sim+'gd'].x[n1:n2]-globals()[sim+'gd'].x[0],
                             hf[grp]['bb-num'][n1:n2],'+',label='numeric',
                             alpha=0.5)
                plt.text(globals()[sim+'gd'].x[n1+1]-globals()[sim+'gd'].x[0],
                         0.5*hf[grp]['bb-num'][n1:n2].mean(),
                         r'$t={:.2}$'.format(hf[grp]['time'][()]))
                plt.xlabel(r'$x$ []')
                plt.ylabel(r'$B_\perp$ []')
                plt.legend(loc='lower left',framealpha=0.5)
                plt.savefig(join(wrkdir,sim,'fig',grp+'Bz'+png))
