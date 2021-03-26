#!/usr/bin/env python3

import os 
from .. import read
import numpy as np
import pylab as plt
from scipy.io import FortranFile

def animate_slices(field='uu1',datadir='data/',proc=-1,extension='xz',format='native',
                tmin=0.,tmax=1.e38,wait=0.,amin=0.,amax=1.,transform='',oldfile=False):
    """
    read 2D slice files and assemble an animation.

    Options:

     field --- which variable to slice
     datadir --- path to data directory
     proc --- an integer giving the processor to read a slice from
     extension --- which plane of xy,xz,yz,Xz. for 2D this should be overwritten.
     format --- endian. one of little, big, or native (default)
     tmin --- start time
     tmax --- end time
     amin --- minimum value for image scaling
     amax --- maximum value for image scaling
     transform --- insert arbitrary numerical code to modify the slice
     wait --- pause in seconds between animation slices
    """
    
    datadir = os.path.expanduser(datadir)
    if proc < 0:
        filename = datadir+'/slice_'+field+'.'+extension
    else:
        filename = datadir+'/proc'+str(proc)+'/slice_'+field+'.'+extension

    # global dim
    #param = read_param(datadir)
    param = read.param(datadir)

    #dim = read_dim(datadir,proc) 
    dim = read.dim(datadir,proc) 
    if dim.precision == 'D':
        precision = 'd'
    else:
        precision = 'f'

    # set up slice plane
    if (extension == 'xy' or extension == 'Xy'):
        hsize = dim.nx
        vsize = dim.ny
    if (extension == 'xz'):
        hsize = dim.nx
        vsize = dim.nz
    if (extension == 'yz'):
        hsize = dim.ny
        vsize = dim.nz
    plane = np.zeros((vsize,hsize),dtype=precision)

    infile = FortranFile(filename)

    ax = plt.axes()
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_ylim

    image = plt.imshow(plane,vmin=amin,vmax=amax)

    # for real-time image display
    manager = plt.get_current_fig_manager()
    manager.show()

    ifirst = True
    islice = 0
    while 1:
        try:
            raw_data = infile.read_record(dtype=precision)
        except ValueError:
            break
        except TypeError:
            break

        if oldfile:
            t = raw_data[-1]
            plane = raw_data[:-1].reshape(vsize,hsize)
        else:
            slice_z2pos = raw_data[-1]
            t = raw_data[-2]
            plane = raw_data[:-2].reshape(vsize,hsize)
        
        if transform:
            exec('plane = plane'+transform)

        if (t > tmin and t < tmax):
            title = 't = %11.3e' % t
            ax.set_title(title)
            image.set_data(plane)
            manager.canvas.draw()
            
            if ifirst:
                print("----islice----------t---------min-------max-------delta")
            print("%10i %10.3e %10.3e %10.3e %10.3e" % (islice,t,plane.min(),plane.max(),plane.max()-plane.min()))
                
            ifirst = False
            islice += 1

            sleep(wait)

    infile.close() 
