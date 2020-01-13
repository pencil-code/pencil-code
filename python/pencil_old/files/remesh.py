# remesh.py
#
# Read VAR files from existing mature simulation located in the source_path 
# directory to be remeshed into newly initialised target simulation of new 
# grid and processor layout located in the target_path directory. 
# 
# First set up new run compile and start it.
# To interpolate the old data into new shape call interp_var with string
# arguments specifying source_path and target_path directories. 
# Then use the resulting fnew array and target_path as arguments for
# distribute_fort to overwrite the var.dat in the process tree
# 
#
# NB: to include particles or additional farray variables the scripts will
#     need to be modified
#     the interpolation does not handle switching the system of coordinates
#     for periodic boundaries the arrays are stretched or shrunk to maintain
#     continuity on the boundaries. This is the default.
#     Not implemented: clipping the data or extending e.g. repeat the array in
#     the periodic directions
#     Not implemented reading or writing farray in hdf5 format
#
# Author: F. Gent (fred.gent.ncl@gmail.com).
#
import numpy as np
import os
from scipy.interpolate import interp1d 
import pencil as pc
from scipy.io import FortranFile as ftn

def interp_var(
                oldvar='var.dat',
                newvar='var.dat',
                source='old_run', target='new_run', 
                source_path=None, target_path=None,
                xlim=None, ylim=None, zlim=None,
                xrepeat=1,
                yrepeat=1,
                zrepeat=1,
                time=None,
                deltay=None,
                arrs=None,
                innghosts=3,
                outghosts=3
              ):
    """ load var file to be interpolated from old simulation and var file from
    started new simulation of correct shape and processor layout
    interpolate the old snapshot onto a mesh of the new shape and old grid 
    limits
    then rescale new grid and deltay to new limits.
    """
    if source_path==None:
        localdir=os.getcwd()+'/'
        if not os.path.exists(localdir+source):
            print('error: source_path must be specified as string')
            return 
        os.chdir(localdir+source)
        print('loading data from '+localdir+source)
        fold=pc.reads.var(oldvar,quiet=True)
    else:
        if not os.path.exists(source_path):
            print('error: source_path does not exist')
            return 
        os.chdir(source_path)
        print('loading data from '+source_path)
        fold=pc.read.var(oldvar,quiet=True)
    if target_path==None:
        if not os.path.exists(localdir+target):
            print('error: target_path must be specified as string')
            return 
        os.chdir(localdir+target)
        print('loading data from '+localdir+target)
        fnew=pc.read.var(newvar,quiet=True)
    else:
        if not os.path.exists(target_path):
            print('error: target_path does not exist')
            return 
        os.chdir(target_path)
        print('loading data from '+target_path)
        fnew=pc.read.var(newvar,quiet=True)
    if xlim==None:
        xlim=[fold.x.min(),fold.x.max()]
    if ylim==None:
        ylim=[fold.y.min(),fold.y.max()]
    if zlim==None:
        zlim=[fold.z.min(),fold.z.max()]
    if arrs==None:
        arrs=['uu','rho','lnrho','ss','aa','shock','netheat','cooling']
    iarr=0
    x=np.linspace(xlim[0],xlim[1],fnew.x.size)
    y=np.linspace(ylim[0],ylim[1],fnew.y.size)
    z=np.linspace(zlim[0],zlim[1],fnew.z.size)
    for arr in arrs:
        if hasattr(fold,arr):
            print('interpolating '+arr)
            tmp = fold.__getattribute__(arr)
            if xrepeat == 1:
                intmpx = interp1d(fold.x,tmp,axis=-1)
                tmp = intmpx(x)
            if yrepeat == 1:
                intmpy = interp1d(fold.y,tmp,axis=-2)
                tmp=intmpy(y)
            if zrepeat == 1:
                intmpz = interp1d(fold.z,tmp,axis=-3)
                tmp = intmpz(z)
                fnew.__getattribute__(arr)[:] = tmp
            if len(tmp.shape)==4:
                fnew.f[iarr:iarr+3] = tmp
                iarr += 3
            else:
                fnew.f[iarr] = tmp
                iarr += 1
    if hasattr(fold,'deltay'):
        fnew.deltay = fold.deltay*(fnew.y.size-2*innghosts)/float(
                                   fold.y.size-2*outghosts)
    if not time==None:
        fnew.t=time
    return fnew

def distribute_fort(
                     farray,
                     target_path,
                     filename='var.dat',
                     arrs=[
                            'uu',
                            'rho',
                            'lnrho',
                            'ss',
                            'lnTT',
                            'aa',
                            'bb',
                            'ecr',
                            'fcr',
                            'shock',
                            'netheat',
                            'cooling'
                          ],
                     persist=None,
                     nghosts=3
                   ):
    """ take new shape farray and write fortran binary files to new processor
    array
    if farray contains variables not listed in arr, submit arr of correct 
    strings listed in farray order
    if fparray to be included the f.write list will need to be extended below
    and parr array of strings inserted to header
    if persist not None submit tuple triplets of value, type integer id
    """
    os.chdir(target_path)
    datadir=os.getcwd()+'/data/'
    if not os.path.exists(datadir):
        print('error: target_path does not exist')
    dim=pc.read.dim()
    if dim.precision=='D':
        prec=np.float64
        intp=np.int64
    else:
        prec=np.float32
        intp=np.int32
    nx=int((farray.x.size-2*nghosts)/dim.nprocx)
    ny=int((farray.y.size-2*nghosts)/dim.nprocy)
    nz=int((farray.z.size-2*nghosts)/dim.nprocz)
    if hasattr(farray,'deltay'):
        tmp_arr=np.zeros(nx+ny+nz+6*nghosts+5)
    else:
        tmp_arr=np.zeros(nx+ny+nz+6*nghosts+4)
    for ipx in range(0,dim.nprocx):
        ix=np.arange(nx+2*nghosts)+ipx*nx
        tmpx=farray.f[:,:,:,ix]
        for ipy in range(0,dim.nprocy):
            iy=np.arange(ny+2*nghosts)+ipy*ny
            tmpy=tmpx[:,:,iy]
            for ipz in range(0,dim.nprocz):
                iproc=ipx+dim.nprocx*ipy+dim.nprocx*dim.nprocy*ipz
                iz=np.arange(nz+2*nghosts)+ipz*nz
                f = ftn(datadir+'proc'+str(iproc)+'/'+filename, 'w') 
                #f = ftn(datadir+'proc'+str(iproc)+'/'+filename, 'w',        
                #        header_dtype=np.uint64)
                print('writing '+datadir+'proc'+str(iproc)+'/'+filename)
                f.write_record(tmpy[:,iz])
                tmp_arr[0]=farray.t
                tmp_arr[1:nx+1+2*nghosts]=farray.x[ix]
                tmp_arr[1+nx+2*nghosts:ny+nx+1+4*nghosts]=farray.y[iy]
                tmp_arr[1+nx+ny+4*nghosts:nx+ny+nz+1+6*nghosts]=farray.z[iz]
                tmp_arr[1+nx+ny+nz+6*nghosts:nx+ny+nz+2+6*nghosts]=farray.dx
                tmp_arr[2+nx+ny+nz+6*nghosts:nx+ny+nz+3+6*nghosts]=farray.dy
                tmp_arr[3+nx+ny+nz+6*nghosts:nx+ny+nz+4+6*nghosts]=farray.dz
                if hasattr(farray,'deltay'):
                    tmp_arr[4+nx+ny+nz+6*nghosts:nx+ny+nz+5+6*nghosts]=farray.deltay
                f.write_record(tmp_arr)
                #f.write_record(np.float32(farray.t))
                #f.write_record(farray.x[ix])
                #f.write_record(farray.y[iy])
                #f.write_record(farray.z[iz])
                #f.write_record(farray.dx)
                #f.write_record(farray.dy)
                #f.write_record(farray.dz)
                #if hasattr(farray,'deltay'):
                #    f.write_record(farray.deltay)
                if not persist==None:
                    for pers in persist:
                        f.write_record(pers[0])
                        f.write_record(pers[1])
                f.close()

#set of persistent variables used by interstellar, listed in record_types.h
#if required replace pco.pers with your list in local_remesh.py below 
pers=[
      (253, 0.0), 
      (254, 0.0), 
      (255, 0.0), 
      (256, 0.0), 
      (257,   1), 
      (258,   1), 
      (260, 0.0), 
      (261, 0.0) 
     ]

""" Example of use to create remesh of mature_run on grid 32x32x64 on 1x2x4
    processors to new_run on grid 40x40x80 on 2x2x4 processors
    > cd /path/to/mature_run
    > pc_newrun /path/to/new_run
    > cd /path/to/new_run
    > vi src/cparam.local
    edit the file to the new grid and processor layout
    > vi start.in
    edit xyz0 and xyz1 or Lxyz if the dimensions of the box are to change
    at this stage do not add new physics, although this can be done after remesh
    > pc_build
    > pc_run start
    > vi local_remesh.py
file to contain:

import pencil_old as pco
fnew=pco.interp_var(
                   target_path='/path/to/mature_run',
                   source_path='/path/to/new_run'
                  )
pco.distribute_fort(fnew,'/path/to/new_run')
##or if the run uses persistent variables comment above line and uncomment 
##below. pco.pers will need to be replaced by the appropriate variables
#pco.distribute_fort(
                    fnew,
                    '/path/to/new_run',
                    persist=pco.pers
                   )

    > python local_remesh.py
    Done!

    If you prefer the last stage can be done in stages using a python gui
"""    

