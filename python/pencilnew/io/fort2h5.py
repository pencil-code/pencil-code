# fort2h5.py
#
# Read existing Fortran unformatted simulation data and write as hdf5.
#
#
# Author: F. Gent (fred.gent.ncl@gmail.com).
#
"""
Contains the functions to read old Fortran binary data and write snapshots
in hdf5 (VAR*.h5 files).
"""

def sim2h5(newdir='.', olddir='.', varfile_names=None,
           todatadir='data/allprocs', fromdatadir='data',
           precision=None, nghost=3, lpersist=False,
           x=None, y=None, z=None, lshear=False,
           lremove_old=False, execute=False
          ):

    """
    Copy a simulation object written in Fortran binary to hdf5.
    The default is to copy all snapshots from/to the current simulation
    directory. Optionally the old files can be removed to 

    call signature:

    sim2h5(newdir='.', olddir='.', varfile_names=None,
           todatadir='data/allprocs', fromdatadir='data',
           precision=None, nghost=3, lpersist=False,
           x=None, y=None, z=None, lshear=False,
           lremove_old=False, execute=False)

    Keyword arguments:

    *newdir*:
      String path to simulation destination directory.
      Path may be relative or absolute.

    *newdir*:
      String path to simulation destination directory.
      Path may be relative or absolute.

    *varfile_names*:
      A list of names of the snapshot files to be written, e.g. VAR0
      If None all varfiles in olddir+'/data/proc0/' will be converted

    *todatadir*:
      Directory to which the data is stored.

    *fromdatadir*:
      Directory from which the data is collected.

    *precision*:
      Single 'f' or double 'd' precision for new data.
      
    *nghost*:
      Number of ghost zones.
      TODO: handle switching size of ghost zones.

    *lpersist*:
      option to include persistent variables from snapshots.

    *xyz*:
      xyz arrays of the domain with ghost zones.
      This will normally be obtained from Grid object, but facility to
      redefine an alternative grid value.

    *lshear*:
      Flag for the shear.

    *lremove_old*:
      If True the old snapshots will be deleted once the new snapshot has
      been saved.
      A warning is given without execution to avoid unintended removal.

    *execute*:
      optional confirmation required if lremove_old.

    """

    import os
    import numpy as np
    import h5py
    import glob
    from .. import read
    from .. import sim

    #test if simulation directories
    os.chdir(olddir)
    if not sim.is_sim_dir():
        print("ERROR: Directory ("+olddir+") needs to be a simulation")
        return -1
    if newdir != olddir:
        os.chdir(newdir)
        if not sim.is_sim_dir():
            print("ERROR: Directory ("+newdir+") needs to be a simulation")
            return -1
    #
    if lremove_old:
        if not execute:
            os.chdir(olddir)
            print("WARNING: Are you sure you wish to remove the Fortran binary"+
                  " files from \n"+
                  os.getcwd()+".\n"+
                  "Set execute=True to proceed.")
            return -1

    os.chdir(olddir)
    if varfile_names == None:
        os.chdir(fromdatadir+'/proc0')
        varfile_names = glob.glob('VAR*') 
        os.chdir(olddir)
    gkeys = ['x', 'y', 'z', 'Lx', 'Ly', 'Lz', 'dx', 'dy', 'dz',
             'dx_1', 'dy_1', 'dz_1', 'dx_tilde', 'dy_tilde', 'dz_tilde',
            ]
    grid=read.grid(quiet=True)
    for key in gkeys:
        if not key in grid.__dict__.keys():
            print("ERROR: key "+key+" missing from grid")
            return -1 
    #obtain the settings from the old simulation
    settings={}
    skeys = ['l1', 'l2', 'm1', 'm2', 'n1', 'n2',
             'nx', 'ny', 'nz', 'mx', 'my', 'mz',
             'nprocx', 'nprocy', 'nprocz',
             'maux', 'mglobal', 'mvar', 'precision',
            ]
    olddim = read.dim()
    for key in skeys:
        settings[key]=olddim.__getattribute__(key)
    settings['nghost']=nghost
    if precision != None:
        settings['precision']=precision
    #obtain physical units from old simulation
    ukeys = ['length', 'velocity', 'density', 'magnetic', 'time',
                 'temperature', 'flux', 'energy', 'mass', 'system',
                ]
    param = read.param()
    param.__setattr__('unit_mass',param.unit_density*param.unit_length**3)
    param.__setattr__('unit_energy',param.unit_mass*param.unit_velocity**2)
    param.__setattr__('unit_time',param.unit_length/param.unit_velocity)
    param.__setattr__('unit_flux',param.unit_mass/param.unit_time**3)
    param.unit_system=param.unit_system.encode()
    #index list for variables in f-array
    indx=read.index()


    #check consistency between Fortran binary and h5 data
    os.chdir(newdir)
    dim = read.dim()
        try:
            dim.mvar == settings['mvar']
            dim.mx   == settings['mx']
            dim.my   == settings['my']
            dim.mz   == settings['mz']
        except ValueError:
            print("ERROR: new simulation dimensions do not match.")
            return -1 
    #move var.h5 out of the way, if it exists for reading binary
    if os.path.exists(todatadir+'/var.h5'):
        cmd='mv '+todatadir+'/var.h5 '+todatadir+'/var.bak'
        os.system(cmd)

    #proceed to copy each snapshot in varfile_names
    for file_name in varfile_names:
        #load Fortran binary snapshot
        os.chdir(olddir)
        var = read.var(file_name, datadir=fromdatadir, quiet=quiet,
                       lpersist=lpersist
                      )
        try:
            var.deltay
            lshear = True
        except
            lshear = False

        if lpersist:
            persist = {}
            for key in read.record_types.keys():
                try:
                    persist[key] = var.__getattribute__(key)[()]
                except:
                    continue 
        else:
            persist = None        
        #write data to h5
        os.chdir(newdir)
        write_h5_snapshot(var.f, file_name=file_name, datadir=todatadir,
                          precision=precision, nghost=nghost, persist=persist,
                          settings=settings, param=param, grid=grid,
                          lghosts=True, indx=indx, t=var.t, x=x, y=y, z=z,
                          lshear=False)
        if lremove_old:
            os.chdir(olddir)
            cmd="rm -f "+os.path.join(fromdatadir, 'proc*', file_name)
            os.system(cmd)
