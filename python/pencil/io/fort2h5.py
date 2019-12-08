# fort2h5.py
#
# Read existing Fortran unformatted simulation data and write as hdf5.
#
#
# Author: F. Gent (fred.gent.ncl@gmail.com).
#
"""
Contains the functions to read old Fortran binary simulation data and
write snapshots in hdf5 (data/allprocs/VAR*.h5 files), video slices to
data/slices/uu1_xy.h5, etc. and averages to data/averages/xy.h5, etc
"""

def var2h5(newdir, olddir, varfile_names, todatadir, fromdatadir,
           precision, lpersist, quiet, nghost, settings, param, grid,
           x, y, z, lshear, lremove_old_snapshots, indx
          ):

    """
    Copy a simulation snapshot set written in Fortran binary to hdf5.

    call signature:

    var2h5(newdir, olddir, varfile_names, todatadir, fromdatadir,
           precision, lpersist, quiet, nghost, settings, param, grid,
           x, y, z, lshear, lremove_old_snapshots, indx
          )

    Keyword arguments:

    *newdir*:
      String path to simulation destination directory.

    *olddir*:
      String path to simulation destination directory.

    *varfile_names*:
      A list of names of the snapshot files to be written, e.g. VAR0.

    *todatadir*:
      Directory to which the data is stored.

    *fromdatadir*:
      Directory from which the data is collected.

    *precision*:
      Single 'f' or double 'd' precision for new data.

    *lpersist*:
      option to include persistent variables from snapshots.

    *quiet*
      Option not to print output.

    *nghost*:
      Number of ghost zones.

    *settings*
      simulation properties.

    *param*
      simulation Param object.

    *grid*
      simulation Grid object.

    *xyz*:
      xyz arrays of the domain with ghost zones.

    *lshear*:
      Flag for the shear.

    *lremove_old_snapshots*:
      If True the old snapshots will be deleted once the new snapshot has
      been saved.

    *indx*
      List of variable indices in the f-array.

    """

    import os
    import numpy as np
    import h5py
    import glob
    from .. import read
    from .. import sim
    from . import write_h5_snapshot

    #move var.h5 out of the way, if it exists for reading binary
    if os.path.exists(todatadir+'/var.h5'):
        cmd='mv '+todatadir+'/var.h5 '+todatadir+'/var.bak'
        os.system(cmd)

    #proceed to copy each snapshot in varfile_names
    for file_name in varfile_names:
        #load Fortran binary snapshot
        print('saving '+file_name)
        os.chdir(olddir)
        var = read.var(file_name, datadir=fromdatadir, quiet=quiet,
                       lpersist=lpersist
                      )
        try:
            var.deltay
            lshear = True
        except:
            lshear = False

        if lpersist:
            persist = {}
            for key in read.record_types.keys():
                try:
                    persist[key] = var.__getattribute__(key)[()]
                    if (type(persist[key][0])==str):
                        persist[key][0] = var.__getattribute__(key)[0].encode()
                except:
                    continue
        else:
            persist = None
        #write data to h5
        os.chdir(newdir)
        write_h5_snapshot(var.f, file_name=file_name, datadir=todatadir,
                          precision=precision, nghost=nghost,
                          persist=persist,
                          settings=settings, param=param, grid=grid,
                          lghosts=True, indx=indx, t=var.t, x=x, y=y, z=z,
                          lshear=lshear)
        if lremove_old_snapshots:
            os.chdir(olddir)
            cmd = "rm -f "+os.path.join(fromdatadir, 'proc*', file_name)
            os.system(cmd)
    os.chdir(olddir)
    var = read.var('var.dat', datadir=fromdatadir, quiet=quiet,
                   lpersist=lpersist
                  )
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
    write_h5_snapshot(var.f, file_name='var', datadir=todatadir,
                      precision=precision, nghost=nghost,
                      persist=persist,
                      settings=settings, param=param, grid=grid,
                      lghosts=True, indx=indx, t=var.t, x=x, y=y, z=z,
                      lshear=lshear)
    if lremove_old_snapshots:
        os.chdir(olddir)
        cmd = "rm -f "+os.path.join(fromdatadir, 'proc*', 'var.dat')
        os.system(cmd)

def slices2h5(newdir, olddir, grid,
              todatadir='data/slices', fromdatadir='data',
              precision='d', quiet=True, lremove_old_slices=False):

    """
    Copy a simulation set of video slices written in Fortran binary to hdf5.

    call signature:

    slices2h5(newdir, olddir, grid,
              todatadir='data/slices', fromdatadir='data',
              precision='d', quiet=True, lremove_old_slices=False)

    Keyword arguments:

    *newdir*:
      String path to simulation destination directory.

    *olddir*:
      String path to simulation destination directory.

    *grid*
      simulation Grid object.

    *todatadir*:
      Directory to which the data is stored.

    *fromdatadir*:
      Directory from which the data is collected.

    *precision*:
      Single 'f' or double 'd' precision for new data.

    *quiet*
      Option not to print output.

    *lremove_old_slices*:
      If True the old video slices will be deleted once the new slices have
      been saved.
    """

    import os
    from .. import read
    from .. import sim
    from . import write_h5_slices

    #copy old video slices to new h5 sim
    os.chdir(olddir)
    vslice = read.slices()
    #identify the coordinates and positions of the slices
    coordinates = {}
    positions = {}
    readlines1 = open('data/slice_position.dat','r').readlines()
    readlines2 = open('data/proc0/slice_position.dat','r').readlines()
    lines1, lines2 = [],[]
    for line in readlines1:
        lines1.append(int(line.split(' ')[-1].split('\n')[0]))
    """In newer binary sims lines2 obtains 7 strings below, but older sims may
    only yield the integer coordinates, so lines2 is hardcoded. The current
    version of the Pencil Code has 7 potential slices, but earlier versions
    may not. If your sim does not conform to this arrangement edit/copy this
    module and set lines1 and lines2 manually from data/slice_position.dat and
    the extensions present in your slice_*.xy  etc.
    """
    #check simulation includes the slice keys in data/proc*/slice_position.dat
    try:
        int(int(readlines2[0].split(' ')[-1].split('\n')[0]))
        lines2=['xy', 'xy2', 'xy3', 'xy4', 'xz', 'xz2', 'yz']
    except:
        for line in readlines2:
            lines2.append(line.split(' ')[-1].split('\n')[0].lower())
    #check if number of slice options as expected
    try:
        len(lines1)==7
    except ValueError:
        print("ERROR: slice keys and positions must be set see lines 212...")
        return -1
    for key, num in zip(lines2, lines1):
        if num > 0:
            if 'xy' in key:
                positions[key] = grid.z[num-1]
            if 'xz' in key:
                positions[key] = grid.y[num-1]
            if 'yz' in key:
                positions[key] = grid.x[num-1]
            coordinates[key] = num
    #write new slices in hdf5
    os.chdir(newdir)
    write_h5_slices(vslice, coordinates, positions, datadir=todatadir,
                    precision=precision, quiet=quiet)
    if lremove_old_slices:
        os.chdir(olddir)
        cmd = "rm -f "+os.path.join(fromdatadir, 'proc*', 'slice_*')
        os.system(cmd)

def aver2h5(newdir, olddir,
            todatadir='data/averages', fromdatadir='data', l2D=True,
            precision='d', quiet=True, lremove_old_averages=False):

    """
    Copy a simulation set of video slices written in Fortran binary to hdf5.

    call signature:

    aver2h5(newdir, olddir,
            todatadir='data/slices', fromdatadir='data', l2D=True,
            precision='d', quiet=True, lremove_old_averages=False)

    Keyword arguments:

    *newdir*:
      String path to simulation destination directory.

    *olddir*:
      String path to simulation destination directory.

    *todatadir*:
      Directory to which the data is stored.

    *fromdatadir*:
      Directory from which the data is collected.

    *l2D*
     Option to include 2D averages if the file sizes are not too large

    *precision*:
      Single 'f' or double 'd' precision for new data.

    *quiet*
      Option not to print output.

    *lremove_old_averages*:
      If True the old averages data will be deleted once the new h5 data
      has been saved.
    """

    import os
    from .. import read
    from .. import sim
    from . import write_h5_averages

    #copy old 1D averages to new h5 sim
    os.chdir(olddir)
    av = read.aver()
    os.chdir(newdir)
    for key in av.__dict__.keys():
        if not key in 't':
            write_h5_averages(av, file_name=key, datadir=todatadir,
                              precision=precision, quiet=quiet)
    if lremove_old_averages:
        os.chdir(olddir)
        cmd = "rm -f "+os.path.join(fromdatadir, '*averages.dat')
        os.system(cmd)
    if l2D:
       plane_list = []
       os.chdir(olddir)
       if os.path.exists('xaver.in'):
           plane_list.append('x')
       if os.path.exists('yaver.in'):
           plane_list.append('y')
       if os.path.exists('zaver.in'):
           plane_list.append('z')
       if len(plane_list) > 0:
           for key in plane_list:
               os.chdir(olddir)
               av = read.aver(plane_list=key)
               os.chdir(newdir)
               write_h5_averages(av, file_name=key, datadir=todatadir,
                                 precision=precision, quiet=quiet)
    if lremove_old_averages:
        os.chdir(olddir)
        cmd = "rm -f "+os.path.join(fromdatadir, 'proc*', '*averages.dat')
        os.system(cmd)

def sim2h5(newdir='.', olddir='.', varfile_names=None,
           todatadir='data/allprocs', fromdatadir='data',
           precision='d', nghost=3, lpersist=False,
           x=None, y=None, z=None, lshear=False,
           lremove_old_snapshots=False, lremove_old_slices=False,
           lremove_old_averages=False, execute=False, quiet=True,
           l2D=True, lvars=True, lvids=True, laver=True
          ):

    """
    Copy a simulation object written in Fortran binary to hdf5.
    The default is to copy all snapshots from/to the current simulation
    directory. Optionally the old files can be removed to 

    call signature:

    sim2h5(newdir='.', olddir='.', varfile_names=None,
           todatadir='data/allprocs', fromdatadir='data',
           precision='d', nghost=3, lpersist=False,
           x=None, y=None, z=None, lshear=False,
           lremove_old_snapshots=False, lremove_old_slices=False,
           lremove_old_averages=False, execute=False, quiet=True,
           l2D=True, lvars=True, lvids=True, laver=True)

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
    from . import write_h5_grid

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
    lremove_old = lremove_old_snapshots or\
                  lremove_old_slices or lremove_old_averages
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
    settings['precision']=precision.encode()
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
    print('precision is ',precision)
    if lvars:
        var2h5(newdir, olddir, varfile_names, todatadir, fromdatadir,
               precision, lpersist, quiet, nghost, settings, param, grid,
               x, y, z, lshear, lremove_old_snapshots, indx)
    #copy old video slices to new h5 sim
    if lvids:
        slices2h5(newdir, olddir, grid,
                  todatadir='data/slices', fromdatadir='data',
                  precision=precision, quiet=quiet,
                  lremove_old_slices=lremove_old_slices)
    #copy old averages data to new h5 sim
    if laver:
        aver2h5(newdir, olddir,
                todatadir='data/averages', fromdatadir='data', l2D=l2D,
                precision=precision, quiet=quiet,
                lremove_old_averages=lremove_old_averages)
    #check some critical sim files are present for new sim without start
    #construct grid.h5 sim information if requied for new h5 sim
    os.chdir(newdir)
    write_h5_grid(file_name='grid', datadir='data', precision=precision,
                  nghost=nghost, settings=settings, param=param, grid=grid,
                  unit=None, quiet=quiet)
    source_file = os.path.join(olddir,fromdatadir,'proc0/varN.list')
    target_file = os.path.join(newdir,todatadir,'varN.list')
    if os.path.exists(source_file):
        cmd='cp '+source_file+' '+target_file
        os.system(cmd)
    items=['def_var.pro', 'index.pro', 'jobid.dat', 'param.nml',
           'particle_index.pro', 'pc_constants.pro', 'pointmass_index.pro',
           'pt_positions.dat', 'sn_series.dat', 'svnid.dat', 'time_series.dat',
           'tsnap.dat', 'tspec.dat', 'tvid.dat', 'var.general',
           'variables.pro', 'varname.dat']
    for item in items:
        source_file = os.path.join(olddir,fromdatadir,item)
        target_file = os.path.join(newdir,fromdatadir,item)
        if os.path.exists(source_file):
            if not os.path.exists(target_file):
                cmd='cp '+source_file+' '+target_file
                os.system(cmd)
