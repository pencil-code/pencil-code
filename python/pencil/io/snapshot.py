# snapshot.py
#
# Generate a snapshot or initial condition from a numpy array.
#
#
# Author: S. Candelaresi (iomsn1@gmail.com).
#
"""
Contains the functions to generate a snapshot (VAR files).
"""


def write_snapshot(snapshot, file_name='VAR0', datadir='data',
                   nprocx=1, nprocy=1, nprocz=1, precision='d', nghost=3,
                   t=None, x=None, y=None, z=None, lshear=False):
    """
    Write a snapshot given as numpy array.

    call signature:

    write_snapshot(snapshot, file_name='VAR0', datadir='data',
                   nprocx=1, nprocy=1, nprocz=1, precision='d', nghost=3)

    Keyword arguments:

    *snapshot*:
      Numpy array containing the snapshot. Must be of shape [nvar, nz, ny, nx]
      (without boundaries).

    *file_name*:
      Name of the snapshot file to be written, e.g. VAR0 or var.dat.

    *datadir*:
      Directory where the data is stored.

    *nproc[xyz]*:
      Number of processors in xyz.

    *precision*:
      Single 'f' or double 'd' precision.

    *nghost*:
      Number of ghost zones.

    *t*:
      Time of the snapshot.

    *xyz*:
      xyz arrays of the domain without ghost zones.

    *lshear*:
      Flag for the shear.
    """

    import os.path
    import numpy as np
    from scipy.io import FortranFile

    # Determine the shape of the input snapshot.
    nx = snapshot.shape[3]
    ny = snapshot.shape[2]
    nz = snapshot.shape[1]

    # Determine the precision used.
    if precision == 'f':
        data_type = np.float32
    elif precision == 'd':
        data_type = np.float64
    else:
        print("ERROR: Precision {0} not understood.".format(precision)+
              " Must be either 'f' or 'd'")

    # Check that the shape does not conflict with the proc numbers.
    if (nx%nprocx > 0) or (ny%nprocy > 0) or (nz%nprocz > 0):
        print('ERROR: Shape of the input array is not compatible with the '+
              'cpu layout. Make sure that nproci devides ni.')
        return -1

    # Check the shape of the xyz arrays and generate them if they don't exist.
    if not x is None:
        if len(x) != nx:
            print(
                'ERROR: x array is incompatible with the shape of snapshot.')
    else:
        x = np.arange(0, nx)
    if not y is None:
        if len(y) != ny:
            print(
                'ERROR: y array is incompatible with the shape of snapshot.')
    else:
        y = np.arange(0, ny)
    if not z is None:
        if len(z) != nz:
            print(
                'ERROR: z array is incompatible with the shape of snapshot.')
    else:
        z = np.arange(0, nz)

    # Add ghost zones to the xyz arrays.
    x_ghost = np.zeros(nx+2*nghost, dtype=data_type)
    x_ghost[nghost:-nghost] = x
    y_ghost = np.zeros(ny+2*nghost, dtype=data_type)
    y_ghost[nghost:-nghost] = y
    z_ghost = np.zeros(nz+2*nghost, dtype=data_type)
    z_ghost[nghost:-nghost] = z

    # Define the deltas.
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    dz = z[1] - z[0]

    # Define a time.
    if t is None:
        t = 0

    # Create the data directories if they don't exist.
    iproc = 0
    datadir_list = os.listdir(datadir)
    for iproc in range(nprocx*nprocy*nprocz):
        if not 'proc{0}'.format(iproc) in datadir_list:
            os.mkdir(os.path.join(datadir, 'proc{0}'.format(iproc)))

    # Determine the precision used.
    snapshot = snapshot.astype(data_type)

    # Add ghost zones.
    snapshot_ghosts = np.zeros(
                      [snapshot.shape[0],
                           nz+2*nghost, ny+2*nghost, nx+2*nghost],
                      dtype=data_type)
    snapshot_ghosts[:, nghost:-nghost,
                       nghost:-nghost,
                       nghost:-nghost] = snapshot

    # Write the data.
    iproc = 0
    for ipz in range(nprocz):
        for ipy in range(nprocy):
            for ipx in range(nprocx):
                destination_file = FortranFile(
                    os.path.join(datadir,
                                 'proc{0}'.format(iproc), file_name), 'w')
                # Construct and write the snapshot for this cpu.
                snapshot_cpu = snapshot_ghosts[:,
                    ipz*int(nz/nprocz):(ipz+1)*int(nz/nprocz)+2*nghost,
                    ipy*int(ny/nprocy):(ipy+1)*int(ny/nprocy)+2*nghost,
                    ipx*int(nx/nprocx):(ipx+1)*int(nx/nprocx)+2*nghost]
                destination_file.write_record(snapshot_cpu)
#                snapshot_cpu.tofile(destination_file.format(iproc))
                # Construct and write the meta data for this cpu.
                x_cpu = x_ghost[ipx*int(nx/nprocx):(ipx+1)*
                                    int(nx/nprocx)+2*nghost]
                y_cpu = y_ghost[ipy*int(ny/nprocy):(ipy+1)*
                                    int(ny/nprocy)+2*nghost]
                z_cpu = z_ghost[ipz*int(nz/nprocz):(ipz+1)*
                                    int(nz/nprocz)+2*nghost]
                meta_data_cpu = np.array([t], dtype=data_type)
                meta_data_cpu = np.append(meta_data_cpu, x_cpu)
                meta_data_cpu = np.append(meta_data_cpu, y_cpu)
                meta_data_cpu = np.append(meta_data_cpu, z_cpu)
                meta_data_cpu = np.append(meta_data_cpu, [dx, dy, dz])
                if lshear:
                    meta_data_cpu = np.append(meta_data_cpu, lshear)
                meta_data_cpu = meta_data_cpu.astype(data_type)
                destination_file.write_record(meta_data_cpu)
                destination_file.close()
                iproc += 1

    return 0

def write_h5_snapshot(snapshot, file_name='VAR0', datadir='data/allprocs',
                   precision='d', nghost=3, persist=None, settings=None,
                   param=None, grid=None, lghosts=False, indx=None,
                   unit=None, t=None, x=None, y=None, z=None,
                   quiet=True, lshear=False, driver=None, comm=None):
    """
    Write a snapshot given as numpy array.
    We assume by default that a run simulation directory has already been
    constructed and start completed successfully in h5 format so that
    files dim, grid and param files are already present.
    If not the contents of these will need to be supplied as dictionaries
    along with persist if included.

    call signature:

    write_h5_snapshot(snapshot, file_name='VAR0', datadir='data/allprocs',
                   precision='d', nghost=3, persist=None, settings=None,
                   param=None, grid=None, lghosts=False, indx=None,
                   unit=None, t=None, x=None, y=None, z=None,
                   quiet=True, lshear=False, driver=None, comm=None)

    Keyword arguments:

    *snapshot*:
      Numpy array containing the snapshot.
      Must be of shape [nvar, nz, ny, nx] without boundaries or.
      Must be of shape [nvar, mz, my, mx] with boundaries for lghosts=True.

    *file_name*:
      Name of the snapshot file to be written, e.g. VAR0 or var.

    *datadir*:
      Directory where the data is stored.

    *precision*:
      Single 'f' or double 'd' precision.

    *persist*:
      optional dictionary of persistent variable.

    *settings*:
      optional dictionary of persistent variable.

    *param*:
      optional Param object.

    *grid*:
      optional Pencil Grid object of grid parameters.

    *nghost*:
      Number of ghost zones.

    *lghosts*:
      If True the snapshot contains the ghost zones.

    *indx*
      Index object of index for each variable in f-array

    *unit*:
      Optional dictionary of simulation units.

    *quiet*:
      Option to print output.

    *t*:
      Time of the snapshot.

    *xyz*:
      xyz arrays of the domain with ghost zones.
      This will normally be obtained from Grid object, but facility to
      redefine an alternative grid value.

    *lshear*:
      Flag for the shear.

    *driver*
      File driver for hdf5 io for use in serial or MPI parallel

    *comm*
      MPI objects supplied if driver is 'mpio'
    """

    import os.path
    import numpy as np
    import h5py
    from .. import read
    from .. import sim


    #test if simulation directory
    if not sim.is_sim_dir():
        print("ERROR: Directory needs to be a simulation")
    if indx == None:
        indx = read.index()
    #
    if settings == None:
        settings = {}
        skeys = ['l1', 'l2', 'm1', 'm2', 'n1', 'n2',
                 'nx', 'ny', 'nz', 'mx', 'my', 'mz',
                 'nprocx', 'nprocy', 'nprocz',
                 'maux', 'mglobal', 'mvar', 'precision',
                ]
        dim = read.dim()
        for key in skeys:
            settings[key] = dim.__getattribute__(key)
        settings['precision'] = precision.encode()
        settings['nghost'] = nghost
        settings['version'] = np.int32(0)
    nprocs = settings['nprocx']*settings['nprocy']*settings['nprocz']
    gkeys = ['x', 'y', 'z', 'Lx', 'Ly', 'Lz', 'dx', 'dy', 'dz',
             'dx_1', 'dy_1', 'dz_1', 'dx_tilde', 'dy_tilde', 'dz_tilde',
            ]
    if grid == None:
        grid = read.grid(quiet=True)
    else:
        gd_err = False
        for key in gkeys:
            if not key in grid.__dict__.keys():
                print("ERROR: key "+key+" missing from grid")
                gd_err = True
        if gd_err:
            print("ERROR: grid incomplete")
    ukeys = ['length', 'velocity', 'density', 'magnetic', 'time',
             'temperature', 'flux', 'energy', 'mass', 'system',
            ]
    if param == None:
        param = read.param(quiet=True)
        param.__setattr__('unit_mass',param.unit_density*param.unit_length**3)
        param.__setattr__('unit_energy',param.unit_mass*param.unit_velocity**2)
        param.__setattr__('unit_time',param.unit_length/param.unit_velocity)
        param.__setattr__('unit_flux',param.unit_mass/param.unit_time**3)
        param.unit_system = param.unit_system.encode()

    #check whether the snapshot matches the simulation shape
    if lghosts:
        try:
            snapshot.shape[0] == settings['mvar']
            snapshot.shape[1] == settings['mx']
            snapshot.shape[2] == settings['my']
            snapshot.shape[3] == settings['mz']
        except ValueError:
            print("ERROR: snapshot shape {} ".format(snapshot.shape)+
                  "does not match simulation dimensions with ghosts.")
    else:
        try:
            snapshot.shape[0] == settings['mvar']
            snapshot.shape[1] == settings['nx']
            snapshot.shape[2] == settings['ny']
            snapshot.shape[3] == settings['nz']
        except ValueError:
            print("ERROR: snapshot shape {} ".format(snapshot.shape)+
                  "does not match simulation dimensions without ghosts.")

    #Determine the precision used and ensure snapshot has correct data_type.
    if precision == 'f':
        data_type = np.float32
        snapshot = np.float32(snapshot)
    elif precision == 'd':
        data_type = np.float64
        snapshot = np.float64(snapshot)
    else:
        print("ERROR: Precision {0} not understood.".format(precision)+
              " Must be either 'f' or 'd'")
        return -1

    # Check that the shape does not conflict with the proc numbers.
    if (settings['nx']%settings['nprocx'] > 0) or\
       (settings['ny']%settings['nprocy'] > 0) or\
       (settings['nz']%settings['nprocz'] > 0):
        print('ERROR: Shape of the input array is not compatible with the '+
              'cpu layout. Make sure that nproci devides ni.')
        return -1

    # Check the shape of the xyz arrays if specified and overwrite grid values.
    if x != None:
        if len(x) != settings['mx']:
            print('ERROR: x array is incompatible with the shape of snapshot.')
            return -1
        grid.x = x
    if y != None:
        if len(y) != settings['my']:
            print('ERROR: y array is incompatible with the shape of snapshot.')
            return -1
        grid.y = y
    if z != None:
        if len(z) != settings['mz']:
            print('ERROR: z array is incompatible with the shape of snapshot.')
            return -1
        grid.z = z

    # Define a time.
    if t is None:
        t = 0

    # Create the data directory if it doesn't exist.
    if not os.path.exists(datadir):
        try:
            os.mkdir(datadir)
        except FileExistsError:
            pass
    #open file for writing data
    filename = os.path.join(datadir,file_name+'.h5')
    with h5py.File(filename, 'w', driver=driver, comm=comm) as ds:
        # Write the data.
        data_grp = ds.create_group('data')
        for key in indx.__dict__.keys():
            if key in ['uu','keys','aa','KR_Frad','uun','gg']:
                continue
            #create ghost zones if required
            if not lghosts:
                tmp_arr = np.zeros([snapshot.shape[1]+2*nghost,
                                   snapshot.shape[2]+2*nghost,
                                   snapshot.shape[3]+2*nghost])
                tmp_arr[dim.n1:dim.n2+1, dim.m1:dim.m2+1, dim.l1:dim.l2+1
                       ] = np.array(snapshot[indx.__getattribute__(key)])
                data_grp.create_dataset(key,
                                        data=(tmp_arr), dtype=data_type)
            else:
                data_grp.create_dataset(
                    key,
                    data=np.array(snapshot[indx.__getattribute__(key)-1]),
                    dtype=data_type)
        # add time data
        ds.create_dataset('time', data=np.array(t), dtype=data_type)
        # add settings
        sets_grp = ds.create_group('settings')
        for key in settings.keys():
            if 'precision' in key:
                sets_grp.create_dataset(key, data=(settings[key],))
            else:
                sets_grp.create_dataset(key, data=(settings[key]))
        # add grid
        grid_grp = ds.create_group('grid')
        for key in gkeys:
            grid_grp.create_dataset(key, data=(grid.__getattribute__(key)))
        grid_grp.create_dataset('Ox',data=(param.__getattribute__('xyz0')[0]))
        grid_grp.create_dataset('Oy',data=(param.__getattribute__('xyz0')[1]))
        grid_grp.create_dataset('Oz',data=(param.__getattribute__('xyz0')[2]))
        # add physical units
        unit_grp = ds.create_group('unit')
        for key in ukeys:
            if 'system' in key:
                unit_grp.create_dataset(key, data=(param.__getattribute__('unit_'+key),))
            else:
                unit_grp.create_dataset(key, data=param.__getattribute__('unit_'+key))
        # add optional pesristent data
        if persist != None:
            pers_grp = ds.create_group('persist')
            for key in persist.keys():
                if not quiet:
                    print(key,type(persist[key][0]))
                arr = np.empty(nprocs,dtype=type(persist[key][0]))
                arr[:] = persist[key][()]
                pers_grp.create_dataset(key, data=(arr))
#    return 0

def write_h5_grid(file_name='grid', datadir='data', precision='d', nghost=3,
                  settings=None, param=None, grid=None, unit=None, quiet=True,
                  driver=None, comm=None):
    """
    Write the grid information as hdf5.
    We assume by default that a run simulation directory has already been
    constructed, but start has not been executed in h5 format so that
    binary sim files dim, grid and param files are already present in the sim
    directroy, or provided from an old binary sim source directory as inputs.

    call signature:

    write_h5_grid(file_name='grid', datadir='data', precision='d', nghost=3,
                  settings=None, param=None, grid=None, unit=None, quiet=True,
                  driver=None, comm=None)

    Keyword arguments:

    *file_name*:
      Prefix of the file name to be written, 'grid'.

    *datadir*:
      Directory where 'grid.h5' is stored.

    *precision*:
      Single 'f' or double 'd' precision.

    *nghost*:
      Number of ghost zones.

    *settings*:
      Optional dictionary of persistent variable.

    *param*:
      Optional Param object.

    *grid*:
      Optional Pencil Grid object of grid parameters.

    *unit*:
      Optional dictionary of simulation units.

    *quiet*:
      Option to print output.
    """

    import os.path
    import numpy as np
    import h5py
    from .. import read
    from .. import sim


    #test if simulation directory
    if not sim.is_sim_dir():
        print("ERROR: Directory needs to be a simulation")
    #
    if settings == None:
        settings = {}
        skeys = ['l1', 'l2', 'm1', 'm2', 'n1', 'n2',
                 'nx', 'ny', 'nz', 'mx', 'my', 'mz',
                 'nprocx', 'nprocy', 'nprocz',
                 'maux', 'mglobal', 'mvar', 'precision',
                ]
        dim = read.dim()
        for key in skeys:
            settings[key] = dim.__getattribute__(key)
        settings['precision'] = precision.encode()
        settings['nghost'] = nghost
        settings['version'] = np.int32(0)
    gkeys = ['x', 'y', 'z', 'Lx', 'Ly', 'Lz', 'dx', 'dy', 'dz',
             'dx_1', 'dy_1', 'dz_1', 'dx_tilde', 'dy_tilde', 'dz_tilde',
            ]
    if grid == None:
        grid = read.grid(quiet=True)
    else:
        gd_err = False
        for key in gkeys:
            if not key in grid.__dict__.keys():
                print("ERROR: key "+key+" missing from grid")
                gd_err = True
        if gd_err:
            print("ERROR: grid incomplete")
    ukeys = ['length', 'velocity', 'density', 'magnetic', 'time',
             'temperature', 'flux', 'energy', 'mass', 'system',
            ]
    if param == None:
        param = read.param(quiet=True)
        param.__setattr__('unit_mass',param.unit_density*param.unit_length**3)
        param.__setattr__('unit_energy',param.unit_mass*param.unit_velocity**2)
        param.__setattr__('unit_time',param.unit_length/param.unit_velocity)
        param.__setattr__('unit_flux',param.unit_mass/param.unit_time**3)
        param.unit_system = param.unit_system.encode()

    # Create the data directory if it doesn't exist.
    if not os.path.exists(datadir):
        try:
            os.mkdir(datadir)
        except FileExistsError:
            pass
    #open file for writing data
    filename = os.path.join(datadir,file_name+'.h5')
    with h5py.File(filename, 'w', driver=driver, comm=comm) as ds:
        # add settings
        sets_grp = ds.create_group('settings')
        for key in settings.keys():
            if 'precision' in key:
                sets_grp.create_dataset(key, data=(settings[key],))
            else:
                sets_grp.create_dataset(key, data=(settings[key]))
        # add grid
        grid_grp = ds.create_group('grid')
        for key in gkeys:
            grid_grp.create_dataset(key, data=(grid.__getattribute__(key)))
        grid_grp.create_dataset('Ox',data=(param.__getattribute__('xyz0')[0]))
        grid_grp.create_dataset('Oy',data=(param.__getattribute__('xyz0')[1]))
        grid_grp.create_dataset('Oz',data=(param.__getattribute__('xyz0')[2]))
        # add physical units
        unit_grp = ds.create_group('unit')
        for key in ukeys:
            if 'system' in key:
                unit_grp.create_dataset(key, data=(param.__getattribute__('unit_'+key),))
            else:
                unit_grp.create_dataset(key, data=param.__getattribute__('unit_'+key))
#    return 0

def write_h5_averages(aver, file_name='xy', datadir='data/averages',
                   precision='d', indx=None, trange=None, quiet=True,
                   append=False, dim=None, driver=None, comm=None):
    """
    Write an hdf5 format averages dataset given as an Averages object.
    We assume by default that a run simulation directory has already been
    constructed and start completed successfully in h5 format so that
    files dim, grid and param files are already present.
    If not the contents of these will need to be supplied as dictionaries
    along with persist if included.

    call signature:

    write_h5_averages(aver, file_name='xy', datadir='data/averages',
                   precision='d', indx=None, trange=None, quiet=True)

    Keyword arguments:

    *aver*:
      Averages object.
      Must be of shape [n_vars, n1] for averages across 'xy', 'xz' or 'yz'.
      Must be of shape [n_vars, n1, n2] for averages across 'y', 'z'.

    *file_name*:
      Name of the snapshot file to be written, e.g. 'xy', 'xz', 'yz', 'y', 'z'.

    *datadir*:
      Directory where the data is stored.

    *precision*:
      Single 'f' or double 'd' precision.

    *indx*
      Restrict iterative range to be written.

    *trange*:
      Restrict time range to be written.

    *append*
      For large binary files the data may need to be appended iteratively.

    *dim*
      Dim object required if the large binary files are supplied in chunks.
    """

    import os.path
    import numpy as np
    import h5py
    from .. import read
    from .. import sim

    #test if simulation directory
    if not sim.is_sim_dir():
        print("ERROR: Directory needs to be a simulation")
        return -1
    if not os.path.exists(datadir):
        try:
            os.mkdir(datadir)
        except FileExistsError:
            pass
    #open file for writing data
    filename = os.path.join(datadir,file_name+'.h5')
    if append:
        state = 'a'
    else:
        state = 'w'
    #number of iterations to record
    nt=aver.t.shape[0]
    print('saving '+filename)
    with h5py.File(filename, state, driver=driver, comm=comm) as ds:
        for key in aver.__getattribute__(file_name).__dict__.keys():
            data=aver.__getattribute__(file_name).__getattribute__(key)[()] 
            nt=min(nt,data.shape[0])
            for it in range(0,nt):
                if not ds.__contains__(str(it)):
                    ds.create_group(str(it))
                if not ds[str(it)].__contains__('time'):
                    ds[str(it)].create_dataset('time', data=aver.t[it])
                if not ds[str(it)].__contains__(key):
                    ds[str(it)].create_dataset(key,
                                   data=data[it])
        if not ds.__contains__('last'):
            ds.create_dataset('last', data=nt-1)

def write_h5_slices(vslice, coordinates, positions, datadir='data/slices', 
                   precision='d', indx=None, trange=None, quiet=True,
                   append=False, dim=None, driver=None, comm=None):
    """
    Write an hdf5 format slices dataset given as an Slices object.
    We assume by default that a run simulation directory has already been
    constructed and start completed successfully in h5 format so that
    files dim, grid and param files are already present.
    If not the contents of these will need to be supplied as dictionaries
    along with persist if included.

    call signature:

    write_h5_slices(vslice, coordinates, positions,
                   datadir='data/slices', 
                   precision='d', indx=None, trange=None, quiet=True)

    Keyword arguments:

    *vslice*:
      Slices object.
      Object with attributes 't', extensions e.g, 'xy', 'xy2', 'xz', 'yz'
      and data fields of shape [nt, n1, n2] e.g 'uu1', 'uu2', 'uu3', ...

    *coordinates*
      Dictionary of lmn indices of all slices in the object n for 'xy', etc.
      Obtained from 'data/positions.dat' in source binary simulation

    *positions*
      Dictionary of xyz values of all slices in the object z for 'xy', etc.
      Obtained from source binary simulation grid at coordinates.

    *datadir*:
      Directory where the data is stored.

    *precision*:
      Single 'f' or double 'd' precision.

    *indx*
      Restrict iterative range to be written.

    *trange*:
      Restrict time range to be written.

    *append*
      For large binary files the data may need to be appended iteratively.

    *dim*
      Dim object required if the large binary files are supplied in chunks.
    """

    import os.path
    import numpy as np
    import h5py
    from .. import read
    from .. import sim

    #test if simulation directory
    if not sim.is_sim_dir():
        print("ERROR: Directory needs to be a simulation")
        return -1
    if not os.path.exists(datadir):
        try:
            os.mkdir(datadir)
        except FileExistsError:
            pass
    #open file for writing data
    nt = vslice.t.shape[0]
    for extension in vslice.__dict__.keys():
        if not extension in 't':
            for field in vslice.__getattribute__(extension).__dict__.keys():
                filename = os.path.join(datadir,field+'_'+extension+'.h5')
                if append:
                    state = 'a'
                else:
                    state = 'w'
                #number of iterations to record
                print('saving '+filename)
                with h5py.File(filename, state, driver=driver, comm=comm) as ds:
                    for it in range(1,nt+1):
                        if not ds.__contains__(str(it)):
                            ds.create_group(str(it))
                        if not ds[str(it)].__contains__('time'):
                            ds[str(it)].create_dataset('time', data=vslice.t[it-1])
                        if not ds[str(it)].__contains__('data'):
                            ds[str(it)].create_dataset('data',
                                           data=vslice.__getattribute__(extension).__getattribute__(field)[it-1])
                        if not ds[str(it)].__contains__('coordinate'):
                            ds[str(it)].create_dataset('coordinate',
                                           data=(np.int32(coordinates[extension]),))
                        if not ds[str(it)].__contains__('position'):
                            ds[str(it)].create_dataset('position',
                                           data=positions[extension])
                    if not ds.__contains__('last'):
                        ds.create_dataset('last', data=(np.int32(nt),))
