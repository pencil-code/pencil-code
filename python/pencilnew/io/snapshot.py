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
                   nprocx=1, nprocy=1, nprocz=1, precision='f', nghost=3,
                   t=None, x=None, y=None, z=None, lshear=False):
    """
    Write a snapshot given as numpy array.

    call signature:

    write_snapshot(snapshot, file_name='VAR0', datadir='data',
                   nprocx=1, nprocy=1, nprocz=1, precision='f', nghost=3)

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
        print("ERROR: Precision {0} not understood. Must be either 'f' or 'd'".format(precision))

    # Check that the shape does not conflict with the proc numbers.
    if (nx%nprocx > 0) or (ny%nprocy > 0) or (nz%nprocz > 0):
        print('ERROR: Shape of the input array is not compatible with the cpu layout. \
              Make sure that nproci devides ni.')
        return -1

    # Check the shape of the xyz arrays and generate them if they don't exist.
    if not x is None:
        if len(x) != nx:
            print('ERROR: x array is incompatible with the shape of snapshot.')
    else:
        x = np.arange(0, nx)
    if not y is None:
        if len(y) != ny:
            print('ERROR: y array is incompatible with the shape of snapshot.')
    else:
        y = np.arange(0, ny)
    if not z is None:
        if len(z) != nz:
            print('ERROR: z array is incompatible with the shape of snapshot.')
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
    snapshot_ghosts = np.zeros([snapshot.shape[0], nz+2*nghost, ny+2*nghost, nx+2*nghost],
                               dtype=data_type)
    snapshot_ghosts[:, nghost:-nghost, nghost:-nghost, nghost:-nghost] = snapshot

    # Write the data.
    iproc = 0
    for ipz in range(nprocz):
        for ipy in range(nprocy):
            for ipx in range(nprocx):
                destination_file = FortranFile(os.path.join(datadir, 'proc{0}'.format(iproc), file_name), 'w')
                # Construct and write the snapshot for this cpu.
                snapshot_cpu = snapshot_ghosts[:, ipz*int(nz/nprocz):(ipz+1)*int(nz/nprocz)+2*nghost,
                                               ipy*int(ny/nprocy):(ipy+1)*int(ny/nprocy)+2*nghost,
                                               ipx*int(nx/nprocx):(ipx+1)*int(nx/nprocx)+2*nghost]
                destination_file.write_record(snapshot_cpu)
#                snapshot_cpu.tofile(destination_file.format(iproc))
                # Construct and write the meta data for this cpu.
                x_cpu = x_ghost[ipx*int(nx/nprocx):(ipx+1)*int(nx/nprocx)+2*nghost]
                y_cpu = y_ghost[ipy*int(ny/nprocy):(ipy+1)*int(ny/nprocy)+2*nghost]
                z_cpu = z_ghost[ipz*int(nz/nprocz):(ipz+1)*int(nz/nprocz)+2*nghost]
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