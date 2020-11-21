# arrays3d_to_vtk.py
#
# Export 3d arrays into vtk format.
#
# Authors:
# A. Schreiber (aschreiber@mpia.de)
"""
Contains the exporting routine for 3d arrays into vtk.
"""


def arrays3d_to_vtk(folder, filename, arrays, names, gridx, gridy, gridz):
    """
    Convert 3-D array in a VTK file.

    call signature:

    arrays3d_to_vtk(folder, filename, arrays, names, gridx, gridy)

    Args:
        - folder:       where to put the vtk
        - fil_ename:    name of the vtk, '.vtk' will be added automatically
        - arrays:       list of 3d arrays to export
        - names:        list of names associated with arrays
        - gridx/y:      grid as array
    """

    import struct
    import os
    import sys
    import numpy as np

    print("## Producing vtk-file from arrays.")

    if not (isinstance(arrays, list) and isinstance(names, list)):
        print('ERROR: Arrays and names must be a list!')

    # Check if the output folder exists.
    if not os.path.exists(folder):
        os.makedirs(folder)

    dimx = len(gridx)
    dimy = len(gridy)
    dimz = len(gridz)
    dim = dimx*dimy*dimz
    dx = (np.max(gridx) - np.min(gridx))/dimx
    dy = (np.max(gridy) - np.min(gridy))/dimy
    dz = (np.max(gridz) - np.min(gridz))/dimz

    # Rescale dimmension to get a square as output 'image'.
#    dy = dx*dimx/dimy

    # Do the vtk output.
    fd = open(os.path.join(folder, filename+'.vtk'), 'wb')
    fd.write('# vtk DataFile Version 2.0\n'.encode('utf-8'))
    fd.write('PENCIL CODE data\n'.encode('utf-8'))
    fd.write('BINARY\n'.encode('utf-8'))
    fd.write('DATASET STRUCTURED_POINTS\n'.encode('utf-8'))
    fd.write('DIMENSIONS {0:9} {1:9} {2:9}\n'.format(dimx, dimy, dimz).encode('utf-8'))
    fd.write('ORIGIN {0:8.12} {1:8.12} {2:8.12}\n'.format(gridx[0], gridy[0], gridz[0]).encode('utf-8'))
    fd.write('SPACING {0:8.12} {1:8.12} {2:8.12}\n'.format(dx, dy, dz).encode('utf-8'))
    fd.write('POINT_DATA {0:9}\n'.format(dim).encode('utf-8'))

    for name, array in zip(names, arrays):
        print('Writing {0}.'.format(name))
        data = array
        if sys.byteorder == 'little':
            data = data.astype(np.float32).byteswap()
        else:
            data = data.astype(np.float32)
        if data.ndim == 4:
            data = np.moveaxis(data, 0, 3)
            fd.write('VECTORS {0} float\n'.format(name).encode('utf-8'))
        else:
            fd.write('SCALARS {0} float\n'.format(name).encode('utf-8'))
            fd.write('LOOKUP_TABLE default\n'.encode('utf-8'))
        fd.write(data.tobytes())

        # fd.write('SCALARS '+name+' float\n')
        # fd.write('LOOKUP_TABLE default\n')
        # for kk in range(dimz):
        #     for jj in range(dimy):
        #         for ll in range(dimx):
        #             fd.write(struct.pack(">f", array[kk, jj, ll]))

#        fd.write(''.encode('utf-8'))

    fd.close()
    print("Done!")
