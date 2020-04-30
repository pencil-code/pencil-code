# arrays2d_to_vtk.py
#
# Export 2d arrays into vtk format.
#
# Authors:
# A. Schreiber (aschreiber@mpia.de)
"""
Contains the exporting routine for 2d arrays into vtk.
"""


def arrays2d_to_vtk(folder, filename, arrays, names, gridx, gridy):
    """
    Convert 2-D array in a VTK file.

    call signature:

    arrays2d_to_vtk(folder, filename, arrays, names, gridx, gridy)

    Args:
        - folder:       where to put the vtk
        - fil_ename:    name of the vtk, '.vtk' will be added automatically
        - arrays:       list of 2d array to export
        - names:        list of names associated with arrays
        - gridx/y:      grid as array
    """

    import struct
    import os
    import numpy as np

    print("## Producing vtk-file from arrays.")

    if not (isinstance(arrays, list) and isinstance(names, list)):
        print('ERROR: Arrays and names must be a list!')

    # Check if the output folder exists.
    if not os.path.exists(folder):
        os.makedirs(folder)

    dimx = len(gridx)
    dimy = len(gridy)
    dim = dimx*dimy
    dx = (np.max(gridx) - np.min(gridx))/dimx
    dy = (np.max(gridy) - np.min(gridy))/dimy

    # Rescale dimmension to get a square as output 'image'.
    dy = dx*dimx/dimy

    # Do the vtk output.
    fd = open(os.path.join(folder, filename+'.vtk'), 'wb')
    fd.write('# vtk DataFile Version 2.0\n')
    fd.write('Pencil Code Data\n')
    fd.write('BINARY\n')
    fd.write('DATASET STRUCTURED_POINTS\n')
    fd.write('DIMENSIONS {0:9} {1:9} 1\n'.format(dimx, dimy))
    fd.write('ORIGIN {0:8.12} {1:8.12} 1\n'.format(gridx[0], gridy[0]))
    fd.write('SPACING {0:8.12} {1:8.12} 1\n'.format(dx, dy))
    fd.write('POINT_DATA {0:9}\n'.format(dim))

    for name, array in zip(names, arrays):
        fd.write('SCALARS '+name+' float\n')
        fd.write('LOOKUP_TABLE default\n')
        for kk in range(dimy):
            for jj in range(dimx):
                fd.write(struct.pack(">f", array[kk, jj]))

        fd.write('')

    fd.close()
    print("Done!")
