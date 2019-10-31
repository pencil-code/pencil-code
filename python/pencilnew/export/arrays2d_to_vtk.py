def arrays2d_to_vtk(folder, filename, arrays, names, gridx, gridy):
    """Convert 2-D array in a VTK file.
    Args:
        - folder:       where to put the vtk
        - filename:     name of the vtk, '.vtk' will be added automatically
        - arrays:       list of 2d array to export
        - names:        list of names associated with arrays
        - gridx/y:      grid as array
    """

    import struct
    import numpy as np
    from os.path import exists, join
    from pencilnew.math import natural_sort
    from pencilnew.io import mkdir

    print("## producing vtk-file from arrays")

    if not (type(arrays) == type(['list']) and type(names) == type(['list'])):
        print('! ERROR: arrays and names must be a list!')

    if not os.path.exists(folder): os.makedirs(folder)        	# check for vtk folder

    dimx = len(gridx); dimy = len(gridy); dim = dimx * dimy
    dx = (np.max(gridx) - np.min(gridx))/(dimx)
    dy = (np.max(gridy) - np.min(gridy))/(dimy)

    dy = dx*dimx/dimy                                          # rescale dimmension to get a square as output 'image'

    ######## do the vtk output
    fd = open(join(folder, filename+'.vtk'), 'wb')
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
                fd.write(struct.pack(">f", array[kk,jj]))

        fd.write('')

    fd.close()
    print("## Done!")
