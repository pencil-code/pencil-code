def array2d_to_vtk(folder, filename, title, array, gridx, gridy):
    """Convert 2-D array in a VTK file.
    Args:
        - folder:       where to put the vtk
        - filename:     name of the vtk, '.vtk' will be added automatically
        - array:        2d array to export
        - gridx/y:      grid as array
    """

    import pencil as pc
    import numpy as np
    import struct
    import glob
    import os
    import sys
    import pen.intern.natural_sort


    print "## producing vtk-file from array"

    if not os.path.exists(folder): os.makedirs(folder)        	# check for vtk folder

    dimx = len(gridx); dimy = len(gridy); dim = dimx * dimy
    dx = (np.max(gridx) - np.min(gridx))/(dimx);   dy = (np.max(gridy) - np.min(gridy))/(dimy)

    dy = dx*dimx/dimy                                          # rescale dimmension to get a square as output 'image'

    ######## do the vtk output
    fd = open(folder + filename + '.vtk', 'wb')
    fd.write('# vtk DataFile Version 2.0\n')
    fd.write('Pencil Code Data from '+title+'\n')
    fd.write('BINARY\n')
    fd.write('DATASET STRUCTURED_POINTS\n')
    fd.write('DIMENSIONS {0:9} {1:9} 1\n'.format(dimx, dimy))
    fd.write('ORIGIN {0:8.12} {1:8.12} 1\n'.format(gridx[0], gridy[0]))
    fd.write('SPACING {0:8.12} {1:8.12} 1\n'.format(dx, dy))
    fd.write('POINT_DATA {0:9}\n'.format(dim))
    ##
    fd.write('SCALARS '+title+' float\n')
    fd.write('LOOKUP_TABLE default\n')
    for kk in range(dimy):
        for jj in range(dimx):
            fd.write(struct.pack(">f", array[kk,jj]))

    fd.write('')

    fd.close()
    print "## Done!"
