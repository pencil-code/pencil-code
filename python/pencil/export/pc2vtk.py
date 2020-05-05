# pc2vtk.py
#
# Convert data from Pencil Code format to vtk. Adopted from pc2vtk.pro.
#
# Author: Simon Candelaresi (iomsn1@googlemail.com).
"""
Contains the routines for converting PencilCode data data into vtk format.
"""

def var2vtk(var_file='var.dat', datadir='data', proc=-1,
            variables=None, b_ext=False, magic=[],
            destination='work', quiet=True, trimall=True, ti=-1, tf=-1):
    """
    Convert data from PencilCode format to vtk.

    call signature::

      var2vtk(var_file='', datadir='data', proc=-1,
             variables='', b_ext=False,
             destination='work', quiet=True, trimall=True, ti=-1, tf=-1)

    Read *var_file* and convert its content into vtk format. Write the result
    in *destination*.

    Keyword arguments:

      *var_file*:
        The original var_file.

      *datadir*:
        Directory where the data is stored.

      *proc*:
        Processor which should be read. Set to -1 for all processors.

      *variables*:
        List of variables which should be written. If None all.

      *b_ext*:
        Add the external magnetic field.

      *destination*:
        Destination file.

      *quiet*:
        Keep quiet when reading the var files.

      *trimall*:
        Trim the data cube to exclude ghost zones.

      *ti, tf*:
        Start and end index for animation. Leave negative for no animation.
        Overwrites variable var_file.
    """

    import numpy as np
    import sys
    from .. import read
    from .. import math

    # Determine of we want an animation.
    if ti < 0 or tf < 0:
        animation = False
    else:
        animation = True

    # If no variables specified collect all by default
    if not variables:
        variables = []
        indx = read.index()
        for key in indx.__dict__.keys(): 
            if 'keys' not in key:
                variables.append(key)
        if 'uu' in variables:
            magic.append('vort')
            variables.append('vort')
        if 'rho' in variables or 'lnrho' in variables:
            if 'ss' in variables:
                magic.append('tt')
                variables.append('tt')
                magic.append('pp')
                variables.append('pp')
        if 'aa' in variables:
            magic.append('bb')
            variables.append('bb')
            magic.append('jj')
            variables.append('jj')
            variables.append('ab')
            variables.append('b_mag')
            variables.append('j_mag')
    else:
        # Convert single variable string into length 1 list of arrays.
        if (len(variables) > 0):
            if (len(variables[0]) == 1):
                variables = [variables]
        if 'tt' in variables:
            magic.append('tt')
        if 'pp' in variables:
            magic.append('pp')
        if 'bb' in variables:
            magic.append('bb')
        if 'jj' in variables:
            magic.append('jj')
        if 'vort' in variables:
            magic.append('vort')
        if 'b_mag' in variables and not 'bb' in magic:
            magic.append('bb')
        if 'j_mag' in variables and not 'jj' in magic:
            magic.append('jj')
        if 'ab' in variables and not 'bb' in magic:
            magic.append('bb')


    for t_idx in range(ti, tf+1):
        if animation:
            var_file = 'VAR' + str(t_idx)

        # Read the PencilCode variables and set the dimensions.
        var = read.var(var_file=var_file, datadir=datadir, proc=proc,
                       magic=magic, trimall=True, quiet=quiet)

        grid = read.grid(datadir=datadir, proc=proc, trim=trimall, quiet=True)

        params = read.param(quiet=True)

        # Add external magnetic field.
        if (b_ext == True):
            B_ext = np.array(params.b_ext)
            var.bb[0, ...] += B_ext[0]
            var.bb[1, ...] += B_ext[1]
            var.bb[2, ...] += B_ext[2]

        dimx = len(grid.x)
        dimy = len(grid.y)
        dimz = len(grid.z)
        dim = dimx*dimy*dimz
        dx = (np.max(grid.x) - np.min(grid.x))/(dimx-1)
        dy = (np.max(grid.y) - np.min(grid.y))/(dimy-1)
        dz = (np.max(grid.z) - np.min(grid.z))/(dimz-1)

        # Write the vtk header.
        if animation:
            fd = open(destination + str(t_idx) + '.vtk', 'wb')
        else:
            fd = open(destination + '.vtk', 'wb')
        fd.write('# vtk DataFile Version 2.0\n'.encode('utf-8'))
        fd.write('VAR files\n'.encode('utf-8'))
        fd.write('BINARY\n'.encode('utf-8'))
        fd.write('DATASET STRUCTURED_POINTS\n'.encode('utf-8'))
        fd.write('DIMENSIONS {0:9} {1:9} {2:9}\n'.format(dimx, dimy, dimz).encode('utf-8'))
        fd.write('ORIGIN {0:8.12} {1:8.12} {2:8.12}\n'.format(grid.x[0], grid.y[0], grid.z[0]).encode('utf-8'))
        fd.write('SPACING {0:8.12} {1:8.12} {2:8.12}\n'.format(dx, dy, dz).encode('utf-8'))
        fd.write('POINT_DATA {0:9}\n'.format(dim).encode('utf-8'))

        # Write the data.
        for v in variables:
            print('Writing {0}.'.format(v))
            # Prepare the data to the correct format.
            if v == 'ab':
                data = math.dot(var.aa, var.bb)
            elif v == 'b_mag':
                data = np.sqrt(math.dot2(var.bb))
            elif v == 'j_mag':
                data = np.sqrt(math.dot2(var.jj))
            else:
                data = getattr(var, v)
            if sys.byteorder == 'little':
                data = data.astype(np.float32).byteswap()
            else:
                data = data.astype(np.float32)
            # Check if we have vectors or scalars.
            if data.ndim == 4:
                data = np.moveaxis(data, 0, 3)
                fd.write('VECTORS {0} float\n'.format(v).encode('utf-8'))
            else:
                fd.write('SCALARS {0} float\n'.format(v).encode('utf-8'))
                fd.write('LOOKUP_TABLE default\n'.encode('utf-8'))
            fd.write(data.tobytes())

        del(var)

        fd.close()



def slices2vtk(field='', extension='', datadir='data', destination='slices', proc=-1):
    """
    Convert slices from PencilCode format to vtk.

    call signature::

      slices2vtk(field='', extension='', datadir='data', destination='slices', proc=-1)

    Read slice files specified by *variables* and convert
    them into vtk format for the specified extensions.
    Write the result in *destination*.
    NB: You need to have called src/read_videofiles.x before using this script.

    Keyword arguments:

      *field*:
        All allowed fields which can be written as slice files, e.g. b2, uu1, lnrho, ...
        See the pencil code manual for more (chapter: "List of parameters for `video.in'").

      *extension*:
        List of slice positions.

      *datadir*:
        Directory where the data is stored.

      *destination*:
        Destination files.

      *proc*:
        Processor which should be read. Set to -1 for all processors.
    """

    import sys
    import numpy as np
    from .. import read

    # Convert single variable string into length 1 list of arrays.
    if (len(field) > 0):
        if (len(field[0]) == 1):
            field = [field]
    if (len(extension) > 0):
        if (len(extension[0]) == 1):
            extension = [extension]

    # Read the grid dimensions.
    grid = read.grid(datadir=datadir, proc=proc, trim=True, quiet=True)

    # Read the dimensions.
    dim = read.dim(datadir=datadir, proc=proc)

    # Read the user given parameters for the slice positions.
    params = read.param(quiet=True)

    # Read the slice file for all specified variables and extensions.
    slices = read.slices(field=field, extension=extension, datadir=datadir, proc=proc)

    # Determine the position of the slices.
    if params.ix != -1:
        x0 = grid.x[params.ix]
    elif params.slice_position == 'm':
        x0 = grid.x[int(len(grid.x)/2)]
    if params.iy != -1:
        y0 = grid.y[params.iy]
    elif params.slice_position == 'm':
        y0 = grid.y[int(len(grid.y)/2)]
    if params.iz != -1:
        z0 = grid.z[params.iz]
    elif params.slice_position == 'm':
        z0 = grid.z[int(len(grid.z)/2)]
    if params.iz2 != -1:
        z02 = grid.z[params.iz]
    elif params.slice_position == 'm':
        z02 = grid.z[int(len(grid.z)/2)]

    for t_idx, t in enumerate(slices.t):
        for ext in extension:
            # Open the destination file for writing.
            fd = open(destination + '_' + ext + '_' + str(t_idx) + '.vtk', 'wb')

            # Write the header.
            fd.write('# vtk DataFile Version 2.0\n'.encode('utf-8'))
            fd.write('slices {0}\n'.format(ext).encode('utf-8'))
            fd.write('BINARY\n'.encode('utf-8'))
            fd.write('DATASET STRUCTURED_POINTS\n'.encode('utf-8'))
            if ext == 'xy':
                fd.write('DIMENSIONS {0:9} {1:9} {2:9}\n'.format(dim.nx, dim.ny, 1).encode('utf-8'))
                fd.write('ORIGIN {0:8.12} {1:8.12} {2:8.12}\n'.format(grid.x[0], grid.y[0], z0).encode('utf-8'))
                fd.write('SPACING {0:8.12} {1:8.12} {2:8.12}\n'.format(grid.dx, grid.dy, 1.).encode('utf-8'))
                dim_p = dim.nx
                dim_q = dim.ny
            if ext == 'xy2':
                fd.write('DIMENSIONS {0:9} {1:9} {2:9}\n'.format(dim.nx, dim.ny, 1).encode('utf-8'))
                fd.write('ORIGIN {0:8.12} {1:8.12} {2:8.12}\n'.format(grid.x[0], grid.y[0], z02).encode('utf-8'))
                fd.write('SPACING {0:8.12} {1:8.12} {2:8.12}\n'.format(grid.dx, grid.dy, 1.).encode('utf-8'))
                dim_p = dim.nx
                dim_q = dim.ny
            if ext == 'xz':
                fd.write('DIMENSIONS {0:9} {1:9} {2:9}\n'.format(dim.nx, 1, dim.nz).encode('utf-8'))
                fd.write('ORIGIN {0:8.12} {1:8.12} {2:8.12}\n'.format(grid.x[0], y0, grid.z[0]).encode('utf-8'))
                fd.write('SPACING {0:8.12} {1:8.12} {2:8.12}\n'.format(grid.dx, 1., grid.dz).encode('utf-8'))
                dim_p = dim.nx
                dim_q = dim.nz
            if ext == 'yz':
                fd.write('DIMENSIONS {0:9} {1:9} {2:9}\n'.format(1, dim.ny, dim.nz).encode('utf-8'))
                fd.write('ORIGIN {0:8.12} {1:8.12} {2:8.12}\n'.format(x0, grid.y[0], grid.z[0]).encode('utf-8'))
                fd.write('SPACING {0:8.12} {1:8.12} {2:8.12}\n'.format(1., grid.dy, grid.dz).encode('utf-8'))
                dim_p = dim.ny
                dim_q = dim.nz
            fd.write('POINT_DATA {0:9}\n'.format(dim_p*dim_q).encode('utf-8'))

            # Write the data.
            for fi in field:
                data = getattr(getattr(slices, ext), fi)
                fd.write(('SCALARS ' + ext + '_' + fi + ' float\n').encode('utf-8'))
                fd.write('LOOKUP_TABLE default\n'.encode('utf-8'))
                if sys.byteorder == 'little':
                    data = data.astype(np.float32).byteswap()
                else:
                    data = data.astype(np.float32)
                fd.write(data[t_idx].tobytes())

            fd.close()



## Convert PencilCode average file to vtk.
#def aver2vtk(varfile = 'xyaverages.dat', datadir = 'data/',
#            destination = 'xyaverages', quiet = 1):
#    """
#    Convert average data from PencilCode format to vtk.
#
#    call signature::
#
#      aver2vtk(varfile = 'xyaverages.dat', datadir = 'data/',
#            destination = 'xyaverages', quiet = 1):
#
#    Read the average file specified in *varfile* and convert the data
#    into vtk format.
#    Write the result in *destination*.
#    
#    Keyword arguments:
#    
#      *varfile*:
#        Name of the average file. This also specifies which dimensions the
#        averages are taken.
#        
#      *datadir*:
#        Directory where the data is stored.
#       
#      *destination*:
#        Destination file.
#               
#    """
#
#    # read the grid dimensions
#    grid = pc.read_grid(datadir = datadir, trim = True, quiet = True)
#    
#    # read the specified average file
#    if varfile[0:2] == 'xy':
#        aver = pc.read_xyaver()
#        line_len = int(np.round(grid.Lz/grid.dz))
#        l0 = grid.z[int((len(grid.z)-line_len)/2)]
#        dl = grid.dz
#    elif varfile[0:2] == 'xz':
#        aver = pc.read_xzaver()
#        line_len = int(np.round(grid.Ly/grid.dy))
#        l0 = grid.y[int((len(grid.y)-line_len)/2)]
#        dl = grid.dy
#    elif varfile[0:2] == 'yz':
#        aver = pc.read_yzaver()
#        line_len = int(np.round(grid.Lx/grid.dx))
#        l0 = grid.x[int((len(grid.x)-line_len)/2)]
#        dl = grid.dx
#    else:
#        print("aver2vtk: ERROR: cannot determine average file\n")
#        print("aver2vtk: The name of the file has to be either xyaver.dat, xzaver.dat or yzaver.dat\n")
#        return -1
#    keys = aver.__dict__.keys()
#    t = aver.t
#    keys.remove('t')
#    
#    # open the destination file
#    fd = open(destination + '.vtk', 'wb')
#    
#    fd.write('# vtk DataFile Version 2.0\n'.encode('utf-8'))
#    fd.write(varfile[0:2] + 'averages\n'.encode('utf-8'))
#    fd.write('BINARY\n'.encode('utf-8'))
#    fd.write('DATASET STRUCTURED_POINTS\n'.encode('utf-8'))
#    fd.write('DIMENSIONS {0:9} {1:9} {2:9}\n'.format(len(t), line_len, 1).encode('utf-8'))
#    fd.write('ORIGIN {0:8.12} {1:8.12} {2:8.12}\n'.format(float(t[0]), l0, 0.).encode('utf-8'))
#    fd.write('SPACING {0:8.12} {1:8.12} {2:8.12}\n'.format(t[1]-t[0], dl, 1.).encode('utf-8'))
#    fd.write('POINT_DATA {0:9}\n'.format(len(t) * line_len))
#                
#    # run through all variables
#    for var in keys:
#        fd.write(('SCALARS ' + var + ' float\n').encode('utf-8'))
#        fd.write('LOOKUP_TABLE default\n'.encode('utf-8'))
#        for j in range(line_len):
#            for i in range(len(t)):
#                    fd.write(struct.pack(">f", aver.__dict__[var][i,j]))
#                                
#    fd.close()
#
#
#
## Convert PencilCode power spectra to vtk.
#def power2vtk(powerfiles = ['power_mag.dat'],
#              datadir = 'data/', destination = 'power', thickness = 1):
#    """
#    Convert power spectra from PencilCode format to vtk.
#
#    call signature::
#    
#      power2vtk(powerfiles = ['power_mag.dat'],
#            datadir = 'data/', destination = 'power.vtk', thickness = 1):
#    
#    Read the power spectra stored in the power*.dat files
#    and convert them into vtk format.
#    Write the result in *destination*.
#    
#    Keyword arguments:
#    
#      *powerfiles*:
#        The files containing the power spectra.
#        
#      *datadir*:
#        Directory where the data is stored.
#       
#      *destination*:
#        Destination file.
#      
#      *thickness*:
#        Dimension in z-direction. Setting it 2 will create n*m*2 dimensional
#        array of data. This is useful in Paraview for visualizing the spectrum
#        in 3 dimensions. Note that this will simply double the amount of data.
#               
#    """
#
#    # this should correct for the case the user types only one variable
#    if (len(powerfiles) > 0):
#        if (len(powerfiles[0]) == 1):
#            powerfiles = [powerfiles]
#            
#    # read the grid dimensions
#    grid = pc.read_grid(datadir = datadir, trim = True, quiet = True)
#    
#    # leave k0 to 1 now, will fix this later
#    k0 = 1.
#    # leave dk to 1 now, will fix this later
#    dk = 1.
#    
#    # open the destination file
#    fd = open(destination + '.vtk', 'wb')
#    
#    # read the first power spectrum
#    t, power = pc.read_power(datadir + powerfiles[0])
#    
#    fd.write('# vtk DataFile Version 2.0\n'.encode('utf-8'))
#    fd.write('power spectra\n'.encode('utf-8'))
#    fd.write('BINARY\n'.encode('utf-8'))
#    fd.write('DATASET STRUCTURED_POINTS\n'.encode('utf-8'))
#    if (thickness == 1):
#        fd.write('DIMENSIONS {0:9} {1:9} {2:9}\n'.format(len(t), power.shape[1], 1).encode('utf-8'))
#    else:
#        fd.write('DIMENSIONS {0:9} {1:9} {2:9}\n'.format(len(t), power.shape[1], 2).encode('utf-8'))
#    fd.write('ORIGIN {0:8.12} {1:8.12} {2:8.12}\n'.format(float(t[0]), k0, 0.).encode('utf-8'))
#    fd.write('SPACING {0:8.12} {1:8.12} {2:8.12}\n'.format(t[1]-t[0], dk, 1.).encode('utf-8'))
#    if (thickness == 1):
#        fd.write('POINT_DATA {0:9}\n'.format(power.shape[0] * power.shape[1]).encode('utf-8'))
#    else:
#        fd.write('POINT_DATA {0:9}\n'.format(power.shape[0] * power.shape[1] * 2).encode('utf-8'))
#    
#    for powfile in powerfiles:
#        # read the power spectrum
#        t, power = pc.read_power(datadir + powfile)
#
#        fd.write(('SCALARS ' + powfile[:-4] + ' float\n').encode('utf-8'))
#        fd.write('LOOKUP_TABLE default\n'.encode('utf-8'))
#        
#        if (thickness == 1):
#            for j in range(power.shape[1]):
#                for i in range(len(t)):
#                        fd.write(struct.pack(">f", power[i,j]))
#        else:
#            for k in [1,2]:
#                for j in range(power.shape[1]):
#                    for i in range(len(t)):
#                            fd.write(struct.pack(">f", power[i,j]))
#                                
#    fd.close()
