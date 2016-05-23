# $Id: pc2vtk.py 13781 2011-05-04 06:25:24Z iomsn $
#
# Convert data from PencilCode format to vtk. Adopted from pc2vtk.pro.
#
# Author: Simon Candelaresi (iomsn@physto.se, iomsn1@googlemail.com).

import pencil as pc
import numpy as np
import struct

# Convert var files into vtk format
def pc2vtk(varfile = 'var.dat', datadir = 'data/', proc = -1,
           variables = ['rho','uu','bb'], magic = [], b_ext = False,
           destination = 'work', quiet = True):
    """
    Convert data from PencilCode format to vtk.

    call signature::
    
      pc2vtk(varfile = 'var.dat', datadir = 'data/', proc = -1,
           variables = ['rho','uu','bb'], magic = [],
           destination = 'work.vtk')
    
    Read *varfile* and convert its content into vtk format. Write the result
    in *destination*.
    
    Keyword arguments:
    
      *varfile*:
        The original varfile.
        
      *datadir*:
        Directory where the data is stored.
       
      *proc*:
        Processor which should be read. Set to -1 for all processors.
      
      *variables* = [ 'rho' , 'lnrho' , 'uu' , 'bb', 'b_mag', 'jj', 'j_mag', 'aa', 'ab', 'TT', 'lnTT', 'cc', 'lncc', 'ss', 'vort' ]
        Variables which should be written.
        
      *magic*: [ 'vort' , 'bb' ]
        Additional variables which should be written.
       
      *b_ext*:
        Add the external magnetic field.
        
      *destination*:
        Destination file.
        
      *quiet*:
        Keep quiet when reading the var files.
    """

    # this should correct for the case the user type only one variable
    if (len(magic) > 0):
        if (len(magic[0]) == 1):
            magic = [magic]

    # make sure magic is set when writing 'vort' or 'bb'
    try:
        index = variables.index('vort')
        magic.append('vort')
    except:
        pass      
    try:
        index = variables.index('bb')
        magic.append('bb')
    except:
        pass
    try:
        index = variables.index('b_mag')
        magic.append('bb')
    except:
        pass
    try:
        index = variables.index('jj')
        magic.append('jj')
    except:
        pass
    try:
        index = variables.index('j_mag')
        magic.append('jj')
    except:
        pass

    # reading pc variables and setting dimensions
    var = pc.read_var(varfile = varfile, datadir = datadir, proc = proc,
                    magic = magic, trimall = True, quiet = quiet)
                    
    grid = pc.read_grid(datadir = datadir, proc = proc, trim = True, quiet = True)
    
    params = pc.read_param(param2 = True, quiet = True)
    B_ext = np.array(params.b_ext)    
    # add external magnetic field
    if (b_ext == True):
        var.bb[0,...] += B_ext[0]
        var.bb[1,...] += B_ext[1]
        var.bb[2,...] += B_ext[2]
       
    dimx = len(grid.x)
    dimy = len(grid.y)
    dimz = len(grid.z)
    dim = dimx * dimy * dimz
    dx = (np.max(grid.x) - np.min(grid.x))/(dimx-1)
    dy = (np.max(grid.y) - np.min(grid.y))/(dimy-1)
    dz = (np.max(grid.z) - np.min(grid.z))/(dimz-1)
    
    fd = open(destination + '.vtk', 'wb')
    fd.write('# vtk DataFile Version 2.0\n'.encode('utf-8'))
    fd.write('VAR files\n'.encode('utf-8'))
    fd.write('BINARY\n'.encode('utf-8'))
    fd.write('DATASET STRUCTURED_POINTS\n'.encode('utf-8'))
    fd.write('DIMENSIONS {0:9} {1:9} {2:9}\n'.format(dimx, dimy, dimz).encode('utf-8'))
    fd.write('ORIGIN {0:8.12} {1:8.12} {2:8.12}\n'.format(grid.x[0], grid.y[0], grid.z[0]).encode('utf-8'))
    fd.write('SPACING {0:8.12} {1:8.12} {2:8.12}\n'.format(dx, dy, dz).encode('utf-8'))
    fd.write('POINT_DATA {0:9}\n'.format(dim).encode('utf-8'))
    
    # this should correct for the case the user type only one variable
    if (len(variables) > 0):
        if (len(variables[0]) == 1):
            variables = [variables]
      
    try:
        index = variables.index('rho')
        print('writing rho')
        fd.write('SCALARS rho float\n'.encode('utf-8'))
        fd.write('LOOKUP_TABLE default\n'.encode('utf-8'))
        for k in range(dimz):
            for j in range(dimy):
                for i in range(dimx):
                    fd.write(struct.pack(">f", var.rho[k,j,i]))
    except:
        pass
    
    try:
        index = variables.index('lnrho')
        print('writing lnrho')
        fd.write('SCALARS lnrho float\n'.encode('utf-8'))
        fd.write('LOOKUP_TABLE default\n'.encode('utf-8'))
        for k in range(dimz):
            for j in range(dimy):
                for i in range(dimx):
                    fd.write(struct.pack(">f", var.lnrho[k,j,i]))
    except:
        pass
                
    try:
        index = variables.index('uu')
        print('writing uu')
        fd.write('VECTORS vfield float\n'.encode('utf-8'))
        for k in range(dimz):
            for j in range(dimy):
                for i in range(dimx):
                    fd.write(struct.pack(">f", var.uu[0,k,j,i]))
                    fd.write(struct.pack(">f", var.uu[1,k,j,i]))
                    fd.write(struct.pack(">f", var.uu[2,k,j,i]))
    except:
        pass
    
    try:
        index = variables.index('bb')
        print('writing bb')
        fd.write('VECTORS bfield float\n'.encode('utf-8'))
        for k in range(dimz):
            for j in range(dimy):
                for i in range(dimx):
                    fd.write(struct.pack(">f", var.bb[0,k,j,i]))
                    fd.write(struct.pack(">f", var.bb[1,k,j,i]))
                    fd.write(struct.pack(">f", var.bb[2,k,j,i]))
    except:
        pass
    
    try:
        index = variables.index('b_mag')
        b_mag = np.sqrt(pc.dot2(var.bb))
        print('writing b_mag')
        fd.write('SCALARS b_mag float\n'.encode('utf-8'))
        fd.write('LOOKUP_TABLE default\n'.encode('utf-8'))
        for k in range(dimz):
            for j in range(dimy):
                for i in range(dimx):
                    fd.write(struct.pack(">f", b_mag[k,j,i]))
    except:
        pass
    
    try:
        index = variables.index('jj')
        print('writing jj')
        fd.write('VECTORS jfield float\n'.encode('utf-8'))
        for k in range(dimz):
            for j in range(dimy):
                for i in range(dimx):
                    fd.write(struct.pack(">f", var.jj[0,k,j,i]))
                    fd.write(struct.pack(">f", var.jj[1,k,j,i]))
                    fd.write(struct.pack(">f", var.jj[2,k,j,i]))
    except:
        pass
        
    try:
        index = variables.index('j_mag')
        j_mag = np.sqrt(pc.dot2(var.jj))
        print('writing j_mag')
        fd.write('SCALARS j_mag float\n'.encode('utf-8'))
        fd.write('LOOKUP_TABLE default\n'.encode('utf-8'))
        for k in range(dimz):
            for j in range(dimy):
                for i in range(dimx):
                    fd.write(struct.pack(">f", j_mag[k,j,i]))
    except:
        pass
    
    try:
        index = variables.index('aa')
        print('writing aa')
        fd.write('VECTORS afield float\n'.encode('utf-8'))
        for k in range(dimz):
            for j in range(dimy):
                for i in range(dimx):
                    fd.write(struct.pack(">f", var.aa[0,k,j,i]))
                    fd.write(struct.pack(">f", var.aa[1,k,j,i]))
                    fd.write(struct.pack(">f", var.aa[2,k,j,i]))
    except:
        pass
    
    try:
        index = variables.index('ab')
        ab = pc.dot(var.aa, var.bb)
        print('writing ab')
        fd.write('SCALARS ab float\n'.encode('utf-8'))
        fd.write('LOOKUP_TABLE default\n'.encode('utf-8'))
        for k in range(dimz):
            for j in range(dimy):
                for i in range(dimx):
                    fd.write(struct.pack(">f", ab[k,j,i]))
    except:
        pass
    
    try:
        index = variables.index('TT')
        print('writing TT')
        fd.write('SCALARS TT float\n'.encode('utf-8'))
        fd.write('LOOKUP_TABLE default\n'.encode('utf-8'))
        for k in range(dimz):
            for j in range(dimy):
                for i in range(dimx):
                    fd.write(struct.pack(">f", var.TT[k,j,i]))
    except:
        pass

    try:
        index = variables.index('lnTT')
        print('writing lnTT')
        fd.write('SCALARS lnTT float\n'.encode('utf-8'))
        fd.write('LOOKUP_TABLE default\n'.encode('utf-8'))
        for k in range(dimz):
            for j in range(dimy):
                for i in range(dimx):
                    fd.write(struct.pack(">f", var.lnTT[k,j,i]))
    except:
        pass
                    
    try:
        index = variables.index('cc')
        print('writing cc')
        fd.write('SCALARS cc float\n'.encode('utf-8'))
        fd.write('LOOKUP_TABLE default\n'.encode('utf-8'))
        for k in range(dimz):
            for j in range(dimy):
                for i in range(dimx):
                    fd.write(struct.pack(">f", var.cc[k,j,i]))
    except:
        pass

    try:
        index = variables.index('lncc')
        print('writing lncc')
        fd.write('SCALARS lncc float\n'.encode('utf-8'))
        fd.write('LOOKUP_TABLE default\n'.encode('utf-8'))
        for k in range(dimz):
            for j in range(dimy):
                for i in range(dimx):
                    fd.write(struct.pack(">f", var.lncc[k,j,i]))                    
    except:
        pass
    
    try:
        index = variables.index('ss')
        print('writing ss')
        fd.write('SCALARS ss float\n'.encode('utf-8'))
        fd.write('LOOKUP_TABLE default\n'.encode('utf-8'))
        for k in range(dimz):
            for j in range(dimy):
                for i in range(dimx):
                    fd.write(struct.pack(">f", var.ss[k,j,i]))                    
    except:
        pass
      
    try:
        index = variables.index('vort')
        print('writing vort')
        fd.write('VECTORS vorticity float\n'.encode('utf-8'))
        for k in range(dimz):
            for j in range(dimy):
                for i in range(dimx):
                    fd.write(struct.pack(">f", var.vort[0,k,j,i]))
                    fd.write(struct.pack(">f", var.vort[1,k,j,i]))
                    fd.write(struct.pack(">f", var.vort[2,k,j,i]))
    except:
        pass

    del(var)
    
    fd.close()



# Convert var files into vtk format and make video
def pc2vtk_vid(ti = 0, tf = 1, datadir = 'data/', proc = -1,
           variables = ['rho','uu','bb'], magic = [], b_ext = False,
           destination = 'animation', quiet = True):
    """
    Convert data from PencilCode format to vtk.

    call signature::
    
      pc2vtk(ti = 0, tf = 1, datadir = 'data/', proc = -1,
           variables = ['rho','uu','bb'], magic = [],
           destination = 'animation')
    
    Read *varfile* and convert its content into vtk format. Write the result
    in *destination*.
    
    Keyword arguments:
    
      *ti*:
        Initial time.
        
      *tf*:
        Final time.
        
      *datadir*:
        Directory where the data is stored.
       
      *proc*:
        Processor which should be read. Set to -1 for all processors.
      
      *variables* = [ 'rho' , 'lnrho' , 'uu' , 'bb', 'b_mag', 'jj', 'j_mag', 'aa', 'ab', 'TT', 'lnTT', 'cc', 'lncc', 'ss', 'vort' ]
        Variables which should be written.
        
      *magic*: [ 'vort' , 'bb' ]
        Additional variables which should be written.
       
      *b_ext*:
        Add the external magnetic field.
        
      *destination*:
        Destination files without '.vtk' extension. 
        
      *quiet*:
        Keep quiet when reading the var files.
    """

    # this should correct for the case the user type only one variable
    if (len(variables) > 0):
        if (len(variables[0]) == 1):
            variables = [variables]
    # this should correct for the case the user type only one variable
    if (len(magic) > 0):
        if (len(magic[0]) == 1):
            magic = [magic]
        
    # make sure magic is set when writing 'vort' or 'bb'
    try:
        index = variables.index('vort')
        magic.append('vort')
    except:
        pass      
    try:
        index = variables.index('bb')
        magic.append('bb')
    except:
        pass
    try:
        index = variables.index('b_mag')
        magic.append('bb')
    except:
        pass
    try:
        index = variables.index('jj')
        magic.append('jj')
    except:
        pass
    try:
        index = variables.index('j_mag')
        magic.append('jj')
    except:
        pass

    for i in range(ti,tf+1):
        varfile = 'VAR' + str(i)
        # reading pc variables and setting dimensions
        var = pc.read_var(varfile = varfile, datadir = datadir, proc = proc,
                        magic = magic, trimall = True, quiet = quiet)
                        
        grid = pc.read_grid(datadir = datadir, proc = proc, trim = True, quiet = True)
        
        params = pc.read_param(param2 = True, quiet = True)
        B_ext = np.array(params.b_ext)    
        # add external magnetic field
        if (b_ext == True):
            var.bb[0,...] += B_ext[0]
            var.bb[1,...] += B_ext[1]
            var.bb[2,...] += B_ext[2]
            
        dimx = len(grid.x)
        dimy = len(grid.y)
        dimz = len(grid.z)
        dim = dimx * dimy * dimz
        dx = (np.max(grid.x) - np.min(grid.x))/(dimx-1)
        dy = (np.max(grid.y) - np.min(grid.y))/(dimy-1)
        dz = (np.max(grid.z) - np.min(grid.z))/(dimz-1)
        
        #fd = open(destination + "{0:1.0f}".format(var.t*1e5) + '.vtk', 'wb')
        fd = open(destination + str(i) + '.vtk', 'wb')
        fd.write('# vtk DataFile Version 2.0\n'.encode('utf-8'))
        fd.write('density + magnetic field\n'.encode('utf-8'))
        fd.write('BINARY\n'.encode('utf-8'))
        fd.write('DATASET STRUCTURED_POINTS\n'.encode('utf-8'))
        fd.write('DIMENSIONS {0:9} {1:9} {2:9}\n'.format(dimx, dimy, dimz).encode('utf-8'))
        fd.write('ORIGIN {0:8.12} {1:8.12} {2:8.12}\n'.format(grid.x[0], grid.y[0], grid.z[0]).encode('utf-8'))
        fd.write('SPACING {0:8.12} {1:8.12} {2:8.12}\n'.format(dx, dy, dz).encode('utf-8'))
        fd.write('POINT_DATA {0:9}\n'.format(dim).encode('utf-8'))
              
        try:
            index = variables.index('rho')
            print('writing rho')
            fd.write('SCALARS rho float\n'.encode('utf-8'))
            fd.write('LOOKUP_TABLE default\n'.encode('utf-8'))
            for k in range(dimz):
                for j in range(dimy):
                    for i in range(dimx):
                        fd.write(struct.pack(">f", var.rho[k,j,i]))
        except:
            pass
        
        try:
            index = variables.index('lnrho')
            print('writing lnrho')
            fd.write('SCALARS lnrho float\n'.encode('utf-8'))
            fd.write('LOOKUP_TABLE default\n'.encode('utf-8'))
            for k in range(dimz):
                for j in range(dimy):
                    for i in range(dimx):
                        fd.write(struct.pack(">f", var.lnrho[k,j,i]))
        except:
            pass
                    
        try:
            index = variables.index('uu')
            print('writing uu')
            fd.write('VECTORS vfield float\n'.encode('utf-8'))
            for k in range(dimz):
                for j in range(dimy):
                    for i in range(dimx):
                        fd.write(struct.pack(">f", var.uu[0,k,j,i]))
                        fd.write(struct.pack(">f", var.uu[1,k,j,i]))
                        fd.write(struct.pack(">f", var.uu[2,k,j,i]))
        except:
            pass
        
        try:
            index = variables.index('bb')
            print('writing bb')
            fd.write('VECTORS bfield float\n'.encode('utf-8'))
            for k in range(dimz):
                for j in range(dimy):
                    for i in range(dimx):
                        fd.write(struct.pack(">f", var.bb[0,k,j,i]))
                        fd.write(struct.pack(">f", var.bb[1,k,j,i]))
                        fd.write(struct.pack(">f", var.bb[2,k,j,i]))
        except:
            pass
        
        try:
            index = variables.index('b_mag')
            b_mag = np.sqrt(pc.dot2(var.bb))
            print('writing b_mag')
            fd.write('SCALARS b_mag float\n'.encode('utf-8'))
            fd.write('LOOKUP_TABLE default\n'.encode('utf-8'))
            for k in range(dimz):
                for j in range(dimy):
                    for i in range(dimx):
                        fd.write(struct.pack(">f", b_mag[k,j,i]))
        except:
            pass
        
        try:
            index = variables.index('jj')
            print('writing jj')
            fd.write('VECTORS jfield float\n'.encode('utf-8'))
            for k in range(dimz):
                for j in range(dimy):
                    for i in range(dimx):
                        fd.write(struct.pack(">f", var.jj[0,k,j,i]))
                        fd.write(struct.pack(">f", var.jj[1,k,j,i]))
                        fd.write(struct.pack(">f", var.jj[2,k,j,i]))
        except:
            pass
        
        try:
            index = variables.index('j_mag')
            j_mag = np.sqrt(pc.dot2(var.jj))
            print('writing j_mag')
            fd.write('SCALARS j_mag float\n'.encode('utf-8'))
            fd.write('LOOKUP_TABLE default\n'.encode('utf-8'))
            for k in range(dimz):
                for j in range(dimy):
                    for i in range(dimx):
                        fd.write(struct.pack(">f", j_mag[k,j,i]))
        except:
            pass
        
        try:
            index = variables.index('aa')
            print('writing aa')
            fd.write('VECTORS afield float\n'.encode('utf-8'))
            for k in range(dimz):
                for j in range(dimy):
                    for i in range(dimx):
                        fd.write(struct.pack(">f", var.aa[0,k,j,i]))
                        fd.write(struct.pack(">f", var.aa[1,k,j,i]))
                        fd.write(struct.pack(">f", var.aa[2,k,j,i]))
        except:
            pass
        
        try:
            index = variables.index('ab')
            ab = pc.dot(var.aa, var.bb)
            print('writing ab')
            fd.write('SCALARS ab float\n'.encode('utf-8'))
            fd.write('LOOKUP_TABLE default\n'.encode('utf-8'))
            for k in range(dimz):
                for j in range(dimy):
                    for i in range(dimx):
                        fd.write(struct.pack(">f", ab[k,j,i]))
        except:
            pass
        
        try:
            index = variables.index('TT')
            print('writing TT')
            fd.write('SCALARS TT float\n'.encode('utf-8'))
            fd.write('LOOKUP_TABLE default\n'.encode('utf-8'))
            for k in range(dimz):
                for j in range(dimy):
                    for i in range(dimx):
                        fd.write(struct.pack(">f", var.TT[k,j,i]))
        except:
            pass

        try:
            index = variables.index('lnTT')
            print('writing lnTT')
            fd.write('SCALARS lnTT float\n'.encode('utf-8'))
            fd.write('LOOKUP_TABLE default\n'.encode('utf-8'))
            for k in range(dimz):
                for j in range(dimy):
                    for i in range(dimx):
                        fd.write(struct.pack(">f", var.lnTT[k,j,i]))
        except:
            pass

        try:
            index = variables.index('cc')
            print('writing cc')
            fd.write('SCALARS cc float\n'.encode('utf-8'))
            fd.write('LOOKUP_TABLE default\n'.encode('utf-8'))
            for k in range(dimz):
                for j in range(dimy):
                    for i in range(dimx):
                        fd.write(struct.pack(">f", var.cc[k,j,i]))
        except:
            pass

        try:
            index = variables.index('lncc')
            print('writing lncc')
            fd.write('SCALARS lncc float\n'.encode('utf-8'))
            fd.write('LOOKUP_TABLE default\n'.encode('utf-8'))
            for k in range(dimz):
                for j in range(dimy):
                    for i in range(dimx):
                        fd.write(struct.pack(">f", var.lncc[k,j,i]))                    
        except:
            pass
        
        try:
            index = variables.index('ss')
            print('writing ss')
            fd.write('SCALARS ss float\n'.encode('utf-8'))
            fd.write('LOOKUP_TABLE default\n'.encode('utf-8'))
            for k in range(dimz):
                for j in range(dimy):
                    for i in range(dimx):
                        fd.write(struct.pack(">f", var.ss[k,j,i]))                    
        except:
            pass

        try:
            index = variables.index('vort')
            print('writing vort')
            fd.write('VECTORS vorticity float\n'.encode('utf-8'))
            for k in range(dimz):
                for j in range(dimy):
                    for i in range(dimx):
                        fd.write(struct.pack(">f", var.vort[0,k,j,i]))
                        fd.write(struct.pack(">f", var.vort[1,k,j,i]))
                        fd.write(struct.pack(">f", var.vort[2,k,j,i]))
        except:
            pass
      
        del(var)
        
        fd.close()



# Convert PencilCode slices to vtk.
def slices2vtk(variables = ['rho'], extensions = ['xy', 'xy2', 'xz', 'yz'],
           datadir = 'data/', destination = 'slices', proc = -1,
           format = 'native'):
    """
    Convert slices from PencilCode format to vtk.

    call signature::
    
      slices2vtk(variables = ['rho'], extensions = ['xy', 'xy2', 'xz', 'yz'],
           datadir = 'data/', destination = 'slices', proc = -1,
           format = 'native'):
    
    Read slice files specified by *variables* and convert
    them into vtk format for the specified extensions.
    Write the result in *destination*.
    NB: You need to have called src/read_videofiles.x before using this script.
    
    Keyword arguments:
    
      *variables*:
        All allowed fields which can be written as slice files, e.g. b2, uu1, lnrho, ...
        See the pencil code manual for more (chapter: "List of parameters for `video.in'").
        
      *extensions*:
        List of slice positions.
      
      *datadir*:
        Directory where the data is stored.
       
      *destination*:
        Destination files.
        
      *proc*:
        Processor which should be read. Set to -1 for all processors.
      
      *format*:
        Endian, one of little, big, or native (default)
       
    """

    # this should correct for the case the user types only one variable
    if (len(variables) > 0):
        if (len(variables[0]) == 1):
            variables = [variables]
    # this should correct for the case the user types only one extension
    if (len(extensions) > 0):
        if (len(extensions[0]) == 1):
            extensions = [extensions]

    # read the grid dimensions
    grid = pc.read_grid(datadir = datadir, proc = proc, trim = True, quiet = True)
    
    # read the user given parameters for the slice positions
    params = pc.read_param(param2 = True, quiet = True)
    
    # run through all specified variables
    for field in variables:
        # run through all specified extensions
        for ext in extensions:
            print("read " + field + ' ' + ext)
            slices, t = pc.read_slices(field=field, datadir=datadir, proc=proc, extension=ext, format=format)
            
            dim_p = slices.shape[2]
            dim_q = slices.shape[1]
            if ext[0] == 'x':
                d_p = (np.max(grid.x) - np.min(grid.x))/(dim_p)                
            else:
                d_p = (np.max(grid.y) - np.min(grid.y))/(dim_p)
            if ext[1] == 'y':
                d_q = (np.max(grid.y) - np.min(grid.y))/(dim_q)
            else:
                d_q = (np.max(grid.z) - np.min(grid.z))/(dim_q)

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

            for i in range(slices.shape[0]):
                # open the destination file for writing
                fd = open(destination + '_' + field + '_' + ext + '_' + str(i) + '.vtk', 'wb')

                # write the header
                fd.write('# vtk DataFile Version 2.0\n')
                fd.write(field + '_' + ext + '\n')
                fd.write('BINARY\n') 
                fd.write('DATASET STRUCTURED_POINTS\n')
                if ext[0:2] == 'xy':
                    x0 = grid.x[0]; y0 = grid.y[0]
                    fd.write('DIMENSIONS {0:9} {1:9} {2:9}\n'.format(dim_p, dim_q, 1))
                    fd.write('ORIGIN {0:8.12} {1:8.12} {2:8.12}\n'.format(x0, y0, z0))
                    fd.write('SPACING {0:8.12} {1:8.12} {2:8.12}\n'.format(grid.dx, grid.dy, 1.))
                elif ext[0:2] == 'xz':
                    x0 = grid.x[0]; z0 = grid.z[0]
                    fd.write('DIMENSIONS {0:9} {1:9} {2:9}\n'.format(dim_p, 1, dim_q))
                    fd.write('ORIGIN {0:8.12} {1:8.12} {2:8.12}\n'.format(x0, y0, z0))
                    fd.write('SPACING {0:8.12} {1:8.12} {2:8.12}\n'.format(grid.dx, 1., grid.dy))
                elif ext[0:2] == 'yz':
                    y0 = grid.y[0]; z0 = grid.z[0]
                    fd.write('DIMENSIONS {0:9} {1:9} {2:9}\n'.format(1, dim_p, dim_q))
                    fd.write('ORIGIN {0:8.12} {1:8.12} {2:8.12}\n'.format(x0, y0, z0))
                    fd.write('SPACING {0:8.12} {1:8.12} {2:8.12}\n'.format(1., grid.dy, grid.dy))                                        
                fd.write('POINT_DATA {0:9}\n'.format(dim_p * dim_q))
        
                fd.write('SCALARS ' + field + '_' + ext + ' float\n')
                fd.write('LOOKUP_TABLE default\n')        
                for j in range(dim_q):
                    for k in range(dim_p):
                            fd.write(struct.pack(">f", slices[i,j,k]))
                                
                fd.close()



# Convert PencilCode average file to vtk.
def aver2vtk(varfile = 'xyaverages.dat', datadir = 'data/',
            destination = 'xyaverages', quiet = 1):
    """
    Convert average data from PencilCode format to vtk.

    call signature::
    
      aver2vtk(varfile = 'xyaverages.dat', datadir = 'data/',
            destination = 'xyaverages', quiet = 1):

    Read the average file specified in *varfile* and convert the data
    into vtk format.
    Write the result in *destination*.
    
    Keyword arguments:
    
      *varfile*:
        Name of the average file. This also specifies which dimensions the
        averages are taken.
        
      *datadir*:
        Directory where the data is stored.
       
      *destination*:
        Destination file.
               
    """

    # read the grid dimensions
    grid = pc.read_grid(datadir = datadir, trim = True, quiet = True)
    
    # read the specified average file
    if varfile[0:2] == 'xy':
        aver = pc.read_xyaver()
        line_len = int(np.round(grid.Lz/grid.dz))
        l0 = grid.z[int((len(grid.z)-line_len)/2)]
        dl = grid.dz
    elif varfile[0:2] == 'xz':
        aver = pc.read_xzaver()
        line_len = int(np.round(grid.Ly/grid.dy))
        l0 = grid.y[int((len(grid.y)-line_len)/2)]
        dl = grid.dy
    elif varfile[0:2] == 'yz':
        aver = pc.read_yzaver()
        line_len = int(np.round(grid.Lx/grid.dx))
        l0 = grid.x[int((len(grid.x)-line_len)/2)]
        dl = grid.dx
    else:
        print("aver2vtk: ERROR: cannot determine average file\n")
        print("aver2vtk: The name of the file has to be either xyaver.dat, xzaver.dat or yzaver.dat\n")
        return -1
    keys = aver.__dict__.keys()
    t = aver.t
    keys.remove('t')
    
    # open the destination file
    fd = open(destination + '.vtk', 'wb')
    
    fd.write('# vtk DataFile Version 2.0\n'.encode('utf-8'))
    fd.write(varfile[0:2] + 'averages\n'.encode('utf-8'))
    fd.write('BINARY\n'.encode('utf-8'))
    fd.write('DATASET STRUCTURED_POINTS\n'.encode('utf-8'))
    fd.write('DIMENSIONS {0:9} {1:9} {2:9}\n'.format(len(t), line_len, 1).encode('utf-8'))
    fd.write('ORIGIN {0:8.12} {1:8.12} {2:8.12}\n'.format(float(t[0]), l0, 0.).encode('utf-8'))
    fd.write('SPACING {0:8.12} {1:8.12} {2:8.12}\n'.format(t[1]-t[0], dl, 1.).encode('utf-8'))
    fd.write('POINT_DATA {0:9}\n'.format(len(t) * line_len))
                
    # run through all variables
    for var in keys:
        fd.write(('SCALARS ' + var + ' float\n').encode('utf-8'))
        fd.write('LOOKUP_TABLE default\n'.encode('utf-8'))
        for j in range(line_len):
            for i in range(len(t)):
                    fd.write(struct.pack(">f", aver.__dict__[var][i,j]))
                                
    fd.close()



# Convert PencilCode power spectra to vtk.
def power2vtk(powerfiles = ['power_mag.dat'],
            datadir = 'data/', destination = 'power', thickness = 1):
    """
    Convert power spectra from PencilCode format to vtk.

    call signature::
    
      power2vtk(powerfiles = ['power_mag.dat'],
            datadir = 'data/', destination = 'power.vtk', thickness = 1):
    
    Read the power spectra stored in the power*.dat files
    and convert them into vtk format.
    Write the result in *destination*.
    
    Keyword arguments:
    
      *powerfiles*:
        The files containing the power spectra.
        
      *datadir*:
        Directory where the data is stored.
       
      *destination*:
        Destination file.
      
      *thickness*:
        Dimension in z-direction. Setting it 2 will create n*m*2 dimensional
        array of data. This is useful in Paraview for visualizing the spectrum
        in 3 dimensions. Note that this will simply double the amount of data.
               
    """

    # this should correct for the case the user types only one variable
    if (len(powerfiles) > 0):
        if (len(powerfiles[0]) == 1):
            powerfiles = [powerfiles]
            
    # read the grid dimensions
    grid = pc.read_grid(datadir = datadir, trim = True, quiet = True)
    
    # leave k0 to 1 now, will fix this later
    k0 = 1.
    # leave dk to 1 now, will fix this later
    dk = 1.
    
    # open the destination file
    fd = open(destination + '.vtk', 'wb')
    
    # read the first power spectrum
    t, power = pc.read_power(datadir + powerfiles[0])
    
    fd.write('# vtk DataFile Version 2.0\n'.encode('utf-8'))
    fd.write('power spectra\n'.encode('utf-8'))
    fd.write('BINARY\n'.encode('utf-8'))
    fd.write('DATASET STRUCTURED_POINTS\n'.encode('utf-8'))
    if (thickness == 1):
        fd.write('DIMENSIONS {0:9} {1:9} {2:9}\n'.format(len(t), power.shape[1], 1).encode('utf-8'))
    else:
        fd.write('DIMENSIONS {0:9} {1:9} {2:9}\n'.format(len(t), power.shape[1], 2).encode('utf-8'))
    fd.write('ORIGIN {0:8.12} {1:8.12} {2:8.12}\n'.format(float(t[0]), k0, 0.).encode('utf-8'))
    fd.write('SPACING {0:8.12} {1:8.12} {2:8.12}\n'.format(t[1]-t[0], dk, 1.).encode('utf-8'))
    if (thickness == 1):
        fd.write('POINT_DATA {0:9}\n'.format(power.shape[0] * power.shape[1]).encode('utf-8'))
    else:
        fd.write('POINT_DATA {0:9}\n'.format(power.shape[0] * power.shape[1] * 2).encode('utf-8'))
    
    for powfile in powerfiles:
        # read the power spectrum
        t, power = pc.read_power(datadir + powfile)

        fd.write(('SCALARS ' + powfile[:-4] + ' float\n').encode('utf-8'))
        fd.write('LOOKUP_TABLE default\n'.encode('utf-8'))
        
        if (thickness == 1):
            for j in range(power.shape[1]):
                for i in range(len(t)):
                        fd.write(struct.pack(">f", power[i,j]))
        else:
            for k in [1,2]:
                for j in range(power.shape[1]):
                    for i in range(len(t)):
                            fd.write(struct.pack(">f", power[i,j]))
                                
    fd.close()
