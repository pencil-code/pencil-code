# $Id: ts.py 13781 2011-05-04 06:25:24Z iomsn $
#
# Convert data from PencilCode format to vtk. Adopted from pc2vtk.pro.
#
# Author: Simon Candelaresi (iomsn@physto.se, iomsn1@googlemail.com).
# 
#

import pencil as pc
import numpy as np
import struct

def pc2vtk(varfile = 'var.dat', datadir = 'data/', proc = -1,
           variables = ['rho','uu','bb'], magic = ['vort','bb'],
           destination = 'work', quiet = False):
    """
    Convert data from PencilCode format to vtk.

    call signature::
    
      pc2vtk(varfile = 'var.dat', datadir = 'data/', proc = -1,
           variables = ['rho','uu','bb'], magic = ['vort','bb'],
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
      
      *variables* = [ 'rho' , 'lnrho' , 'uu' , 'bb', 'b_mag', 'jj', 'j_mag', 'aa', 'TT', 'lnTT', 'cc', 'lncc' ]
        Variables which should be written.
        
      *magic*: [ 'vort' , 'bb' ]
        Additional variables which should be written.
       
      *destination*:
        Destination file.
    """

    # reading pc variables and setting dimensions
    var = pc.read_var(varfile = varfile, datadir = datadir, proc = proc,
                    magic = magic, trimall = True, quiet = quiet)
                    
    grid = pc.read_grid(datadir = datadir, proc = proc, trim = True, quiet = True)
    
    dimx = len(grid.x)
    dimy = len(grid.y)
    dimz = len(grid.z)
    dim = dimx * dimy * dimz
    dx = (np.max(grid.x) - np.min(grid.x))/(dimx)
    dy = (np.max(grid.y) - np.min(grid.y))/(dimy)
    dz = (np.max(grid.z) - np.min(grid.z))/(dimz)
    
    fd = open(destination + '.vtk', 'wb')
    fd.write('# vtk DataFile Version 2.0\n')
    fd.write('VAR files\n')
    fd.write('BINARY\n')
    fd.write('DATASET STRUCTURED_POINTS\n')
    fd.write('DIMENSIONS {0:9} {1:9} {2:9}\n'.format(dimx, dimy, dimz))
    fd.write('ORIGIN {0:8.12} {1:8.12} {2:8.12}\n'.format(grid.x[0], grid.y[0], grid.z[0]))
    fd.write('SPACING {0:8.12} {1:8.12} {2:8.12}\n'.format(dx, dy, dz))
    fd.write('POINT_DATA {0:9}\n'.format(dim))
    
    # this should correct for the case the user type only one variable
    if (len(variables) > 0):
        if (len(variables[0]) == 1):
            variables = [variables]
    # this should correct for the case the user type only one variable
    if (len(magic) > 0):
        if (len(magic[0]) == 1):
            magic = [magic]
            
    try:
        index = variables.index('rho')
        print 'writing rho'
        fd.write('SCALARS rho float\n')
        fd.write('LOOKUP_TABLE default\n')    
        for k in range(dimz):
            for j in range(dimy):
                for i in range(dimx):
                    fd.write(struct.pack(">f", var.rho[k,j,i]))
    except:
        pass
    
    try:
        index = variables.index('lnrho')
        print 'writing lnrho'
        fd.write('SCALARS lnrho float\n')
        fd.write('LOOKUP_TABLE default\n')    
        for k in range(dimz):
            for j in range(dimy):
                for i in range(dimx):
                    fd.write(struct.pack(">f", var.lnrho[k,j,i]))
    except:
        pass
                
    try:
        index = variables.index('uu')
        print 'writing uu'
        fd.write('VECTORS vfield float\n')
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
        print 'writing bb'
        fd.write('VECTORS bfield float\n')
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
        print 'writing b_mag'
        fd.write('SCALARS b_mag float\n')
        fd.write('LOOKUP_TABLE default\n')        
        for k in range(dimz):
            for j in range(dimy):
                for i in range(dimx):
                    fd.write(struct.pack(">f", b_mag[k,j,i]))
    except:
        pass
    
    try:
        index = variables.index('jj')
        print 'writing jj'
        fd.write('VECTORS jfield float\n')
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
        print 'writing j_mag'
        fd.write('SCALARS j_mag float\n')
        fd.write('LOOKUP_TABLE default\n')    
        for k in range(dimz):
            for j in range(dimy):
                for i in range(dimx):
                    fd.write(struct.pack(">f", j_mag[k,j,i]))
    except:
        pass
    
    try:
        index = variables.index('aa')
        print 'writing aa'
        fd.write('VECTORS afield float\n')
        for k in range(dimz):
            for j in range(dimy):
                for i in range(dimx):
                    fd.write(struct.pack(">f", var.aa[0,k,j,i]))
                    fd.write(struct.pack(">f", var.aa[1,k,j,i]))
                    fd.write(struct.pack(">f", var.aa[2,k,j,i]))
    except:
        pass
    
    try:
        index = variables.index('TT')
        print 'writing TT'
        fd.write('SCALARS TT float\n')
        fd.write('LOOKUP_TABLE default\n')    
        for k in range(dimz):
            for j in range(dimy):
                for i in range(dimx):
                    fd.write(struct.pack(">f", var.TT[k,j,i]))
    except:
        pass

    try:
        index = variables.index('lnTT')
        print 'writing lnTT'
        fd.write('SCALARS lnTT float\n')
        fd.write('LOOKUP_TABLE default\n')    
        for k in range(dimz):
            for j in range(dimy):
                for i in range(dimx):
                    fd.write(struct.pack(">f", var.lnTT[k,j,i]))
    except:
        pass
                    
    try:
        index = variables.index('cc')
        print 'writing cc'
        fd.write('SCALARS cc float\n')
        fd.write('LOOKUP_TABLE default\n')    
        for k in range(dimz):
            for j in range(dimy):
                for i in range(dimx):
                    fd.write(struct.pack(">f", var.cc[k,j,i]))
    except:
        pass

    try:
        index = variables.index('lncc')
        print 'writing lncc'
        fd.write('SCALARS lncc float\n')
        fd.write('LOOKUP_TABLE default\n')    
        for k in range(dimz):
            for j in range(dimy):
                for i in range(dimx):
                    fd.write(struct.pack(">f", var.lncc[k,j,i]))                    
    except:
        pass
    
    fd.close()




def pc2vtk_vid(ti = 0, tf = 1, datadir = 'data/', proc = -1,
           variables = ['rho','uu','bb'], magic = ['vort','bb'],
           destination = 'animation', quiet = True):
    """
    Convert data from PencilCode format to vtk.

    call signature::
    
      pc2vtk(ti = 0, tf = 1, datadir = 'data/', proc = -1,
           variables = ['rho','uu','bb'], magic = ['vort','bb'],
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
      
      *variables* = [ 'rho' , 'lnrho' , 'uu' , 'bb', 'b_mag', 'jj', 'j_mag', 'aa', 'TT', 'lnTT' ]
        Variables which should be written.
        
      *magic*: [ 'vort' , 'bb' ]
        Additional variables which should be written.
       
      *destination*:
        Destination files without '.vtk' extension.        
    """
    
    for i in range(ti,tf+1):
        varfile = 'VAR' + str(i)
        # reading pc variables and setting dimensions
        var = pc.read_var(varfile = varfile, datadir = datadir, proc = proc,
                        magic = magic, trimall = True, quiet = quiet)
                        
        grid = pc.read_grid(datadir = datadir, proc = proc, trim = True, quiet = True)
    
        dimx = len(grid.x)
        dimy = len(grid.y)
        dimz = len(grid.z)
        dim = dimx * dimy * dimz
        dx = (np.max(grid.x) - np.min(grid.x))/(dimx)
        dy = (np.max(grid.y) - np.min(grid.y))/(dimy)
        dz = (np.max(grid.z) - np.min(grid.z))/(dimz)
        
        fd = open(destination + str(i) + '.vtk', 'wb')
        fd.write('# vtk DataFile Version 2.0\n')
        fd.write('density + magnetic field\n')
        fd.write('BINARY\n')        
        fd.write('DATASET STRUCTURED_POINTS\n')
        fd.write('DIMENSIONS {0:9} {1:9} {2:9}\n'.format(dimx, dimy, dimz))
        fd.write('ORIGIN {0:8.12} {1:8.12} {2:8.12}\n'.format(grid.x[0], grid.y[0], grid.z[0]))
        fd.write('SPACING {0:8.12} {1:8.12} {2:8.12}\n'.format(dx, dy, dz))
        fd.write('POINT_DATA {0:9}\n'.format(dim))
        
        # this should correct for the case the user type only one variable
        if (len(variables) > 0):
            if (len(variables[0]) == 1):
                variables = [variables]
        # this should correct for the case the user type only one variable
        if (len(magic) > 0):
            if (len(magic[0]) == 1):
                magic = [magic]
            
        try:
            index = variables.index('rho')
            print 'writing rho'
            fd.write('SCALARS rho float\n')
            fd.write('LOOKUP_TABLE default\n')    
            for k in range(dimz):
                for j in range(dimy):
                    for i in range(dimx):
                        fd.write(struct.pack(">f", var.rho[k,j,i]))
        except:
            pass
        
        try:
            index = variables.index('lnrho')
            print 'writing lnrho'
            fd.write('SCALARS lnrho float\n')
            fd.write('LOOKUP_TABLE default\n')    
            for k in range(dimz):
                for j in range(dimy):
                    for i in range(dimx):
                        fd.write(struct.pack(">f", var.lnrho[k,j,i]))
        except:
            pass
                    
        try:
            index = variables.index('uu')
            print 'writing uu'
            fd.write('VECTORS vfield float\n')
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
            print 'writing bb'
            fd.write('VECTORS bfield float\n')
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
            print 'writing b_mag'
            fd.write('SCALARS b_mag float\n')
            fd.write('LOOKUP_TABLE default\n')        
            for k in range(dimz):
                for j in range(dimy):
                    for i in range(dimx):
                        fd.write(struct.pack(">f", b_mag[k,j,i]))
        except:
            pass
        
        try:
            index = variables.index('jj')
            print 'writing jj'
            fd.write('VECTORS jfield float\n')
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
            print 'writing j_mag'
            fd.write('SCALARS j_mag float\n')
            fd.write('LOOKUP_TABLE default\n')    
            for k in range(dimz):
                for j in range(dimy):
                    for i in range(dimx):
                        fd.write(struct.pack(">f", j_mag[k,j,i]))
        except:
            pass
        
        try:
            index = variables.index('aa')
            print 'writing aa'
            fd.write('VECTORS afield float\n')
            for k in range(dimz):
                for j in range(dimy):
                    for i in range(dimx):
                        fd.write(struct.pack(">f", var.aa[0,k,j,i]))
                        fd.write(struct.pack(">f", var.aa[1,k,j,i]))
                        fd.write(struct.pack(">f", var.aa[2,k,j,i]))
        except:
            pass
        
        try:
            index = variables.index('TT')
            print 'writing TT'
            fd.write('SCALARS TT float\n')
            fd.write('LOOKUP_TABLE default\n')    
            for k in range(dimz):
                for j in range(dimy):
                    for i in range(dimx):
                        fd.write(struct.pack(">f", var.TT[k,j,i]))
        except:
            pass

        try:
            index = variables.index('lnTT')
            print 'writing lnTT'
            fd.write('SCALARS lnTT float\n')
            fd.write('LOOKUP_TABLE default\n')    
            for k in range(dimz):
                for j in range(dimy):
                    for i in range(dimx):
                        fd.write(struct.pack(">f", var.lnTT[k,j,i]))
        except:
            pass

        try:
            index = variables.index('cc')
            print 'writing cc'
            fd.write('SCALARS cc float\n')
            fd.write('LOOKUP_TABLE default\n')    
            for k in range(dimz):
                for j in range(dimy):
                    for i in range(dimx):
                        fd.write(struct.pack(">f", var.cc[k,j,i]))
        except:
            pass

        try:
            index = variables.index('lncc')
            print 'writing lncc'
            fd.write('SCALARS lncc float\n')
            fd.write('LOOKUP_TABLE default\n')    
            for k in range(dimz):
                for j in range(dimy):
                    for i in range(dimx):
                        fd.write(struct.pack(">f", var.lncc[k,j,i]))                    
        except:
            
        pass
        
        fd.close()
