# Convert power spectra from PencilCode format to vtk. Adopted from pc2vtk.py.
#
# Author: Simon Candelaresi (iomsn@physto.se, iomsn1@googlemail.com).
# 
#

import pencil as pc
import numpy as np
import struct

def power2vtk(powerFiles = ['mag_spec.dat'], destination = 'spectra.vtk', mulz = 2):
    """
    Convert slices from PencilCode format to vtk.

    call signature::
    
      power2vtk(specFiles = ['mag_spec'], destination = 'spectra.vtk')
    
    Read the power files and convert their content into vtk format.
    Write the result in *destination*.
    NB: The z-dimension is used for the time.
    
    Keyword arguments:
    
      *powerFiles*:
        The files containing the spectra. The files have to be in the 'data' directory.
        
      *destination*:
        Destination file.
        
      *mulz*:
        Multiplication of the data in z-direction.
        Use mulz=2 if you want to do fancy animations in Paraview.
        
    """
    
    # open the destination file for writing
    fd = open(destination, 'wb')
    
    # write the header
    fd.write('# vtk DataFile Version 2.0\n')
    fd.write('power spectra\n')
    fd.write('BINARY\n')

    # rad the first power spectrum
    if (len(powerFiles[0]) > 1):        
        pfile = powerFiles[0]
    else:
        pfile = powerFiles        
    t, power = pc.read_power('data/'+pfile)
    dimk = len(power[0,:])
    dimt = len(t)
    dt = t[1]-t[0]
    
    fd.write('DATASET STRUCTURED_POINTS\n')
    fd.write('DIMENSIONS {0:9} {1:9} {2:9}\n'.format(dimk, dimt, mulz))
    fd.write('ORIGIN {0:8.12} {1:8.12} {2:8.12}\n'.format(0.0, 0.0, 0.0))
    fd.write('SPACING {0:8.12} {1:8.12} {2:8.12}\n'.format(1.0, dt, 1.0))
    fd.write('POINT_DATA {0:9}\n'.format(np.size(power)*mulz))

    print 'writing ', pfile[:-4]
    fd.write('SCALARS '+pfile[:-4]+' float\n')
    fd.write('LOOKUP_TABLE default\n')
    for k in range(mulz):
        for j in range(dimt):
            for i in range(dimk):
                fd.write(struct.pack(">f", power[j,i]))
                            
    # run through all power files
    if (len(powerFiles[0]) > 1):        
        for pfile in powerFiles[1:]:        
            t, power = pc.read_power('data/'+pfile)
            print 'writing ', pfile[:-4]
            fd.write('SCALARS '+pfile[:-4]+' float\n')
            fd.write('LOOKUP_TABLE default\n')        
            for k in range(mulz):
                for j in range(dimt):
                    for i in range(dimk):
                        fd.write(struct.pack(">f", power[j,i]))
            
    fd.close()
