# $Id: $
#
# Convert var to xml format vtk file using PyEVTK library
#
# https://bitbucket.org/pauloh/pyevtk
#
##!/usr/bin/env python2.7

#import sys
#sys.path.append('/home/cmcnally/vc/pencil-hg/python/')
#sys.path.append('/home/cmcnally/pyevtk/lib/python2.7/site-packages/')
import pencil as pc
import numpy as np
import struct
from evtk.vtk import VtkFile, VtkStructuredGrid



# Convert var files into vtk format
def pc2vtkxml(varfile = 'var.dat', datadir = 'data/', proc = -1,
           variables = ['rho','uu','bb'], magic = [],
           destination = 'work', quiet = True):
    """
    Convert data from PencilCode format to XML vtk.
    Write .vts Structured Grid, not Rectilinear Grid as VisIt screws up reading Rectilinear Grid.
    However, this is set to write large grids in VTK XML, which is not yet suported by VisIt anyways. Use ParaView.

    call signature::
    
      pc2xmlvtk(varfile = 'var.dat', datadir = 'data/', proc = -1,
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
      
      *variables* = [ 'rho' , 'lnrho' , 'uu' , 'bb', 'b_mag', 'jj', 'j_mag', 'aa', 'tt', 'lnTT', 'cc', 'lncc', 'ss', 'vort', 'eth' ]
        Variables which should be written.
        
      *magic*: [ 'vort' , 'bb' ]
        Additional variables which should be written.
       
      *destination*:
        Destination file.
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
        index = variables.index('tt')
        magic.append('tt')
    except:
        pass

    # get endian format of the data 
    format = pc.get_format(datadir = datadir)

    # reading pc variables and setting dimensions
    var = pc.read_var(varfile = varfile, datadir = datadir, proc = proc,
                    magic = magic, trimall = True, quiet = quiet, format = format)
                    
    grid = pc.read_grid(datadir = datadir, proc = proc, trim = True, quiet = True, format = format)

    
    dimx = len(grid.x)
    dimy = len(grid.y)
    dimz = len(grid.z)
    dim = dimx * dimy * dimz

    scalardata = {}
    if ('rho' in variables) :
      rho = np.transpose(var.rho.copy()) 
      scalardata['rho'] = rho
    if ('lnrho' in variables) :
      lnrho = np.transpose(var.lnrho.copy()) 
      scalardata['lnrho'] = lnrho
    if ('tt' in variables) :
      tt = np.transpose(var.tt.copy()) 
      scalardata['tt'] = tt
    if ('lntt' in variables) :
      lntt = np.transpose(var.lntt.copy()) 
      scalardata['lntt'] = lntt
    if ('cc' in variables) :
      cc = np.transpose(var.cc.copy()) 
      scalardata['cc'] = cc
    if ('lncc' in variables) :
      lncc = np.transpose(var.lncc.copy()) 
      scalardata['lncc'] = lncc
    if ('ss' in variables) :
      ss = np.transpose(var.ss.copy()) 
      scalardata['ss'] = ss
    if ('eth' in variables) :
      eth = np.transpose(var.eth.copy()) 
      scalardata['eth'] = eth

    vectordata = {}
    if ('uu' in variables) :
      uu1 = np.transpose(var.uu[0,:,:,:].copy()) 
      uu2 = np.transpose(var.uu[1,:,:,:].copy()) 
      uu3 = np.transpose(var.uu[2,:,:,:].copy()) 
      vectordata['uu'] = (uu1,uu2,uu3)
    if ('bb' in variables) :
      bb1 = np.transpose(var.bb[0,:,:,:].copy()) 
      bb2 = np.transpose(var.bb[1,:,:,:].copy()) 
      bb3 = np.transpose(var.bb[2,:,:,:].copy()) 
      vectordata['bb'] = (bb1,bb2,bb3)
    if ('jj' in variables) :
      jj1 = np.transpose(var.jj[0,:,:,:].copy()) 
      jj2 = np.transpose(var.jj[1,:,:,:].copy()) 
      jj3 = np.transpose(var.jj[2,:,:,:].copy()) 
      vectordata['jj'] = (jj1,jj2,jj3)
    if ('aa' in variables) :
      aa1 = np.transpose(var.aa[0,:,:,:].copy()) 
      aa2 = np.transpose(var.aa[1,:,:,:].copy()) 
      aa3 = np.transpose(var.aa[2,:,:,:].copy()) 
      vectordata['aa'] = (aa1,aa2,aa3)
    if ('vort' in variables) :
      vort1 = np.transpose(var.vort[0,:,:,:].copy()) 
      vort2 = np.transpose(var.vort[1,:,:,:].copy()) 
      vort3 = np.transpose(var.vort[2,:,:,:].copy()) 
      vectordata['vort'] = (vort1,vort2,vort3)



    X = np.zeros([dimx,dimy,dimz])
    Y = np.zeros([dimx,dimy,dimz])
    Z = np.zeros([dimx,dimy,dimz])
    for k in range(dimz):
      for j in range(dimy):
        for i in range(dimx):
          X[i,j,k] = grid.x[i] 
          Y[i,j,k] = grid.y[j]
          Z[i,j,k] = grid.z[k] 

    start = (0,0,0)
    end  = (dimx-1, dimy-1, dimz-1)

    time = np.array([var.t])

    w = VtkFile(destination, VtkStructuredGrid,largeFile=True)


    w.openGrid(start = start, end = end)

    #this s for wirting Time in VisIt files. However, when usign large grid Visit does not work anyways.
    #w.openFieldData()
    #w.addTuple('TIME', time.dtype.name,len(time))
    #w.closeFieldData()

    w.openPiece(start = start, end = end)
    w.openElement("Points")
    w.addData("points", (X,Y,Z))
    w.closeElement("Points")

    w.openData("Point", scalars = scalardata.keys(), vectors = vectordata.keys())
    for key in scalardata:
      w.addData(key,scalardata[key])
    for key in vectordata:
      w.addData(key,vectordata[key])
    w.closeData("Point")

    w.closePiece()
    w.closeGrid()

    #w.appendData( time )
    w.appendData( (X,Y,Z) )
    for key in scalardata:
      w.appendData(data = scalardata[key])
    for key in vectordata:
      w.appendData(data = vectordata[key])
    w.save()


