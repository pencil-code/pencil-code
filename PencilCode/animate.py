#!/usr/bin/python3
# Last Modification: $Id$
#=======================================================================
# animate.py
#
# Facilities for animating the Pencil Code data.
#
# Chao-Chin Yang, 2015-04-22
#=======================================================================
def slices(field, datadir='./data', **kwarg):
    """Animate the video slices.

    Positional Argument:
        field
            Field name.

    Keyword Arguments:
        datadir
            Data directory.
        **kwarg
            Keyword arguments passed to the respective animator.
    """
    # Chao-Chin Yang, 2015-04-22
    from . import read
    # Read the slices.
    t, s = read.slices(field, datadir=datadir)
    # Get the dimensions, parameters, and grid.
    dim = read.dimensions(datadir=datadir)
    par = read.parameters(datadir=datadir)
    grid = read.grid(datadir=datadir, trim=True)
    # Dispatch.
    ndim = (dim.nxgrid > 1) + (dim.nygrid > 1) + (dim.nzgrid > 1)
    if ndim == 1:
        raise NotImplementedError("1D run")
    elif ndim == 2:
        raise NotImplementedError("2D run")
    elif ndim == 3:
        raise NotImplementedError("3D run")
#=======================================================================
###### Local Functions ######
#=======================================================================
