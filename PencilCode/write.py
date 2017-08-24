#! /usr/bin/env python3
# Last Modification: $Id$
#=======================================================================
# write.py
#
# Facilities for writing the Pencil Code data in VTK.
#
# Chao-Chin Yang, 2017-08-24
#=======================================================================
def var(**kwarg):
    """Writes one VAR data file in VTK format under allprocs/.

    Keyword Arguments
        **kwarg
            Keywords passed to module read.
    """
    # Author: Chao-Chin Yang
    # Created: 2017-08-24
    # Last Modified: 2017-08-24
    from . import read
    from evtk.hl import gridToVTK
    from Toolbox import get

    # Read the parameters.
    datadir = kwarg.pop("datadir", "./data")
    par = kwarg.pop("par", None)
    if par is None: par = read.parameters(datadir=datadir)

    # Currently only works for rectilinear grid.
    if par.coord_system != "cartesian":
        raise NotImplementedError("Non-rectilinear grid")

    # Read the variable names.
    varnames = read.varname(datadir=datadir)

    # Read the VAR file.
    f = read.var(datadir=datadir, par=par, **kwarg)

    # Organize the data.
    pointData = {}
    for var in varnames:
        pointData[var] = getattr(f, var)

    # Get the file name of the snapshot.
    if "ivar" in kwarg:
        varfile = "VAR{}".format(kwarg["ivar"])
    elif "varfile" in kwarg:
        varfile = kwarg["varfile"]
    else:
        varfile = "var"

    # Write the data.
    print("Writing", varfile, "in VTK under allprocs/...")
    path = datadir + "/allprocs/" + varfile
    gridToVTK(path, f.x, f.y, f.z, pointData=pointData)
    print("Done. ")
#=======================================================================
def var_all(**kwarg):
    """Writes all VAR files in VTK format under allprocs/.

    Keyword Arguments
        **kwarg
            Keywords passed to var().
    """
    # Author: Chao-Chin Yang
    # Created: 2017-08-24
    # Last Modified: 2017-08-24
    from . import read

    # Read the list of files.
    datadir = kwarg.setdefault("datadir", "./data")
    varNlist = read.varname(datadir=datadir, filename="proc0/varN.list")

    # Read the parameters if necessary and pass it forward.
    if kwarg.setdefault("par") is None:
        kwarg["par"] = read.parameters(datadir=datadir)

    # Process each VAR.
    for varfile in varNlist:
        var(varfile=varfile, **kwarg)
