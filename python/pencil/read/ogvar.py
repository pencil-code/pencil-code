# ogvar.py
#
# Read the ogvar files.
# This routine inherents the class properties from DataCube in
# var.py, used to read var files.
# NB: the f array returned is C-ordered: f[nvar, nz, ny, nx]
#     NOT Fortran as in Pencil (& IDL):  f[nx, ny, nz, nvar]
#
# Author:
# J. Aarnes (jorgenaarnes@gmail.com)
"""
Contains the read class for the OGVAR file reading,
some simulation attributes and the data cube.
"""


from pencil.read.varfile import DataCube


def ogvar(*args, **kwargs):
    """
    Read OGVAR files from Pencil Code. If proc < 0, then load all data
    and assemble, otherwise load OGVAR file from specified processor.

    The file format written by output() (and used, e.g. in ogvar.dat)
    consists of the followinig Fortran records:
    1. data(mx, my, mz, nvar)
    2. t(1), x(mx), y(my), z(mz), dx(1), dy(1), dz(1), deltay(1)
    Here nvar denotes the number of slots, i.e. 1 for one scalar field, 3
    for one vector field, 8 for ogvar.dat in the case of MHD with entropy.
    but, deltay(1) is only there if lshear is on! need to know parameters.

    call signature:

    ogvar(var_file='', datadir='data', proc=-1, iogvar=-1,
          quiet=True, trimall=False,
          magic=None, sim=None, precision='f')

    Keyword arguments:
        var_file:   Name of the OGVAR file.
        datadir:    Directory where the data is stored.
        proc:       Processor to be read. If -1 read all and assemble to one array.
        iogvar:     Index of the OGVAR file, if var_file is not specified.
        quiet:      Flag for switching off output.
        trimall:    Trim the data cube to exclude ghost zones.
        magic:      Values to be computed from the data, e.g. B = curl(A).
        sim:        Simulation sim object.
        precision:  Precision of the data. Either float 'f' or double 'd'.
    """

    from pencil.sim import __Simulation__

    started = None

    for a in args:
        if type(a, __Simulation__):
            started = a.started()
            break
    else:
        if "sim" in kwargs.keys():
            started = kwargs["sim"].started()
        elif "datadir" in kwargs.keys():
            from os.path import join, exists

            if exists(join(kwargs["datadir"], "time_series.dat")):
                started = True

    if "ivar" not in kwargs.keys():
        ivar = 0
    if not started and ivar != 0:
        print("!! ERROR: Simulation has not jet started. There are not ogvar files.")
        return False

    if "var_file" in kwargs:
        if isinstance(kwargs["var_file"], __Simulation__):
            kwargs["var_file"] = "ogvar.dat"
    else:
        if "varfile" in kwargs:
            kwargs["var_file"] = kwargs["varfile"]
        elif "ivar" in kwargs:
            if kwargs["ivar"] < 0:
                kwargs["var_file"] = "ogvar.dat"
            else:
                kwargs["var_file"] = "OGVAR" + str(kwargs["ivar"])
        else:
            kwargs["var_file"] = "ogvar.dat"

    if kwargs["var_file"][0:2].lower() != "og":
        print("!! ERROR: Read procedure not called with ogvar-file.")
        print("          Did you mean to call read.var() instead?")
        return False

    ogvar_tmp = ogDataCube()
    ogvar_tmp.read(*args, **kwargs)
    if "trim_all" in kwargs:
        trim_all = kwargs["trim_all"]
    else:
        trim_all = True

    ogvar_tmp.transform(trim_all)
    return ogvar_tmp


class ogDataCube(DataCube):
    """
    DataCube -- holds Pencil Code OGVAR file data.
    """

    def __init__(self):
        """
        Initialize the og data cube as the var data cube.
        """

        super(ogDataCube, self).__init__()

    def keys(self):
        for i in self.__dict__.keys():
            print(i)

    def transform(self, trim_all):
        """
        Transform velocity coordinates from ur and uth to ux and uy.
        """

        from pencil.math.transform import pospolar2cart, velpolar2cart

        if trim_all:
            zcoord = 0
        else:
            zcoord = 3

        self.r = self.x
        self.th = self.y
        self.x, self.y = pospolar2cart(self.r, self.th)
        self.ur = self.ux
        self.uth = self.uy
        self.ux, self.uy = velpolar2cart(self.ur, self.uth, self.r, self.th, zcoord)
