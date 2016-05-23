# index.py
#
# Read the index.pro file.
#
# Authors: S. Candelaresi (iomsn1@gmail.com), T. Gastine (tgastine@ast.obs-mip.fr).
"""
Contains the index information.
"""


def index(*args, **kwargs):
    """Read index.pro and return an Index object.

    Params:
    ------
    data_dir='data/' (optionnal)
    """
    return Index(*args, **kwargs)


class Index(object):
    """
    Index -- holds pencil code index data.
    """

    def __init__(self):
        """
        Fill members with default values.
        """

        self.keys = []


    def read(self, data_dir='data/', param=None, dim=None):
        """
        Read Pencil Code index data.

        call signature:

        read(self, data_dir='data/', param=None, dim=None)

        Keyword arguments:

        *file_name*:
          Name of the time series file.

        *data_dir*:
          Directory where the data is stored.

        *quiet*
          Flag for switching of output.

        *comment_char*
          Comment character in the time series file.
        """




    """
     read index.pro and returns a read_index class which is composed
     of a python dictionnary

     Author: T. Gastine (tgastine@ast.obs-mip.fr)
    """
    def __init__(self, data_dir='data/', param=None, dim=None):
        """Constructor:
         -----------

         Params:
         ------
            data_dir='data/' (optionnal)

         Returns:
         -------
            a read_index class
        """
        import os
        from pencil import read_param, read_dim
        
        if param is None:
            param = read_param(data_dir=data_dir, quiet=True)
        if dim is None:
            dim = read_dim(data_dir=data_dir)

        if param.lwrite_aux:
            totalvars = dim.mvar+dim.maux
        else:
            totalvars = dim.mvar

        f = open(os.path.join(data_dir,'index.pro'))
        for line in f.readlines():
            clean = line.strip()
            name=clean.split('=')[0].strip().replace('[','').replace(']','')
            if (clean.split('=')[1].strip().startswith('intarr(370)')):
                continue
            val=int(clean.split('=')[1].strip())
            #            print name,val
            # need to compare val to totalvars as global indices
            # may be present in index.pro
            #            if (val != 0 and val <= totalvars and \
            if (val != 0  and val <= totalvars \
                and not name.startswith('i_') and name.startswith('i')):
                name=name.lstrip('i')
                if (name == 'lnTT' and param.ltemperature_nolog):
                    name = 'tt'
                self[name] = val
