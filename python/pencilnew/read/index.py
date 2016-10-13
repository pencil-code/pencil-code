# index.py
#
# Read the index.pro file.
#
# Authors:
# T. Gastine (tgastine@ast.obs-mip.fr)
# S. Candelaresi (iomsn1@gmail.com)
"""
Contains the index information.
"""


def index(*args, **kwargs):
    """
    Read Pencil Code index data from index.pro.

    call signature:

    read(data_dir='data', param=None, dim=None)

    Keyword arguments:

    *data_dir*:
      Directory where the data is stored.

    *param*
      Parameter object.

    *dim*
      Dimension object.
    """

    index_tmp = Index()
    index_tmp.read(*args, **kwargs)
    return index_tmp


class Index(object):
    """
    Index -- holds pencil code index data.
    """

    def __init__(self):
        """
        Fill members with default values.
        """

        self.keys = []


    def read(self, data_dir='data', param=None, dim=None):
        """
        Read Pencil Code index data from index.pro.

        call signature:

        read(self, data_dir='data/', param=None, dim=None)

        Keyword arguments:

        *data_dir*:
          Directory where the data is stored.

        *param*
          Parameter object.

        *dim*
          Dimension object.
        """

        import os
        import pencilnew.read as read

        if param is None:
            param = read.param(data_dir=data_dir, quiet=True)
        if dim is None:
            dim = read.dim(data_dir=data_dir)

        if param.lwrite_aux:
            totalvars = dim.mvar + dim.maux
        else:
            totalvars = dim.mvar

        index_file = open(os.path.join(data_dir, 'index.pro'))
        for line in index_file.readlines():
            clean = line.strip()
            name = clean.split('=')[0].strip().replace('[', '').replace(']', '')
            if (clean.split('=')[1].strip().startswith('intarr(370)')):
                continue
            val = int(clean.split('=')[1].strip())

            if (val != 0  and val <= totalvars \
                and not name.startswith('i_') and name.startswith('i')):
                name = name.lstrip('i')
                if (name == 'lnTT' and param.ltemperature_nolog):
                    name = 'tt'
                setattr(self, name, val)
