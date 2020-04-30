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

    read(datadir='data', param=None, dim=None)

    Keyword arguments:

    *datadir*:
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


    def read(self, datadir='data', param=None, dim=None):
        """
        Read Pencil Code index data from index.pro.

        call signature:

        read(self, datadir='data', param=None, dim=None)

        Keyword arguments:

        *datadir*:
          Directory where the data is stored.

        *param*
          Parameter object.

        *dim*
          Dimension object.
        """

        import os
        from .. import read

        if param is None:
            param = read.param(datadir=datadir, quiet=True)
        if dim is None:
            dim = read.dim(datadir=datadir)

        if param.lwrite_aux:
            totalvars = dim.mvar + dim.maux
        else:
            totalvars = dim.mvar

        index_file = open(os.path.join(datadir, 'index.pro'))
        ntestfield, ntestflow, ntestlnrho, ntestscalar = 0, 0, 0, 0
        for line in index_file.readlines():
            clean = line.strip()
            name = clean.split('=')[0].strip().replace('[', '').replace(']', '')
            if clean.split('=')[1].strip().startswith('intarr(370)'):
                continue
            val = int(clean.split('=')[1].strip())

            if val != 0  and val <= totalvars \
                and not name.startswith('i_') and name.startswith('i'):
                name = name.lstrip('i')
                if name == 'lnTT' and param.ltemperature_nolog:
                    name = 'tt'
                if name == 'aatest':
                    iaatest = val
                if name == 'uutest':
                    iuutest = val
                if name == 'hhtest':
                    ihhtest = val
                if name == 'cctest':
                    icctest = val
                setattr(self, name, val)
                
            elif name == 'ntestfield':
                ntestfield = val
            elif name == 'ntestflow':
                ntestflow = val
            elif name == 'ntestlnrho':
                ntestlnrho = val
            elif name == 'ntestscalar':
                ntestscalar = val
        if ntestfield > 0:
            self.__delattr__('aatest') 
            for i in range(1,ntestfield+1):
                setattr(self, 'aatest'+str(i), iaatest-1+i)
        if ntestflow > 0:
            self.__delattr__('uutest') 
            for i in range(1,ntestflow+1):
                setattr(self, 'uutest'+str(i), iuutest-1+i)
        if ntestlnrho > 0:
            self.__delattr__('hhtest') 
            for i in range(1,ntestlnrho+1):
                setattr(self, 'hhtest'+str(i), ihhtest-1+i)
        if ntestscalar > 0:
            self.__delattr__('cctest') 
            for i in range(1,ntestscalar+1):
                setattr(self, 'cctest'+str(i), icctest-1+i)
