#$Id$
import os
import numpy as np
from .. import read_param, read_dim


def read_index(*args, **kwargs):
    """Read index.pro and return an Index object.

    Params:
    ------
    datadir='data/' (optionnal)
    """
    return Index(*args, **kwargs)


class Index(dict):
    """
     read index.pro and returns a read_index class which is composed
     of a python dictionnary

     Author: T. Gastine (tgastine@ast.obs-mip.fr)
    """
    def __init__(self, datadir='data/', param=None, dim=None, down=False):
        """Constructor:
         -----------

         Params:
         ------
            datadir='data/' (optionnal)

         Returns:
         -------
            a read_index class
        """
        if param is None:
            param = read_param(datadir=datadir, quiet=True)
        if dim is None:
            dim = read_dim(datadir=datadir,down=down)

        if param.lwrite_aux:
            totalvars = dim.mvar+dim.maux
        else:
            totalvars = dim.mvar

        f = open(os.path.join(datadir,'index.pro'))
        for line in f.readlines():
            clean = line.strip()
            name=clean.split('=')[0].strip().replace('[','').replace(']','')
            if (clean.split('=')[1].strip().startswith('intarr(')):
                continue
            if (clean.split('=')[1].strip().startswith('indgen(')):
                val=int(clean.split('=')[1].strip().replace('indgen(','').split(')')[0])
                app=clean.split('=')[1].strip().split('+')
                val=np.arange(val)
                if (len(app) > 1 ):
                    val=val+int(app[1])
                if (all(val != 0)  and all(val <= totalvars) \
                    and not name.startswith('i_') and name.startswith('i')):
                    name=name.lstrip('i')
                    if (name == 'lnTT' and param.ltemperature_nolog):
                        name = 'tt'
                    self[name] = val
            else:
                val=int(clean.split('=')[1].strip())
                if (val != 0  and val <= totalvars \
                    and not name.startswith('i_') and name.startswith('i')):
                        name=name.lstrip('i')
                        if (name == 'lnTT' and param.ltemperature_nolog):
                                name = 'tt'
                        self[name] = val

