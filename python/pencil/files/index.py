#$Id$
import os
from pencil import read_param, read_dim

class read_index(dict):
    """
     read index.pro and returns a read_index class which is composed
     of a python dictionnary

     Author: T. Gastine (tgastine@ast.obs-mip.fr)
    """
    def __init__(self, datadir='data/', param=None, dim=None):
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
            dim = read_dim(datadir=datadir)

        if param.lwrite_aux:
            totalvars = dim.mvar+dim.maux
        else:
            totalvars = dim.mvar

        f = open(os.path.join(datadir,'index.pro'))
        for line in f.readlines():
            clean = line.strip()
            name=clean.split('=')[0].strip()
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

