#$Id$

from param import read_param

class read_index(dict):
    """
     read index.pro and returns a read_index class which is composed
     of a python dictionnary

     Author: T. Gastine (tgastine@ast.obs-mip.fr)
    """
    def __init__(self, datadir='data/', param=None):
        """Constructor:
         -----------

         Params:
         ------
            datadir='data/' (optionnal)
         
         Returns:
         -------      
            a read_index class
        """
        if datadir.endswith('/'):
            datadir += '/'

        if param is None:
            param = read_param(datadir=datadir, quiet=True)

        f = open(datadir+'index.pro')
        for line in f.readlines():
            clean = line.strip()

            if (not clean.endswith('0') and not clean.startswith('i_') \
               and clean.startswith('i')):
                val = clean.lstrip('i').split('=')
                name = val[0]
                if (name == 'lnTT' and param.ltemperature_nolog):
                    name = 'tt'
                self[name] = int(val[1])
