#$Id: index.py,v 1.1 2008-03-12 15:57:59 tgastine Exp $

class read_index:
   """
     read index.pro and returns a read_index class which is composed
     of a python dictionnary

     Author: T. Gastine (tgastine@ast.obs-mip.fr)
   """
   def __init__(self, datadir='data/'):
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

      f = open(datadir+'index.pro')
      self.index={}
      for line in f.readlines():
        clean = line.strip()

        if (not clean.endswith('0') and not clean.startswith('i_') and clean.startswith('i')):
           val=clean.split('=')
           self.index[val[0].lstrip('i')]=int(val[1])

