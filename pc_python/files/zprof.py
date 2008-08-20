#$Id: zprof.py,v 1.2 2008-03-13 14:16:27 dintrans Exp $

import numpy as N
import pc_python as pc

class read_zprof:
  """ 
    13-mar-2008/dintrans: coded
    f = read_zprof(varname,datadir='data',dim=dim)
    Read a vertical profile written in data/proc*/zprof_varname.dat
  """

  def __init__(self,varname,datadir='data/',dim=None):
    """Constructor:
         -----------

         Params:
         ------
            varname = we want to read data/proc*/zprof_varname.dat files
            datadir ='data/' (optionnal)
            dim     = None (optional)
         Returns:
         -------
            a read_zprof class with z and prof(z)
    """

    if (dim==None): dim=pc.read_dim()
    nz=dim.nzgrid/dim.nprocz
    self.z=N.zeros(dim.nzgrid,'f')
    self.prof=N.zeros(dim.nzgrid,'f')
#
#  loop over all processors and records in file
#
    izcount=0
    for iprocz in range(0,dim.nprocz):
      procname='/proc'+str(iprocz)
      filename=datadir+procname+'/zprof_'+varname+'.dat'
      file=open(filename,'r')
      for i in range(nz):
        st=file.readline()
        data=N.asarray(st.split()).astype('f')
        self.z[izcount]=data[0]
        self.prof[izcount]=data[1]
        izcount=izcount+1
      file.close()

