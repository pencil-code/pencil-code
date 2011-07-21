#$Id$

from numpy import zeros, asarray
from pencil import read_dim

class read_zprof:
    """ 
    13-mar-2008/dintrans: coded
    f = read_zprof(fname='zprof_hcond.dat',datadir='data/',dim=None,nfield=1)
    Read vertical profiles written in data/proc*/fname
    """

    def __init__(self,fname='zprof_hcond.dat',datadir='data/',dim=None,nfield=1):
        """Constructor:
            -----------
            Params:
            ------
            fname   = we want to read data/proc*/fname files
            datadir = 'data/' (optionnal)
            dim     = None (optional)
            nfield  = number of fields to be read (optional)
            Returns:
            -------
            a class with z and profiles(z)
        """

        if (dim==None): dim=read_dim()
        nz=dim.nzgrid/dim.nprocz
        self.z=zeros(dim.nzgrid,'f')
        self.prof=zeros((nfield,dim.nzgrid),'f')
#
#  loop over all processors and records in file
#
        izcount=0
        for iprocz in range(0,dim.nprocz):
            procname='/proc'+str(iprocz)
            filename=datadir+procname+'/'+fname
            file=open(filename,'r')
            for i in range(nz):
                st=file.readline()
                data=asarray(st.split()).astype('f')
                self.z[izcount]=data[0]
                for j in range(nfield):
                    self.prof[j,izcount]=data[j+1]
                izcount=izcount+1
        file.close()

