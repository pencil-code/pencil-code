#$Id$

from numpy import zeros, asarray
from pencil import read_dim

class read_zprof:
    """ 
    13-mar-2008/dintrans: coded
    f = read_zprof(varname,datadir='data/',dim=None,nfield=1)
    Read vertical profiles written in data/proc*/zprof_varname.dat
    """

    def __init__(self,varname,datadir='data/',dim=None,nfield=1):
        """Constructor:
            -----------
            Params:
            ------
            varname = we want to read data/proc*/zprof_varname.dat files
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
        if nfield>1:
            self.prof=zeros((nfield,dim.nzgrid),'f')
        else:
            self.prof=zeros(dim.nzgrid,'f')
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
                data=asarray(st.split()).astype('f')
                self.z[izcount]=data[0]
                if nfield>1:
                    for j in range(nfield):
                        self.prof[j,izcount]=data[j+1]
                else:
                    self.prof[izcount]=data[1]
                izcount=izcount+1
        file.close()

