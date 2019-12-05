#$Id$

from numpy import zeros, asarray
from .. import read_dim


def read_zprof(*args, **kwargs):
    """Read vertical profiles written in data/proc*/zprof_varname.dat

    Params:
    ------
    varname = we want to read data/proc*/zprof_varname.dat files
    datadir = 'data/' (optionnal)
    dim     = None (optional)
    nfield  = number of fields to be read (optional)

    Returns:
    -------
    A ZProfile object with z and profiles(z)
    """
    return ZProfile(*args, **kwargs)


class ZProfile(object):
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
        nz=int(dim.nzgrid/dim.nprocz)
        self.z=zeros(nz*dim.nprocz,'f')
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
#  when reading a zprof_once_X file, the first dim.nghostz gridpoints are
#  not saved
            if (varname.find('once')!=-1):
                for i in range(dim.nghostz): st=file.readline()
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
