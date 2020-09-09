# derived_h5.py
#
# 05-may-20
# Author: F. Gent (fred.gent.ncl@gmail.com).
#
""" Derive auxilliary data and other diagnostics from var.h5 file and 
    save to new h5 file
 
    uses:
      compute 'data' arrays of size [nz,ny,nx] as required
      store 'time' of snapshot
      compute 'masks' for example by temperature phase
      compute summary statistics 'stats' 
      compute 'structure' functions as required
"""
import numpy as np
from scipy.interpolate import interp1d
from ..math import dot, dot2, natural_sort, helmholtz_fft, cpu_optimal
from ..math.derivatives import curl, div, curl2, grad
from ..calc import fluid_reynolds, magnetic_reynolds
from ..io import open_h5, group_h5, dataset_h5
from fileinput import input
from sys import stdout
import subprocess as sub
from .. import read 

def is_vector(key):
    """Check if the variable denoted by the label key is a vector.
    """
    vec = False
    for tag in ('bb', 'jj', 'upot', 'urot', 'vort'):
        if key in tag:
            vec = True
    return vec

def derive_data(sim_path, src, dst, magic=['pp','tt'], par=[],
                gd=[], overwrite=False, rank=0, size=1, nghost=3,
                chunksize = 1000.0, dtype=np.float64, quiet=True  
               ):

    if isinstance(par, list):
        os.chdir(sim_path)
        par = read.param(quiet=True,conflicts_quiet=True)
    if isinstance(gd, list):
        os.chdir(sim_path)
        gd = read.grid(quiet=True)
    #get data dimensions
    nx, ny, nz = src['settings']['nx'][0],\
                 src['settings']['ny'][0],\
                 src['settings']['nz'][0]
    mx, my, mz = src['settings']['mx'][0],\
                 src['settings']['my'][0],\
                 src['settings']['mz'][0]
    #split data into manageable memory chunks
    dstchunksize = 8*nx*ny*nz/1024*1024
    if dstchunksize > chunksize:
        nchunks = cpu_optimal(nx,ny,nz,mvar=1,maux=0,MBmin=chunksize,nmin=8)[1]
        print('nchunks {}'.format(nchunks)) 
    # for mpi split chunks across processes
    if size > 1:
        allchunks = np.arange(nchunks[0]*nchunks[1]*nchunks[2])
        locchunks = np.array_split(allchunks,size)
        ichunks = locchunks[rank]
    else:
        ichunks = np.arange(nchunks[0]*nchunks[1]*nchunks[2])
    indx = np.array_split(np.arange(nx)+nghost,nchunks[0]) 
    indy = np.array_split(np.arange(ny)+nghost,nchunks[1]) 
    indz = np.array_split(np.arange(nz)+nghost,nchunks[2])
    # save time
    if not dst.__contains__('time'):
        dst.create_dataset('time', data=src['time'][()])
    elif overwrite:
        dst.__delitem__('time')
        dst.create_dataset('time', data=src['time'][()])
    # ensure derived variables are in a list
    if isinstance(magic, list):
        magic = magic
    else:
        magic = [magic]
    # initialise group 
    group = group_h5(dst, 'data', status='a', overwrite=overwrite)
    for key in magic:
        if is_vector(key):
            dset = dataset_h5(group, key, status=dst.mode, shape=[3,mz,my,mx],
                          overwrite=True, dtype=dtype)
            print('writing '+key+' shape {}'.format([3,mz,my,mx]))
        else:
            dset = dataset_h5(group, key, status=dst.mode, shape=[mz,my,mx],
                          overwrite=True, dtype=dtype)
            print('writing '+key+' shape {}'.format([mz,my,mx]))
        for iz in np.mod(ichunks,nchunks[2]):
            n1, n2 = indz[iz][ 0]-nghost,\
                     indz[iz][-1]+nghost+1
            n1out = n1+nghost
            n2out = n2-nghost
            varn1 =  nghost
            varn2 = -nghost
            if iz == 0:
                n1out = 0
                varn1 = 0
            if iz == nchunks[2]-1:
                n2out = n2
                varn2 = n2
            for iy in np.mod(ichunks,nchunks[1]):
                m1, m2 = indy[iy][ 0]-nghost,\
                         indy[iy][-1]+nghost+1
                m1out = m1+nghost
                m2out = m2-nghost
                varm1 =  nghost
                varm2 = -nghost
                if iy == 0:
                    m1out = 0
                    varm1 = 0
                if iy == nchunks[1]-1:
                    m2out = m2
                    varm2 = m2
                for ix in np.mod(ichunks,nchunks[0]):
                    l1, l2 = indx[ix][ 0]-nghost,\
                             indx[ix][-1]+nghost+1
                    l1out = l1+nghost
                    l2out = l2-nghost
                    varl1 =  nghost
                    varl2 = -nghost
                    if ix == 0:
                        l1out = 0
                        varl1 = 0
                    if ix == nchunks[0]-1:
                        l2out = l2
                        varl2 = l2
                    if not quiet:
                        print('remeshing '+key+' chunk {}'.format(
                               [iz,iy,ix]))
                    var = calc_derived_data(src['data'], dst['data'], key, par, gd, 
                                   l1, l2, m1, m2, n1, n2)
                    print('var shape {}'.format(var.shape))
                    if not quiet:
                        print('writing '+key+
                                       ' shape {} chunk {}'.format(
                                             var.shape, [iz,iy,ix]))
                    dset[n1out:n2out,
                         m1out:m2out,
                         l1out:l2out] = dtype(var[
                                                varn1:varn2,
                                                varm1:varm2,
                                                varl1:varl2])
#==============================================================================
def calc_derived_data(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2):
    """ 
    compute from src data and existing dst data derived data
    """
    #==========================================================================


    def pressure(src, par, gd, l1, l2, m1, m2, n1, n2):
        if 'rho' in src.keys():
            rho = src['rho'][n1:n2,m1:m2,l1:l2]
        elif 'lnrho' in src.keys():
            rho = np.exp(src['lnrho'][n1:n2,m1:m2,l1:l2])
        else: 
            print('no density used setting rho=1 in pressure calculation')
            rho = 1

        if 'ss' in src.keys():
            ss = src['ss'][n1:n2,m1:m2,l1:l2]
            var = np.exp(par.gamma*(ss + np.log(rho)))
        elif 'tt' in dst.keys():
            tt = dst['tt'][n1:n2,m1:m2,l1:l2]
            if not par.gamma == 1:
                cv = par.cp/par.gamma
            else:
                cv = 0
            var = (par.cp - cv)*tt*rho
        else:
            if 'rho' is src.keys() or 'lnrho' in src.keys():
                print('no entropy or temperature using cs^2'+
                  ' in pressure calculation')
                var = rho*par.cs0**2
            else:
                print('no density or temperature,'+
                      ' pressure cannot be calculated')
                return 1
        return var
     
    #==========================================================================
    def temperature(src, par, gd, l1, l2, m1, m2, n1, n2):
        if 'rho' in src.keys():
            rho = src['rho'][n1:n2,m1:m2,l1:l2]
        elif 'lnrho' in src.keys():
            rho = np.exp(src['lnrho'][n1:n2,m1:m2,l1:l2])
        else: 
            print('no density used setting rho=1 in temperature calculation')
            rho = 1
        if 'ss' in src.keys():
            ss = src['ss'][n1:n2,m1:m2,l1:l2]
            lnrho0 = np.log(par.rho0)
            if not par.gamma == 1:
                lnTT0 = np.log(par.cs0**2/(par.cp*(par.gamma-1)))
                lnTT = lnTT0 + par.gamma/par.cp*ss +\
                           (par.gamma-1)*(np.log(rho)-lnrho0)
            else:
                lnTT0 = np.log(par.cs0**2/par.cp)
                lnTT = lnTT0 + par.gamma/par.cp*ss
        else:
            lnTT0 = np.log(par.cs0**2/(par.cp*(par.gamma-1)))
            lnTT = (par.gamma-1)*(np.log(rho)-lnrho0)+lnTT0
        return np.exp(lnTT)

    #==========================================================================
    def Re_number(src, par, gd, l1, l2, m1, m2, n1, n2):
        uu = np.array([
                      src['ux'][n1:n2,m1:m2,l1:l2],
                      src['uy'][n1:n2,m1:m2,l1:l2],
                      src['uz'][n1:n2,m1:m2,l1:l2]
                      ])
        if 'rho' in src.keys():
            lnrho = np.log(src['rho'][n1:n2,m1:m2,l1:l2])
        elif 'lnrho' in src.keys():
            lnrho = src['lnrho'][n1:n2,m1:m2,l1:l2]
        else: 
            lnrho = list()
        if 'shock' in src.keys():
            shock = src['shock'][n1:n2,m1:m2,l1:l2]
        else:
            shock = list()
        var = fluid_reynolds(uu, par, gd, lnrho=lnrho, shock=shock)
        return var
        
    #==========================================================================
    def Rm_number(src, par, gd, l1, l2, m1, m2, n1, n2):
        uu = np.array([
                      src['ux'][n1:n2,m1:m2,l1:l2],
                      src['uy'][n1:n2,m1:m2,l1:l2],
                      src['uz'][n1:n2,m1:m2,l1:l2]
                      ])
        aa = np.array([
                      src['ax'][n1:n2,m1:m2,l1:l2],
                      src['ay'][n1:n2,m1:m2,l1:l2],
                      src['az'][n1:n2,m1:m2,l1:l2]
                      ])
        if 'bb' in dst.keys():
            bb = dst['bb'][:,n1:n2,m1:m2,l1:l2]
        else:
            bb = bfield(src, par, gd, l1, l2, m1, m2, n1, n2)
        if 'jj' in dst.keys():
            jj = dst['jj'][:,n1:n2,m1:m2,l1:l2]
        else:
            jj = current(src, par, gd, l1, l2, m1, m2, n1, n2)
        var = magnetic_reynolds(uu, par, gd, aa=aa, bb=bb, jj=jj)
        return var

    #==========================================================================
    def Pm_number(src, par, gd, l1, l2, m1, m2, n1, n2):
        if 'Re' in dst.keys():
            Re = dst['Re'][n1:n2,m1:m2,l1:l2]
        else:
            Re = Re_number(src, par, gd, l1, l2, m1, m2, n1, n2)
        if 'Rm' in dst.keys():
            Rm = dst['Rm'][n1:n2,m1:m2,l1:l2]
        else:
            Rm = Rm_number(src, par, gd, l1, l2, m1, m2, n1, n2)
        return Rm/Re


    #==========================================================================
    def rot_flow(src, par, gd, l1, l2, m1, m2, n1, n2):
        uu = np.array([
                      src['ux'][n1:n2,m1:m2,l1:l2],
                      src['uy'][n1:n2,m1:m2,l1:l2],
                      src['uz'][n1:n2,m1:m2,l1:l2]
                      ])
        var = helmholtz_fft(uu, gd, par, rot=True, pot=False)
        print('helmholtz shape {}'.format(var.shape))
        return var

    #==========================================================================
    def pot_flow(src, par, gd, l1, l2, m1, m2, n1, n2):
        uu = np.array([
                      src['ux'][n1:n2,m1:m2,l1:l2],
                      src['uy'][n1:n2,m1:m2,l1:l2],
                      src['uz'][n1:n2,m1:m2,l1:l2]
                      ])
        var = helmholtz_fft(uu, gd, par, pot=True, rot=False)
        return var

    #==========================================================================
    def Mach_cs(src, par, gd, l1, l2, m1, m2, n1, n2):
        uu = np.array([
                      src['ux'][n1:n2,m1:m2,l1:l2],
                      src['uy'][n1:n2,m1:m2,l1:l2],
                      src['uz'][n1:n2,m1:m2,l1:l2]
                      ])
        if 'tt' in dst.keys():
            tt = dst['tt'][n1:n2,m1:m2,l1:l2]
        else:
            tt = temperature(src, par, gd, l1, l2, m1, m2, n1, n2)
        if not par.gamma == 1:
            cs2 = par.cp*(par.gamma-1)*tt
        else:
            cs2 = par.cp*tt
        var = np.sqrt(dot2(uu)/cs2)
        return var

    #==========================================================================
    def Mach_Av(src, par, gd, l1, l2, m1, m2, n1, n2):
        if 'bb' in dst.keys():
            bb = dst['bb'][:,n1:n2,m1:m2,l1:l2]
        else:
            bb = bfield(src, par, gd, l1, l2, m1, m2, n1, n2)
        if 'rho' in src.keys():
            rho = src['rho'][n1:n2,m1:m2,l1:l2]
        elif 'lnrho' in src.keys():
            rho = np.exp(src['lnrho'][n1:n2,m1:m2,l1:l2])
        else: 
            print('no density used setting rho=1 in pressure calculation')
            rho = 1
        var = np.sqrt(dot2(bb)/(par.mu0*rho))
        return var

    #==========================================================================
    def mag_pressure(src, par, gd, l1, l2, m1, m2, n1, n2):
        if 'bb' in dst.keys():
            bb = dst['bb'][:,n1:n2,m1:m2,l1:l2]
        else:
            bb = bfield(src, par, gd, l1, l2, m1, m2, n1, n2)
        var = 0.5*dot2(bb)/par.mu0
        return var

    #==========================================================================
    def vorticity(src, par, gd, l1, l2, m1, m2, n1, n2):
        uu = np.array([
                      src['ux'][n1:n2,m1:m2,l1:l2],
                      src['uy'][n1:n2,m1:m2,l1:l2],
                      src['uz'][n1:n2,m1:m2,l1:l2]
                      ])
        var = curl(uu, gd.dx, gd.dy, gd.dz)
        return var

    #==========================================================================
    def bfield(src, par, gd, l1, l2, m1, m2, n1, n2):
        aa = np.array([
                      src['ax'][n1:n2,m1:m2,l1:l2],
                      src['ay'][n1:n2,m1:m2,l1:l2],
                      src['az'][n1:n2,m1:m2,l1:l2]
                      ])
        var = curl(aa, gd.dx, gd.dy, gd.dz)
        return var

    #==========================================================================
    def current(src, par, gd, l1, l2, m1, m2, n1, n2):
        aa = np.array([
                      src['ax'][n1:n2,m1:m2,l1:l2],
                      src['ay'][n1:n2,m1:m2,l1:l2],
                      src['az'][n1:n2,m1:m2,l1:l2]
                      ])
        var = curl2(aa, gd.dx, gd.dy, gd.dz)
        return var

    #==========================================================================
    def kin_helicity(src, par, gd, l1, l2, m1, m2, n1, n2):
        uu = np.array([
                      src['ux'][n1:n2,m1:m2,l1:l2],
                      src['uy'][n1:n2,m1:m2,l1:l2],
                      src['uz'][n1:n2,m1:m2,l1:l2]
                      ])
        if 'vort' in dst.keys():
            oo = dst['vort'][:,n1:n2,m1:m2,l1:l2]
        else:
            oo = vorticity(src, par, gd, l1, l2, m1, m2, n1, n2)
        var = dot(uu, oo)
        return var

    #==========================================================================
    def mag_helicity(src, par, gd, l1, l2, m1, m2, n1, n2):
        aa = np.array([
                      src['ax'][n1:n2,m1:m2,l1:l2],
                      src['ay'][n1:n2,m1:m2,l1:l2],
                      src['az'][n1:n2,m1:m2,l1:l2]
                      ])
        if 'bb' in dst.keys():
            bb = dst['bb'][:,n1:n2,m1:m2,l1:l2]
        else:
            bb = bfield(src, par, gd, l1, l2, m1, m2, n1, n2)
        var = dot(aa, bb)
        return var

    #==========================================================================
    def calc_derived_item(key):
        case = {
                'urot': rot_flow(src, par, gd, l1, l2, m1, m2, n1, n2),
                'upot': pot_flow(src, par, gd, l1, l2, m1, m2, n1, n2),
                'vort': vorticity(src, par, gd, l1, l2, m1, m2, n1, n2),
                'ou'  : kin_helicity(src, par, gd, l1, l2, m1, m2, n1, n2),
                'tt'  : temperature(src, par, gd, l1, l2, m1, m2, n1, n2),
                'pp'  : pressure(src, par, gd, l1, l2, m1, m2, n1, n2),
                'Re'  : Re_number(src, par, gd, l1, l2, m1, m2, n1, n2),
                'Ms'  : Mach_cs(src, par, gd, l1, l2, m1, m2, n1, n2),
                'bb'  : bfield(src, par, gd, l1, l2, m1, m2, n1, n2),
                'pb'  : mag_pressure(src, par, gd, l1, l2, m1, m2, n1, n2),
                'ab'  : mag_helicity(src, par, gd, l1, l2, m1, m2, n1, n2),
                'jj'  : current(src, par, gd, l1, l2, m1, m2, n1, n2),
                'Ma'  : Mach_Av(src, par, gd, l1, l2, m1, m2, n1, n2),
                'Rm'  : Rm_number(src, par, gd, l1, l2, m1, m2, n1, n2),
                'Pm'  : Pm_number(src, par, gd, l1, l2, m1, m2, n1, n2),
               }
        func = case.get(key, lambda: 'No function for '+key) 
        return func
    #==========================================================================
    return calc_derived_item(key)


#    print('end at {} after {} seconds'.format(
#                                     time.ctime(end_time),end_time-start_time))
# remains to copy other files and edit param files
