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
import os

def is_vector(key):
    """Check if the variable denoted by the label key is a vector.
    """
    vec = False
    for tag in ('aa', 'uu', 'bb', 'jj', 'upot', 'urot', 'vort'):
        if key in tag:
            vec = True
    return vec

def der_limits(n1,n2,m1,m2,l1,l2,nghost):
    if n1 >= nghost:
        n1shift = n1 - nghost
    else:
        n1shift = n1
    if m1 >= nghost:
        m1shift = m1 - nghost
    else:
        m1shift = m1
    if l1 >= nghost:
        l1shift = l1 - nghost
    else:
        l1shift = l1
    n2shift = n2 + nghost
    m2shift = m2 + nghost
    l2shift = l2 + nghost
    return n1shift,n2shift,m1shift,m2shift,l1shift,l2shift

def under_limits(n1,m1,l1,n1shift,m1shift,l1shift,nghost):
    if n1shift == n1:
        n1r = 0
    else:
        n1r = nghost
    if m1shift == m1:
        m1r = 0
    else:
        m1r = nghost
    if l1shift == l1:
        l1r = 0
    else:
        l1r = nghost
    return n1r,m1r,l1r

def derive_data(sim_path, src, dst, magic=['pp','tt'], par=[], comm=None,
                gd=[], overwrite=False, rank=0, size=1, nghost=3,status='a',
                chunksize = 1000.0, dtype=np.float64, quiet=True, nmin=32
               ):

    if comm:
        overwrite = False
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
        nchunks = cpu_optimal(nx,ny,nz,quiet=quiet,
                              mvar=src['settings/mvar'][0],
                              maux=src['settings/maux'][0],
                              MBmin=chunksize,nmin=nmin,size=size)[1]
    else:
        nchunks = [1,1,1]
    print('nchunks {}'.format(nchunks))
    # for mpi split chunks across processes
    if size > 1:
        locindx = np.array_split(np.arange(nx)+nghost,nchunks[0])
        locindy = np.array_split(np.arange(ny)+nghost,nchunks[1])
        locindz = np.array_split(np.arange(nz)+nghost,nchunks[2])
        indx = [locindx[np.mod(rank+int(rank/nchunks[2])
                                   +int(rank/nchunks[1]),nchunks[0])]]
        indy = [locindy[np.mod(rank+int(rank/nchunks[2]),nchunks[1])]]
        indz = [locindz[np.mod(rank,nchunks[2])]]
        allchunks = 1
    else:
        locindx = np.array_split(np.arange(nx)+nghost,nchunks[0])
        locindy = np.array_split(np.arange(ny)+nghost,nchunks[1])
        locindz = np.array_split(np.arange(nz)+nghost,nchunks[2])
        indx = np.array_split(np.arange(nx)+nghost,nchunks[0])
        indy = np.array_split(np.arange(ny)+nghost,nchunks[1])
        indz = np.array_split(np.arange(nz)+nghost,nchunks[2])
        allchunks = nchunks[0]*nchunks[1]*nchunks[2]
    # save time
    dataset_h5(dst, 'time', status=status, data=src['time'][()],
                          comm=comm, size=size, rank=rank,
                          overwrite=overwrite, dtype=dtype)
    # ensure derived variables are in a list
    if isinstance(magic, list):
        magic = magic
    else:
        magic = [magic]
    # initialise group
    group = group_h5(dst, 'data', status='a', overwrite=overwrite,
                     comm=comm, rank=rank, size=size)
    for key in magic:
        if is_vector(key):
            dataset_h5(group, key, status=status, shape=[3,mz,my,mx],
                          comm=comm, size=size, rank=rank,
                          overwrite=overwrite, dtype=dtype)
            print('writing '+key+' shape {}'.format([3,mz,my,mx]))
        else:
            dataset_h5(group, key, status=status, shape=[mz,my,mx],
                          comm=comm, size=size, rank=rank,
                          overwrite=overwrite, dtype=dtype)
            print('writing '+key+' shape {}'.format([mz,my,mx]))
        for ichunk in range(allchunks):
            for iz in [indz[np.mod(ichunk,nchunks[2])]]:
                n1, n2 = iz[ 0]-nghost,\
                         iz[-1]+nghost+1
                n1out = n1+nghost
                n2out = n2-nghost
                varn1 =  nghost
                varn2 = -nghost
                if iz[0] == locindz[0][0]:
                    n1out = 0
                    varn1 = 0
                if iz[-1] == locindz[-1][-1]:
                    n2out = n2
                    varn2 = n2
                for iy in [indy[np.mod(ichunk+
                                   int(ichunk/nchunks[2]),nchunks[1])]]:
                    m1, m2 = iy[ 0]-nghost,\
                             iy[-1]+nghost+1
                    m1out = m1+nghost
                    m2out = m2-nghost
                    varm1 =  nghost
                    varm2 = -nghost
                    if iy[0] == locindy[0][0]:
                        m1out = 0
                        varm1 = 0
                    if iy[-1] == locindy[-1][-1]:
                        m2out = m2
                        varm2 = m2
                    for ix in [indx[np.mod(ichunk+int(ichunk/nchunks[2])
                                   +int(ichunk/nchunks[1]),nchunks[0])]]:
                        l1, l2 = ix[ 0]-nghost,\
                                 ix[-1]+nghost+1
                        l1out = l1+nghost
                        l2out = l2-nghost
                        varl1 =  nghost
                        varl2 = -nghost
                        if ix[0] == locindx[0][0]:
                            l1out = 0
                            varl1 = 0
                        if ix[-1] == locindx[-1][-1]:
                            l2out = l2
                            varl2 = l2
                        if not quiet:
                            print('remeshing '+key+' chunk {}'.format(
                                   [iz,iy,ix]))
                        var = calc_derived_data(src['data'], dst['data'],
                              key, par, gd, l1, l2, m1, m2, n1, n2, nghost=nghost)
                        #print('var shape {}'.format(var.shape))
                        #if not quiet:
                        #    print('writing '+key+
                        #                   ' shape {} chunk {}'.format(
                        #                         var.shape, [iz,iy,ix]))
                        if is_vector(key):
                             dst['data'][key][:,n1out:n2out,
                                                m1out:m2out,
                                                l1out:l2out] = dtype(var[:,
                                                         varn1:varn2,
                                                         varm1:varm2,
                                                         varl1:varl2])
                        else:
                            dst['data'][key][n1out:n2out,
                                             m1out:m2out,
                                             l1out:l2out] = dtype(var[
                                                         varn1:varn2,
                                                         varm1:varm2,
                                                         varl1:varl2])
#==============================================================================
def calc_derived_data(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2,
                      nghost=3):
    """
    compute from src data and existing dst data derived data
    """
    #==========================================================================
    def pressure(src, dst, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if 'rho' in src.keys():
            rho = src['rho'][n1:n2,m1:m2,l1:l2]
        elif 'lnrho' in src.keys():
            rho = np.exp(src['lnrho'][n1:n2,m1:m2,l1:l2])
        else:
            print('no density used setting rho=1 in pressure calculation')
            rho = 1
        if not par.gamma == 1:
            lnTT0 = np.log(par.cs0**2/(par.cp*(par.gamma-1)))
            cv = par.cp/par.gamma
        else:
            lnTT0 = np.log(par.cs0**2/par.cp)
            cv = 1
        lnrho0 = np.log(par.rho0)
        cv1 = 1./cv
        gamma_m1 = par.gamma-1.

        if 'ss' in src.keys():
            ss = src['ss'][n1:n2,m1:m2,l1:l2]
            var = (par.cp - cv)*np.exp(cv1*ss +
                   par.gamma*np.log(rho)-gamma_m1*lnrho0)
        elif 'tt' in dst.keys():
            tt = dst['tt'][n1:n2,m1:m2,l1:l2]
            var = (par.cp - cv)*tt*rho
        else:
            if 'rho' in src.keys() or 'lnrho' in src.keys():
                print('no entropy or temperature using cs^2'+
                  ' in pressure calculation')
                var = rho*par.cs0**2
            else:
                print('no density or temperature,'+
                      ' pressure cannot be calculated')
                return 1
        return var

    #==========================================================================
    def temperature(src, dst, par, gd, l1, l2, m1, m2, n1, n2, nghost):
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

    #======================================================================
    def Re_number(src, dst, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        n1shift,n2shift,m1shift,m2shift,l1shift,l2shift=der_limits(
                                                 n1,n2,m1,m2,l1,l2,nghost)
        uu = np.array([
                      src['ux'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift],
                      src['uy'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift],
                      src['uz'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift]
                      ])
        if 'rho' in src.keys():
            lnrho = np.log(src['rho'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift])
        elif 'lnrho' in src.keys():
            lnrho = src['lnrho'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift]
        else:
            lnrho = list()
        if 'shock' in src.keys():
            shock = src['shock'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift]
        else:
            shock = list()
        var = fluid_reynolds(uu, par, gd, lnrho=lnrho, shock=shock)
        n1r,m1r,l1r = under_limits(n1,m1,l1,n1shift,m1shift,l1shift,nghost)
        return var[n1r:n2-n1+n1r,m1r:m2-m1+m1r,l1r:l2-l1+l1r]

    #======================================================================
    def Rm_number(src, dst, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        n1shift,n2shift,m1shift,m2shift,l1shift,l2shift=der_limits(
                                                 n1,n2,m1,m2,l1,l2,nghost)
        uu = np.array([
                      src['ux'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift],
                      src['uy'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift],
                      src['uz'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift]
                      ])
        aa = np.array([
                      src['ax'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift],
                      src['ay'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift],
                      src['az'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift]
                      ])
        if 'bb' in dst.keys():
            bb = dst['bb'][:,n1shift:n2shift,m1shift:m2shift,l1shift:l2shift]
        else:
            bb = bfield(src, dst, par, gd, l1shift, l2shift, m1shift, m2shift, n1shift, n2shift, nghost)
        if 'jj' in dst.keys():
            jj = dst['jj'][:,n1shift:n2shift,m1shift:m2shift,l1shift:l2shift]
        else:
            jj = current(src, dst, par, gd, l1shift, l2shift, m1shift, m2shift, n1shift, n2shift, nghost)
        var = magnetic_reynolds(uu, par, gd, aa=aa, bb=bb, jj=jj)
        n1r,m1r,l1r = under_limits(n1,m1,l1,n1shift,m1shift,l1shift,nghost)
        return var[n1r:n2-n1+n1r,m1r:m2-m1+m1r,l1r:l2-l1+l1r]

    #======================================================================
    def Pm_number(src, dst, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if 'Re' in dst.keys():
            Re = dst['Re'][n1:n2,m1:m2,l1:l2]
        else:
            Re = Re_number(src, dst, par, gd, l1, l2, m1, m2, n1, n2, nghost)
        if 'Rm' in dst.keys():
            Rm = dst['Rm'][n1:n2,m1:m2,l1:l2]
        else:
            Rm = Rm_number(src, dst, par, gd, l1, l2, m1, m2, n1, n2, nghost)
        if Re.max() > 0:
            Re[np.where(Re==0)] = Re[np.where(Re>0)].min()
        else:
            if Rm.max() > 0:
                Re[:] = Rm.min()
            else:
                Re[:] = 1
        return Rm/Re

    #======================================================================
    def rot_flow(src, dst, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        uu = np.array([
                      src['ux'][n1:n2,m1:m2,l1:l2],
                      src['uy'][n1:n2,m1:m2,l1:l2],
                      src['uz'][n1:n2,m1:m2,l1:l2]
                      ])
        var = helmholtz_fft(uu, gd, par, rot=True, pot=False)
        #print('helmholtz shape {}'.format(var.shape))
        return var

    #======================================================================
    def pot_flow(src, dst, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        uu = np.array([
                      src['ux'][n1:n2,m1:m2,l1:l2],
                      src['uy'][n1:n2,m1:m2,l1:l2],
                      src['uz'][n1:n2,m1:m2,l1:l2]
                      ])
        var = helmholtz_fft(uu, gd, par, pot=True, rot=False)
        return var

    #======================================================================
    def Mach_cs(src, dst, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        uu = np.array([
                      src['ux'][n1:n2,m1:m2,l1:l2],
                      src['uy'][n1:n2,m1:m2,l1:l2],
                      src['uz'][n1:n2,m1:m2,l1:l2]
                      ])
        if 'tt' in dst.keys():
            tt = dst['tt'][n1:n2,m1:m2,l1:l2]
        else:
            tt = temperature(src, dst, par, gd, l1, l2, m1, m2, n1, n2, nghost)
        if not par.gamma == 1:
            cs2 = par.cp*(par.gamma-1)*tt
        else:
            cs2 = par.cp*tt
        print('tt min {} max {}'.format(tt.min(),tt.max()))
        var = np.sqrt(dot2(uu)/cs2)
        return var

    #======================================================================
    def Mach_Av(src, dst, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if 'bb' in dst.keys():
            bb = dst['bb'][:,n1:n2,m1:m2,l1:l2]
        else:
            bb = bfield(src, dst, par, gd, l1, l2, m1, m2, n1, n2, nghost)
        if 'rho' in src.keys():
            rho = src['rho'][n1:n2,m1:m2,l1:l2]
        elif 'lnrho' in src.keys():
            rho = np.exp(src['lnrho'][n1:n2,m1:m2,l1:l2])
        else:
            print('no density used setting rho=1 in pressure calculation')
            rho = 1
        var = np.sqrt(dot2(bb)/(par.mu0*rho))
        return var

    #======================================================================
    def mag_pressure(src, dst, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if 'bb' in dst.keys():
            bb = dst['bb'][:,n1:n2,m1:m2,l1:l2]
        else:
            bb = bfield(src, dst, par, gd, l1, l2, m1, m2, n1, n2, nghost)
        var = 0.5*dot2(bb)/par.mu0
        return var

    #======================================================================
    def vorticity(src, dst, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        n1shift,n2shift,m1shift,m2shift,l1shift,l2shift=der_limits(
                                                 n1,n2,m1,m2,l1,l2,nghost)
        uu = np.array([
                      src['ux'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift],
                      src['uy'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift],
                      src['uz'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift]
                      ])
        var = curl(uu, gd.dx, gd.dy, gd.dz)
        n1r,m1r,l1r = under_limits(n1,m1,l1,n1shift,m1shift,l1shift,nghost)
        return var[:,n1r:n2-n1+n1r,m1r:m2-m1+m1r,l1r:l2-l1+l1r]

    #======================================================================
    def bfield(src, dst, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        n1shift,n2shift,m1shift,m2shift,l1shift,l2shift=der_limits(
                                                 n1,n2,m1,m2,l1,l2,nghost)
        aa = np.array([
                      src['ax'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift],
                      src['ay'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift],
                      src['az'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift]
                      ])
        var = curl(aa, gd.dx, gd.dy, gd.dz)
        n1r,m1r,l1r = under_limits(n1,m1,l1,n1shift,m1shift,l1shift,nghost)
        return var[:,n1r:n2-n1+n1r,m1r:m2-m1+m1r,l1r:l2-l1+l1r]

    #======================================================================
    def current(src, dst, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        n1shift,n2shift,m1shift,m2shift,l1shift,l2shift=der_limits(
                                                 n1,n2,m1,m2,l1,l2,nghost)
        aa = np.array([
                      src['ax'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift],
                      src['ay'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift],
                      src['az'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift]
                      ])
        var = curl2(aa, gd.dx, gd.dy, gd.dz)
        n1r,m1r,l1r = under_limits(n1,m1,l1,n1shift,m1shift,l1shift,nghost)
        return var[:,n1r:n2-n1+n1r,m1r:m2-m1+m1r,l1r:l2-l1+l1r]

    #======================================================================
    def kin_helicity(src, dst, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        uu = np.array([
                      src['ux'][n1:n2,m1:m2,l1:l2],
                      src['uy'][n1:n2,m1:m2,l1:l2],
                      src['uz'][n1:n2,m1:m2,l1:l2]
                      ])
        if 'vort' in dst.keys():
            oo = dst['vort'][:,n1:n2,m1:m2,l1:l2]
        else:
            oo = vorticity(src, dst, par, gd, l1, l2, m1, m2, n1, n2, nghost)
        var = dot(uu, oo)
        return var

    #======================================================================
    def mag_helicity(src, dst, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        aa = np.array([
                      src['ax'][n1:n2,m1:m2,l1:l2],
                      src['ay'][n1:n2,m1:m2,l1:l2],
                      src['az'][n1:n2,m1:m2,l1:l2]
                      ])
        if 'bb' in dst.keys():
            bb = dst['bb'][:,n1:n2,m1:m2,l1:l2]
        else:
            bb = bfield(src, dst, par, gd, l1, l2, m1, m2, n1, n2, nghost)
        var = dot(aa, bb)
        return var

    #==========================================================================
    def calc_derived_item(key):
        case = {
                'urot': rot_flow(    src, dst, par, gd, l1, l2, m1, m2, n1, n2, nghost),
                'upot': pot_flow(    src, dst, par, gd, l1, l2, m1, m2, n1, n2, nghost),
                'vort': vorticity(   src, dst, par, gd, l1, l2, m1, m2, n1, n2, nghost),
                'ou'  : kin_helicity(src, dst, par, gd, l1, l2, m1, m2, n1, n2, nghost),
                'tt'  : temperature( src, dst, par, gd, l1, l2, m1, m2, n1, n2, nghost),
                'pp'  : pressure(    src, dst, par, gd, l1, l2, m1, m2, n1, n2, nghost),
                'Re'  : Re_number(   src, dst, par, gd, l1, l2, m1, m2, n1, n2, nghost),
                'Ms'  : Mach_cs(     src, dst, par, gd, l1, l2, m1, m2, n1, n2, nghost),
                'bb'  : bfield(      src, dst, par, gd, l1, l2, m1, m2, n1, n2, nghost),
                'pb'  : mag_pressure(src, dst, par, gd, l1, l2, m1, m2, n1, n2, nghost),
                'ab'  : mag_helicity(src, dst, par, gd, l1, l2, m1, m2, n1, n2, nghost),
                'jj'  : current(     src, dst, par, gd, l1, l2, m1, m2, n1, n2, nghost),
                'Ma'  : Mach_Av(     src, dst, par, gd, l1, l2, m1, m2, n1, n2, nghost),
                'Rm'  : Rm_number(   src, dst, par, gd, l1, l2, m1, m2, n1, n2, nghost),
                'Pm'  : Pm_number(   src, dst, par, gd, l1, l2, m1, m2, n1, n2, nghost),
               }
        func = case.get(key, lambda: 'No function for '+key)
        return func
        #======================================================================
    return calc_derived_item(key)


#    print('end at {} after {} seconds'.format(
#                                     time.ctime(end_time),end_time-start_time))
# remains to copy other files and edit param files
