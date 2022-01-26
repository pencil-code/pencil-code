# rhs_terms.py
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
from pencil.math import dot, dot2, cross, cpu_optimal
from pencil.math.derivatives import curl, div, curl2, grad, del2
from scipy.ndimage import gaussian_filter as gf
try:
     from pencil.calc import grav_profile
except:
     print('grav_profile not defined: setting to zero')
     def grav_profile(grav, x, y, z):
         if 'z' in grav:
             return np.zeros_like(z)
         if 'y' in grav:
             return np.zeros_like(y)
         if 'x' in grav:
             return np.zeros_like(x)

from pencil.ism_dyn import calc_derived_data
from pencil.io import group_h5, dataset_h5
from pencil import read
import os

def is_vector(key):
    """Check if the variable denoted by the label key is a vector.
    """
    vec = False
    for tag in ('aa', 'uu', 'bb', 'jj', 'upot', 'urot', 'vort', 'uxb', 
                'etadel2a', 'curluxb', 'curletadel2a', 'advec_force',
                'fvisc', 'grav', 'gradp', 'shear', 'coriolis', 'lorentz' ):
        if key == tag:
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

def rhs_data(sim_path, src, dst, magic=['uxb','etadel2a'], par=[], comm=None,
                gd=[], grp_overwrite=False, overwrite=False, 
                rank=0, size=1, nghost=3,status='a',
                chunksize = 1000.0, dtype=np.float64, quiet=True, nmin=32,
                Reynolds_shock=False, lmix=False
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
    ## save time
    #dataset_h5(dst, 'time', status=status, data=src['time'][()],
    #                      comm=comm, size=size, rank=rank,
    #                      overwrite=overwrite, dtype=dtype)
    # ensure derived variables are in a list
    if isinstance(magic, list):
        magic = magic
    else:
        magic = [magic]
    # confirm exists group
    group_h5(dst, 'data', status='a', overwrite=grp_overwrite,
                     comm=comm, rank=rank, size=size)
    # initialise group
    group = group_h5(dst, 'calc', status='a', overwrite=grp_overwrite,
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
                        var = calc_rhs_data(src, dst,
                              key, par, gd, l1, l2, m1, m2, n1, n2,
                              nghost=nghost, Reynolds_shock=Reynolds_shock,
                              lmix=lmix)
                        #print(key,'var shape {}'.format(var.shape))
                        #if not quiet:
                        #    print('writing '+key+
                        #                   ' shape {} chunk {}'.format(
                        #                         var.shape, [iz,iy,ix]))
                        if is_vector(key):
                            dst['calc'][key][:,n1out:n2out,
                                                m1out:m2out,
                                                l1out:l2out] = dtype(var[:,
                                                         varn1:varn2,
                                                         varm1:varm2,
                                                         varl1:varl2])
                        else:
                            dst['calc'][key][n1out:n2out,
                                             m1out:m2out,
                                             l1out:l2out] = dtype(var[
                                                         varn1:varn2,
                                                         varm1:varm2,
                                                         varl1:varl2])
#==============================================================================
def calc_rhs_data(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2,
                      nghost=3, Reynolds_shock=False, lmix=False):
    """ 
    compute from src data and existing dst data derived data
    """
    #==========================================================================
    def induction(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if key=='uxb':
            if 'bb' in dst['data'].keys():
                bb = dst['data/bb'][:,n1:n2,m1:m2,l1:l2]
            else:
                calc_derived_data(src['data'], dst['data'], 'bb', par,
                                    gd, l1, l2, m1, m2, n1, n2, nghost)
                bb = dst['data/bb'][:,n1:n2,m1:m2,l1:l2]
            if 'ux' in src['data'].keys():
                uu = np.array([
                      src['data/ux'][n1:n2,m1:m2,l1:l2],
                      src['data/uy'][n1:n2,m1:m2,l1:l2],
                      src['data/uz'][n1:n2,m1:m2,l1:l2]
                      ])
                var = cross(uu,bb)
            else:
                print('no velocity used setting uxb=0 in induction calculation')
                var = np.zeros_like(bb)
            return var
     
    #==========================================================================
    def resistive(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if key == 'etadel2a':
            n1shift,n2shift,m1shift,m2shift,l1shift,l2shift=der_limits(
                                                     n1,n2,m1,m2,l1,l2,nghost)
            if 'ax' in src['data'].keys():
                aa = np.array([
                      src['data/ax'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift],
                      src['data/ay'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift],
                      src['data/az'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift]
                      ])
                if par.eta == 0.:
                    print('no resistivity used setting etadel2a=0 in resistive')
                    var = np.zeros_like(aa)
                else:
                    var = par.eta*del2(aa, gd.dx, gd.dy, gd.dz,x=gd.x[l1shift:l2shift],y=gd.y[m1shift:m2shift],coordinate_system=par.coord_system)
            else:
                print('no vector potential for resistive calculation')
                for key in src['data'].keys():
                    if len(src['data'][key].shape) == 3:
                        var = np.zeros([3,src['data'][key].shape[0],
                                          src['data'][key].shape[1],
                                          src['data'][key].shape[2] ])
                        break
            n1r,m1r,l1r = under_limits(n1,m1,l1,n1shift,m1shift,l1shift,nghost)
            return var[:,n1r:n2-n1+n1r,m1r:m2-m1+m1r,l1r:l2-l1+l1r]

    #======================================================================
    def curluxb(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if key == 'curluxb':
            n1shift,n2shift,m1shift,m2shift,l1shift,l2shift=der_limits(
                                                     n1,n2,m1,m2,l1,l2,nghost)
            if 'uxb' in dst['calc'].keys():
                tmp = dst['calc/uxb'][:,n1shift:n2shift,
                                        m1shift:m2shift,
                                        l1shift:l2shift]
            else:
                tmp = induction(src, dst, 'uxb', par, gd, n1shift, n2shift,
                                m1shift, m2shift, l1shift, l2shift, nghost)
            var = curl(tmp, gd.dx, gd.dy, gd.dz,x=gd.x[l1shift:l2shift],y=gd.y[m1shift:m2shift],coordinate_system=par.coord_system)
            n1r,m1r,l1r = under_limits(n1,m1,l1,n1shift,m1shift,l1shift,nghost)
            return var[:,n1r:n2-n1+n1r,m1r:m2-m1+m1r,l1r:l2-l1+l1r]
        
    #======================================================================
    def curletadel2a(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if key == 'curletadel2a':
            n1shift,n2shift,m1shift,m2shift,l1shift,l2shift=der_limits(
                                                     n1,n2,m1,m2,l1,l2,nghost)
            if 'etadel2a' in dst['calc'].keys():
                tmp = dst['calc/etadel2a'][:,n1shift:n2shift,
                                        m1shift:m2shift,
                                        l1shift:l2shift]
            else:
                tmp = induction(src, dst, 'etadel2a', par, gd, n1shift, n2shift,
                                m1shift, m2shift, l1shift, l2shift, nghost)
            var = curl(tmp, gd.dx, gd.dy, gd.dz,x=gd.x[l1shift:l2shift],y=gd.y[m1shift:m2shift],coordinate_system=par.coord_system)
            n1r,m1r,l1r = under_limits(n1,m1,l1,n1shift,m1shift,l1shift,nghost)
            return var[:,n1r:n2-n1+n1r,m1r:m2-m1+m1r,l1r:l2-l1+l1r]
        
    #======================================================================
    def bcurluxb(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if key == 'bcurluxb':
            if 'curluxb' in dst['calc'].keys():
                tmp = dst['calc/curluxb'][:,n1:n2, m1:m2, l1:l2]
            else:
                tmp = curluxb(src, dst, 'curluxb', par, gd, n1, n2,
                                m1, m2, l1, l2, nghost)
            if 'bb' in dst['data'].keys():
                bb = dst['data/bb'][:,n1:n2,m1:m2,l1:l2]
            else:
                calc_derived_data(src['data'], dst['data'], 'bb', par,
                                    gd, l1, l2, m1, m2, n1, n2, nghost)
                bb = dst['data/bb'][:,n1:n2,m1:m2,l1:l2]
            var = dot(bb,tmp)
            return var
        
    #======================================================================
    def bcurletadel2a(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if key == 'bcurletadel2a':
            if 'curletadel2a' in dst['calc'].keys():
                tmp = dst['calc/curletadel2a'][:,n1:n2, m1:m2, l1:l2]
            else:
                tmp = curletadel2a(src, dst, 'curletadel2a', par, gd, n1, n2,
                                m1, m2, l1, l2, nghost)
            if 'bb' in dst['data'].keys():
                bb = dst['data/bb'][:,n1:n2,m1:m2,l1:l2]
            else:
                calc_derived_data(src['data'], dst['data'], 'bb', par,
                                    gd, l1, l2, m1, m2, n1, n2, nghost)
                bb = dst['data/bb'][:,n1:n2,m1:m2,l1:l2]
            var = dot(bb,tmp)
            return var

    #==========================================================================
    def urand(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if key == 'u2rand':
            sigma = 0.02
            if 'sigma' in par.keys:
                sigma = par.sigma
            n1shift,n2shift,m1shift,m2shift,l1shift,l2shift=der_limits(
                                                     n1,n2,m1,m2,l1,l2,nghost)
            uu = np.array([
                              src['data/ux'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift],
                              src['data/uy'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift],
                              src['data/uz'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift]
                              ])
            umean = uu.copy()
            urand = uu.copy()
            dx = min(gd.dx,gd.dy,gd.dz)
            sigma_dx = sigma/dx
            print('sigma_dx',sigma_dx)
            for j in range(3):
                umean[j] = gf(uu[j], sigma_dx, mode='reflect')
                urand[j] = uu[j] - umean[j]
            var = dot2(urand)
            n1r,m1r,l1r = under_limits(n1,m1,l1,n1shift,m1shift,l1shift,nghost)
            return var[n1r:n2-n1+n1r,m1r:m2-m1+m1r,l1r:l2-l1+l1r]
        
    #======================================================================
    def gradp(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if key == 'gradp':
            n1shift,n2shift,m1shift,m2shift,l1shift,l2shift=der_limits(
                                                     n1,n2,m1,m2,l1,l2,nghost)
            if 'pp' in dst['data'].keys():
                pp = dst['data/pp'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift]
            else:
                calc_derived_data(src['data'], dst['data'], 'pp', par,
                                    gd, l1, l2, m1, m2, n1, n2, nghost)
                pp = dst['data/pp'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift]
            var = -grad(pp, gd.dx, gd.dy, gd.dz,x=gd.x[l1shift:l2shift],y=gd.y[m1shift:m2shift],coordinate_system=par.coord_system)
            if 'rho' in src['data'].keys():
                rho = dst['data/rho'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift]
            elif 'lnrho' in src['data'].keys():
                rho = np.exp(dst['data/lnrho'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift])
            else:
                rho = 1.
            var /= rho
            n1r,m1r,l1r = under_limits(n1,m1,l1,n1shift,m1shift,l1shift,nghost)
            return var[:,n1r:n2-n1+n1r,m1r:m2-m1+m1r,l1r:l2-l1+l1r]

    #======================================================================
    def grav(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if key == 'grav':
            grav = np.zeros_like(np.array([
                          src['data/ux'][n1:n2,m1:m2,l1:l2],
                          src['data/uy'][n1:n2,m1:m2,l1:l2],
                          src['data/uz'][n1:n2,m1:m2,l1:l2]
                          ]))
            for j, x in zip((1,2,3),('z','y','x')):
                if 'grav{}_profile'.format(x) in par.keys:
                    gg = grav_profile(par.__getattribute__('grav{}_profile'.format(x)), gd.x[l1:l2], gd.y[m1:m2], gd.z[n1:n2], par=par)
                    grav[j] = gg
                if 'grav{}'.format(x) in par.keys:
                    if par.__getattribute__('grav{}'.format(x))>0:
                        grav[j] += par.__getattribute__('grav{}'.format(x))
            if 'potself' in src['data'].keys():
                grav += src['data/potself'][:,n1:n2,m1:m2,l1:l2]
            if 'global_gg' in src['data'].keys():
                grav += src['data/global_gg'][:,n1:n2,m1:m2,l1:l2]
            return grav

    #======================================================================
    def shear_flow(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if key == 'shear':
            uu = np.array([
                          np.zeros_like(src['data/ux'][n1:n2,m1:m2,l1:l2]),
                                        src['data/ux'][n1:n2,m1:m2,l1:l2],
                          np.zeros_like(src['data/ux'][n1:n2,m1:m2,l1:l2])
                          ])
            var = -par.sshear*uu
            return var

    #======================================================================
    def coriolis_flow(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if key == 'coriolis':
            uu = np.array([
                          src['data/ux'][n1:n2,m1:m2,l1:l2],
                          src['data/uy'][n1:n2,m1:m2,l1:l2],
                          src['data/uz'][n1:n2,m1:m2,l1:l2]
                          ])
            omega = np.zeros_like(uu)
            omega[2] = -2*par.omega
            var = cross(omega,uu)
            return var

    #======================================================================
    def lorentz_force(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if key == 'lorentz':
            if 'bb' in dst['data'].keys():
                bb = dst['data/bb'][:,n1:n2,m1:m2,l1:l2]
            else:
                calc_derived_data(src['data'], dst['data'], 'bb', par,
                                    gd, l1, l2, m1, m2, n1, n2, nghost)
                bb = dst['data/bb'][:,n1:n2,m1:m2,l1:l2]
            if 'jj' in dst['data'].keys():
                jj = dst['data/jj'][:,n1:n2,m1:m2,l1:l2]
            else:
                calc_derived_data(src['data'], dst['data'], 'jj', par,
                                    gd, l1, l2, m1, m2, n1, n2, nghost)
                jj = dst['data/jj'][:,n1:n2,m1:m2,l1:l2]
            if 'rho' in src['data'].keys():
                rho = dst['data/rho'][n1:n2,m1:m2,l1:l2]
            elif 'lnrho' in src['data'].keys():
                rho = np.exp(dst['data/lnrho'][n1:n2,m1:m2,l1:l2])
            else:
                rho = 1.
            var = cross(jj,bb)/rho
            return var

    #======================================================================
    def rostrain(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if key == 'fvisc':
            n1shift,n2shift,m1shift,m2shift,l1shift,l2shift=der_limits(
                                                     n1,n2,m1,m2,l1,l2,nghost)
            if par.nu > 0:
                uu = np.array([
                              src['data/ux'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift],
                              src['data/uy'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift],
                              src['data/uz'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift]
                              ])
                #viscous forces
                th2 = 2./3
                th1 = 1./3
                fvisc = np.zeros_like(uu)
                del2u = np.zeros_like(uu)
                for j in range(0,3):
                    del2u[j] = del2(uu[j],gd.dx,gd.dy,gd.dz,x=gd.x[l1shift:l2shift],y=gd.y[m1shift:m2shift],coordinate_system=par.coord_system)
                fvisc += param.nu*del2u
                del(del2u)
                divu = div(uu,grid.dx,grid.dy,grid.dz,x=grid.x[l1shift:l2shift],y=grid.y[m1shift:m2shift],coordinate_system=param.coord_system)
                graddivu = grad(divu,grid.dx,grid.dy,grid.dz,x=grid.x[l1shift:l2shift],y=grid.y[m1shift:m2shift],coordinate_system=param.coord_system)
                fvisc += th1*par.nu*graddivu
                del(graddivu)
                if 'rho' in src['data'].keys():
                    lnrho = np.log(src['data/rho'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift])
                    lrho = True
                elif 'lnrho' in src['data'].keys():
                    lnrho = src['data/lnrho'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift]
                    lrho = True
                else:
                    lrho = False
                if lrho:
                    tmp0 = grad(uu[0],grid.dx,grid.dy,grid.dz,x=grid.x[l1shift:l2shift],y=grid.y[m1shift:m2shift],coordinate_system=param.coord_system)
                    tmp1 = grad(uu[1],grid.dx,grid.dy,grid.dz,x=grid.x[l1shift:l2shift],y=grid.y[m1shift:m2shift],coordinate_system=param.coord_system)
                    tmp2 = grad(uu[2],grid.dx,grid.dy,grid.dz,x=grid.x[l1shift:l2shift],y=grid.y[m1shift:m2shift],coordinate_system=param.coord_system)
                    gradlnrho = grad(lnrho,grid.dx,grid.dy,grid.dz,x=grid.x[l1shift:l2shift],y=grid.y[m1shift:m2shift],coordinate_system=param.coord_system)
                    Sglnrho = np.zeros_like(uu)
                    Sglnrho[0] = dot(tmp0,gradlnrho) +\
                                    (tmp0[0]+tmp1[0]+tmp2[0]-th2*divu)*gradlnrho[0]
                    Sglnrho[1] = dot(tmp1,gradlnrho) +\
                                    (tmp0[1]+tmp1[1]+tmp2[1]-th2*divu)*gradlnrho[1]
                    Sglnrho[2] = dot(tmp2,gradlnrho) +\
                                    (tmp0[2]+tmp1[2]+tmp2[2]-th2*divu)*gradlnrho[2]
                    fvisc += par.nu*Sglnrho
                    del(gradlnrho,Sglnrho)
                del(divu)
            n1r,m1r,l1r = under_limits(n1,m1,l1,n1shift,m1shift,l1shift,nghost)
            return fvisc[:,n1r:n2-n1+n1r,m1r:m2-m1+m1r,l1r:l2-l1+l1r]

    #==========================================================================
    def advec_force(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if key == 'uadvec':
            n1shift,n2shift,m1shift,m2shift,l1shift,l2shift=der_limits(
                                                     n1,n2,m1,m2,l1,l2,nghost)
            uu = np.array([
                              src['data/ux'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift],
                              src['data/uy'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift],
                              src['data/uz'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift]
                              ])
            tmp0 = grad(uu[0],grid.dx,grid.dy,grid.dz,x=grid.x[l1shift:l2shift],y=grid.y[m1shift:m2shift],coordinate_system=param.coord_system)
            tmp1 = grad(uu[1],grid.dx,grid.dy,grid.dz,x=grid.x[l1shift:l2shift],y=grid.y[m1shift:m2shift],coordinate_system=param.coord_system)
            tmp2 = grad(uu[2],grid.dx,grid.dy,grid.dz,x=grid.x[l1shift:l2shift],y=grid.y[m1shift:m2shift],coordinate_system=param.coord_system)
            advec = np.zeros_like(uu)
            advec[0] = -dot(uu,tmp0)
            advec[1] = -dot(uu,tmp1)
            advec[2] = -dot(uu,tmp2)
            del(tmp0,tmp1,tmp2)
            n1r,m1r,l1r = under_limits(n1,m1,l1,n1shift,m1shift,l1shift,nghost)
            return advec[:,n1r:n2-n1+n1r,m1r:m2-m1+m1r,l1r:l2-l1+l1r]

    #==========================================================================
    def advec_heat(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if key == 'hadvec':
            n1shift,n2shift,m1shift,m2shift,l1shift,l2shift=der_limits(
                                                     n1,n2,m1,m2,l1,l2,nghost)
            uu = np.array([
                              src['data/ux'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift],
                              src['data/uy'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift],
                              src['data/uz'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift]
                              ])
            if 'ss' in dst['data'].keys():
                ss = src['ss'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift]
            sgrad = grad(ss,grid.dx,grid.dy,grid.dz,x=grid.x[l1shift:l2shift],y=grid.y[m1shift:m2shift],coordinate_system=param.coord_system)
            if 'tt' in dst['data'].keys():
                tt = dst['tt'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift]
            else:
                calc_derived_data(src['data'], dst['data'], 'tt', par,
                                    gd, l1, l2, m1, m2, n1, n2, nghost)
                tt = dst['data/tt'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift]
            if 'rho' in dst['data'].keys():
                rho = src['data/rho'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift]
            elif 'lnrho' in dst['data'].keys():
                rho = np.exp(src['data/lnrho'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift])
            else:
                rho=1
            advec = -dot(uu,sgrad)/(rho*tt)
            n1r,m1r,l1r = under_limits(n1,m1,l1,n1shift,m1shift,l1shift,nghost)
            return advec[n1r:n2-n1+n1r,m1r:m2-m1+m1r,l1r:l2-l1+l1r]

    #==========================================================================
    def ohmic_heat(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if key == 'ohmic':
            if 'jj' in dst['data'].keys():
                jj = dst['data/jj'][:,n1:n2,m1:m2,l1:l2]
            else:
                calc_derived_data(src['data'], dst['data'], 'jj', par,
                                    gd, l1, l2, m1, m2, n1, n2, nghost)
                jj = dst['data/jj'][:,n1:n2,m1:m2,l1:l2]
            if 'tt' in dst['data'].keys():
                tt = dst['tt'][n1:n2,m1:m2,l1:l2]
            else:
                calc_derived_data(src['data'], dst['data'], 'tt', par,
                                    gd, l1, l2, m1, m2, n1, n2, nghost)
                tt = dst['data/tt'][n1:n2,m1:m2,l1:l2]
            if 'rho' in dst['data'].keys():
                rho = src['data/rho'][n1:n2,m1:m2,l1:l2]
            elif 'lnrho' in dst['data'].keys():
                rho = np.exp(src['data/lnrho'][n1:n2,m1:m2,l1:l2])
            else:
                rho=1
            var = par.eta*par.mu0*dot2(jj)/(rho*tt)
            return var

    #==========================================================================
    def visc_heat(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if key == 'hvisc':
            th1 = 1./3
            n1shift,n2shift,m1shift,m2shift,l1shift,l2shift=der_limits(
                                                     n1,n2,m1,m2,l1,l2,nghost)
            uu = np.array([
                          src['data/ux'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift],
                          src['data/uy'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift],
                          src['data/uz'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift]
                          ])
            tmp0 = grad(uu[0],grid.dx,grid.dy,grid.dz,x=grid.x[l1shift:l2shift],y=grid.y[m1shift:m2shift],coordinate_system=param.coord_system)
            tmp1 = grad(uu[1],grid.dx,grid.dy,grid.dz,x=grid.x[l1shift:l2shift],y=grid.y[m1shift:m2shift],coordinate_system=param.coord_system)
            tmp2 = grad(uu[2],grid.dx,grid.dy,grid.dz,x=grid.x[l1shift:l2shift],y=grid.y[m1shift:m2shift],coordinate_system=param.coord_system)
            var  = tmp0[1]**2
            var += tmp0[2]**2
            var += tmp1[2]**2
            if 'tt' in dst['data'].keys():
                tt = dst['tt'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift]
            else:
                calc_derived_data(src['data'], dst['data'], 'tt', par,
                                    gd, l1, l2, m1, m2, n1, n2, nghost)
                tt = dst['data/tt'][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift]
            var *= 2*par.nu*var/tt
            n1r,m1r,l1r = under_limits(n1,m1,l1,n1shift,m1shift,l1shift,nghost)
            return var[n1r:n2-n1+n1r,m1r:m2-m1+m1r,l1r:l2-l1+l1r]

    #==========================================================================
    def heatcond(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if key == 'hcond':
            n1shift,n2shift,m1shift,m2shift,l1shift,l2shift=der_limits(
                                                     n1,n2,m1,m2,l1,l2,nghost)
            if 'tt' in dst['data'].keys():
                tt = dst['data/tt'][n1shift,n2shift,m1shift,m2shift,l1shift,l2shift]
            else:
                calc_derived_data(src['data'], dst['data'], 'tt', par,
                                    gd, l1, l2, m1, m2, n1, n2, nghost)
                tt = dst['data/tt'][n1shift,n2shift,m1shift,m2shift,l1shift,l2shift]
                lntt = np.log(dst['data/tt'][n1shift,n2shift,m1shift,m2shift,l1shift,l2shift])
                gradlnT = grad(lntt,grid.dx,grid.dy,grid.dz,x=grid.x[l1shift:l2shift],y=grid.y[m1shift:m2shift],coordinate_system=param.coord_system)
                del2T = del2(tt,grid.dx,grid.dy,grid.dz,x=grid.x[l1shift:l2shift],y=grid.y[m1shift:m2shift],coordinate_system=param.coord_system)
            if 'rho' in src['data'].keys():
                lnrho = np.log(src['rho'][n1shift,n2shift,m1shift,m2shift,l1shift,l2shift])
                gradlnrho = grad(lnrho,grid.dx,grid.dy,grid.dz,x=grid.x[l1shift:l2shift],y=grid.y[m1shift:m2shift],coordinate_system=param.coord_system)
                lrho = True
            elif 'lnrho' in src['data'].keys():
                lnrho = src['lnrho'][n1shift,n2shift,m1shift,m2shift,l1shift,l2shift]
                gradlnrho = grad(lnrho,grid.dx,grid.dy,grid.dz,x=grid.x[l1shift:l2shift],y=grid.y[m1shift:m2shift],coordinate_system=param.coord_system)
                lrho = True
            else:
                lrho = False
            var = par.cp*par.chi*del2T/tt
            if lrho:
                var += par.cp*par.chi*dot(gradlnrho,gradlnT)
            n1r,m1r,l1r = under_limits(n1,m1,l1,n1shift,m1shift,l1shift,nghost)
            return var[n1r:n2-n1+n1r,m1r:m2-m1+m1r,l1r:l2-l1+l1r]


    #======================================================================
    def calc_rhs_item(key):
        case = {
                'etadel2a'     : resistive(     src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost),
                'uxb'          : induction(     src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost),
                'curluxb'      : curluxb(       src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost),
                'curletadel2a' : curletadel2a(  src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost),
                'bcurluxb'     : bcurluxb(      src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost),
                'bcurletadel2a': bcurletadel2a( src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost),
                'gradp'        : gradp(         src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost),
                'grav'         : grav(          src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost),
                'shear'        : shear_flow(    src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost),
                'coriolis'     : coriolis_flow( src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost),
                'lorentz'      : lorentz_force( src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost),
                'fvisc'        : rostrain(      src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost),
                'uadvec'       : advec_force(   src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost),
                'hadvec'       : advec_heat(    src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost),
                'ohmic'        : ohmic_heat(    src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost),
                'hvisc'        : visc_heat(     src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost),
                'hcond'        : heatcond(      src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost),
                'u2rand'       : urand(         src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost),
               }
        func = case.get(key, lambda: 'No function for '+key)
        return func
        #======================================================================
    return calc_rhs_item(key)


#    print('end at {} after {} seconds'.format(
#                                     time.ctime(end_time),end_time-start_time))
# remains to copy other files and edit param files
