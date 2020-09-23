#!/usr/bin/env python
#coding: utf-8
# Helmholtz.py
# Written by Fred Gent (fred.gent.ncl@gmail.com)
"""
Perform a Helmholtz decomposition on a vector field returning a pair of 
orthogonal vectors with zero divergence and zero curl respectively.
"""

import numpy as np
from ..math import dot, dot2, cross
from ..math.derivatives import div, curl, grad

def helmholtz_fft(tot_field, grid, params, nghost=3, pot=True, rot=True,
                  lno_mean=False, nonperi_bc=None, field_scalar=[], s=None,
                  quiet=True):
    """
    Creates the decomposition vector pair for the supplied vector field.

    call signature:

    helmholz_fft(field, grid, params, )
 
    Keyword arguments:

     *tot_field*:
       Vector field of dimension [3,mz,my,mx], which is decomposed.

     *grid*:
       Grid object with grid spacing dx, dy, dz.

     *params*:
       Simulation Params object with domain dimensions Lx, Ly and Lz.

     *nghost*:
       Number of ghost zones to exclude from the fft.

     *lno_mean*:
       Exclude any mean flow from the decomposition - should drop anyway.

     *nonperi_bc:
       String if not None with boundary condition label.
       How to apply the bc needs to be implemented as required.

     *field_scalar
       Scalar field (density) as debug tool for energy comparison.

     *s*:
      List of three integers if not None for fft dimension.
      If none the dimension of the field [nz,ny,nx] is used.

    """
    if lno_mean:
        # exclude volume mean flows
        field = np.zeros_like(tot_field)
        for j in range(0,3):
            field[j] = tot_field[j] - tot_field[j].mean()
    else:
        field = tot_field
    #use mean speed and grid spacing in normalization of div/curl check
    amp_field_1 = 1./np.sqrt(dot2(field)).mean()
    nz,ny,nx = field[:,nghost:-nghost,nghost:-nghost,nghost:-nghost].shape[-3],\
               field[:,nghost:-nghost,nghost:-nghost,nghost:-nghost].shape[-2],\
               field[:,nghost:-nghost,nghost:-nghost,nghost:-nghost].shape[-1]
    invs = [nz,ny,nx]
    if not s:
        s =  [nz,ny,nx]
    knz,kny,knx = s[0], s[1], s[2]
    nz2,ny2,nx2 = int(knz/2),int(kny/2),int(knx/2)
    #derive wavenumbers k scaled to dimension of simulation domain
    kk = np.empty(shape=[3,s[0],s[1],s[2]])
    k0 = np.arange(knx); k0[nx2:] = -k0[nx2-1::-1]-1
    k1 = np.arange(kny); k1[ny2:] = -k1[ny2-1::-1]-1
    k2 = np.arange(knz); k2[nz2:] = -k2[nz2-1::-1]-1
    for j in range(0,k0.size):
        kk[0,:,:,j] = k0[j]*2*np.pi/params.lxyz[0]
    for j in range(0,k1.size):
        kk[1,:,j,:] = k1[j]*2*np.pi/params.lxyz[1]
    for j in range(0,k2.size):
        kk[2,j,:,:] = k2[j]*2*np.pi/params.lxyz[2]
    knorm = dot2(kk)
    #apply fast Fourier transform to the vector field
    kfield = np.empty(shape=[3,s[0],s[1],s[2]], dtype=complex)
    for j in range(0,3):
        kfield[j] = np.fft.fftn(
                    field[j,nghost:-nghost,nghost:-nghost,nghost:-nghost],s=s
                                )
    if pot:
        #reverse fft to obtain the scalar potential
        pfield = -1j*dot(kk,kfield)
        pfield[np.where(knorm==0)] = 0.
    if rot:
        #reverse fft to obtain the vector potential
        rfield =  1j*cross(kk,kfield)
        for j in range(3):
            rfield[j][np.where(knorm==0)] = 0.
    #avoid division by zero
    knorm[np.where(knorm==0)] = 1.
    if pot:
        pfield /= knorm
        pot_field = np.zeros_like(field)
    if rot:
        for j in range(3):
            rfield[j] /= knorm
        rot_field = np.zeros_like(field)
    if nonperi_bc:
        print('Please implement new nonperi_bc not yet implemented.\n',
              'Applying periodic boundary conditions for now.')
    for j in range(0,3):
        if pot:
            pot_field[j,nghost:-nghost,nghost:-nghost,nghost:-nghost] =\
                                      np.fft.ifftn(1j*pfield*kk[j],s=invs).real
        if rot:
            rot_field[j,nghost:-nghost,nghost:-nghost,nghost:-nghost] =\
                               np.fft.ifftn(cross(1j*kk,rfield)[j],s=invs).real
    #apply the periodic boundary conditions for the ghost zones:
    for j in range(0,3):
        if pot:
            pot_field[j, :nghost,:,:] = pot_field[j,-2*nghost:-nghost,:,:]
            pot_field[j,-nghost:,:,:] = pot_field[j, nghost: 2*nghost,:,:]
            pot_field[j,:, :nghost,:] = pot_field[j,:,-2*nghost:-nghost,:]
            pot_field[j,:,-nghost:,:] = pot_field[j,:, nghost: 2*nghost,:]
            pot_field[j,:,:, :nghost] = pot_field[j,:,:,-2*nghost:-nghost]
            pot_field[j,:,:,-nghost:] = pot_field[j,:,:, nghost: 2*nghost]
        if rot:
            rot_field[j, :nghost,:,:] = rot_field[j,-2*nghost:-nghost,:,:]
            rot_field[j,-nghost:,:,:] = rot_field[j, nghost: 2*nghost,:,:]
            rot_field[j,:, :nghost,:] = rot_field[j,:,-2*nghost:-nghost,:]
            rot_field[j,:,-nghost:,:] = rot_field[j,:, nghost: 2*nghost,:]
            rot_field[j,:,:, :nghost] = rot_field[j,:,:,-2*nghost:-nghost]
            rot_field[j,:,:,-nghost:] = rot_field[j,:,:, nghost: 2*nghost]
    #compare internal energy of original and sum of decomposed vectors
    if pot:
        pot2 = dot2(pot_field)[nghost:-nghost,nghost:-nghost,nghost:-nghost]
    if rot:
        rot2 = dot2(rot_field)[nghost:-nghost,nghost:-nghost,nghost:-nghost]
    field2 = dot2(field)[nghost:-nghost,nghost:-nghost,nghost:-nghost]
    if len(field_scalar) > 0:
        #compare kinetic energy of original and sum of decomposed vectors
        field2 *= field_scalar[nghost:-nghost,nghost:-nghost,nghost:-nghost]
        if pot:
            pot2 *= field_scalar[nghost:-nghost,nghost:-nghost,nghost:-nghost]
        if rot:
            rot2 *= field_scalar[nghost:-nghost,nghost:-nghost,nghost:-nghost]
    if rot and not pot:
        if not quiet:
            print(
            'mean total field energy {} mean rotational energy {}'.format(
            np.mean(field2), np.mean(rot2)))
    elif pot and not rot:
        if not quiet:
            print(
            'mean total field energy {} mean irrotational energy {}'.format(
            np.mean(field2), np.mean(pot2)))
    elif rot and pot: 
        if not quiet:
            print(
            'mean total field energy {} mean summed component energy {}'.format(
            np.mean(field2), np.mean(rot2+pot2)))
    #check div and curl approximate/equal zero
    if pot:
        if not quiet:
            print('Max {} and mean {} of abs(curl(pot field))'.format(
              max(grid.dx, grid.dy, grid.dz)*amp_field_1*
              np.sqrt(dot2(curl(pot_field, grid.dx, grid.dy, grid.dz))
                 )[nghost:-nghost,nghost:-nghost,nghost:-nghost].max(),
              max(grid.dx, grid.dy, grid.dz)*amp_field_1*
              np.sqrt(dot2(curl(pot_field, grid.dx, grid.dy, grid.dz))
                 )[nghost:-nghost,nghost:-nghost,nghost:-nghost].mean()))
    if rot:
        if not quiet:
            print('Max {} and mean {} of abs(div(rot field))'.format(
              max(grid.dx, grid.dy, grid.dz)*amp_field_1*
              np.abs(div(rot_field, grid.dx, grid.dy, grid.dz)
             )[nghost:-nghost,nghost:-nghost,nghost:-nghost].max(),
              max(grid.dx, grid.dy, grid.dz)*amp_field_1*
              np.abs(div(rot_field, grid.dx, grid.dy, grid.dz)
             )[nghost:-nghost,nghost:-nghost,nghost:-nghost].mean()))
    if rot and not pot:
        ret_opt = rot_field
    elif pot and not rot:
        ret_opt = pot_field
    elif rot and pot: 
        ret_opt = [rot_field, pot_field]
    else:
        print('pot and/or rot must be True, returning ones')
        ret_opt = np.ones_like(tot_field)
    return ret_opt
