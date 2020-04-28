#!/usr/bin/env python
#coding: utf-8
# Helmholz.py
# Written by Fred Gent (fred.gent.ncl@gmail.com)
"""
Perform a Helmholz decomposition on a vector field returning a pair of 
orthogonal vectors with zero divergence and zero curl respectively.
"""

import numpy as np
from ..math import dot, dot2, cross
from ..math.derivatives import div, curl, grad

def helmholz_fft(field, grid, params, nghost=3, nonperi_bc=None):
    """
    Creates the decomposition vector pair for the supplied vector field.

    call signature:

    helmholz_fft(field, grid, params, )
 
    Keyword arguments:

     *field*:
       Vector field of dimension [3,mz,my,mx], which is decomposed.

     *grid*:
       Grid object with grid spacing dx, dy, dz.

     *params*:
       Simulation Params object with domain dimensions Lx, Ly and Lz.

    """
    nz,ny,nx = field[:,nghost:-nghost,nghost:-nghost,nghost:-nghost].shape[-3],\
               field[:,nghost:-nghost,nghost:-nghost,nghost:-nghost].shape[-2],\
               field[:,nghost:-nghost,nghost:-nghost,nghost:-nghost].shape[-1]
    nz2,ny2,nx2 = int(nz/2),int(ny/2),int(nx/2)
    #derive wavenumbers k scaled to dimension of simulation domain
    kk = np.empty(shape=
                  field[:,nghost:-nghost,nghost:-nghost,nghost:-nghost].shape)
    k0 = np.arange(nx); k0[nx2:] = -k0[nx2-1::-1]-1
    k1 = np.arange(ny); k1[ny2:] = -k1[ny2-1::-1]-1
    k2 = np.arange(nz); k2[nz2:] = -k2[nz2-1::-1]-1
    for j in range(0,k0.size):
        kk[0,:,:,j] = k0[j]*2*np.pi/params.lxyz[0]
    for j in range(0,k1.size):
        kk[1,:,j,:] = k1[j]*2*np.pi/params.lxyz[1]
    for j in range(0,k2.size):
        kk[2,j,:,:] = k2[j]*2*np.pi/params.lxyz[2]
    knorm = dot2(kk)
    #apply fast Fourier transform to the vector field
    kfield = np.empty(shape=kk.shape, dtype=complex)
    for j in range(0,3):
        kfield[j] = np.fft.fftn(
                          field[j,nghost:-nghost,nghost:-nghost,nghost:-nghost])
    #avoid division by zero
    knorm[np.where(knorm==0)] = 1.
    #reverse fft to obtain the scalar potential
    pfield = -1j*dot(kk,kfield)/knorm
    scalar_pot = np.zeros_like(field[0])
    scalar_pot[nghost:-nghost,nghost:-nghost,nghost:-nghost] =\
                                                       np.fft.ifftn(pfield).real
    #reverse fft to obtain the vector potential
    rfield =  1j*cross(kk,kfield)/knorm
    vector_pot = np.zeros_like(field[:])
    for j in range(0,3):
        vector_pot[j,nghost:-nghost,nghost:-nghost,nghost:-nghost] =\
                                                    np.fft.ifftn(rfield[j]).real
    #apply the periodic boundary conditions for the ghost zones:
    scalar_pot[:nghost] = scalar_pot[-2*nghost:-nghost]
    scalar_pot[-nghost:] = scalar_pot[nghost:2*nghost]
    scalar_pot[:,:nghost] = scalar_pot[:,-2*nghost:-nghost]
    scalar_pot[:,-nghost:] = scalar_pot[:,nghost:2*nghost]
    scalar_pot[:,:,:nghost] = scalar_pot[:,:,-2*nghost:-nghost]
    scalar_pot[:,:,-nghost:] = scalar_pot[:,:,nghost:2*nghost]
    vector_pot[:,:nghost] = vector_pot[:,-2*nghost:-nghost]
    vector_pot[:,-nghost:] = vector_pot[:,nghost:2*nghost]
    vector_pot[:,:,:nghost] = vector_pot[:,:,-2*nghost:-nghost]
    vector_pot[:,:,-nghost:] = vector_pot[:,:,nghost:2*nghost]
    vector_pot[:,:,:,:nghost] = vector_pot[:,:,:,-2*nghost:-nghost]
    vector_pot[:,:,:,-nghost:] = vector_pot[:,:,:,nghost:2*nghost]
    if nonperi_bc:
        print('Please implement new nonperi_bc not yet implemented.\n',
              'Applying periodic boundary conditions for now.')
    pot_field = np.zeros_like(field)
    rot_field = np.zeros_like(field)
    pot_field[:,nghost:-nghost,nghost:-nghost,nghost:-nghost] =\
       grad(scalar_pot, grid.dx, grid.dy, grid.dz
                               )[:,nghost:-nghost,nghost:-nghost,nghost:-nghost]
    rot_field[:,nghost:-nghost,nghost:-nghost,nghost:-nghost] =\
       curl(vector_pot, grid.dx, grid.dy, grid.dz
                               )[:,nghost:-nghost,nghost:-nghost,nghost:-nghost]
    #apply the periodic boundary conditions for the ghost zones:
    for j in range(0,3):
        pot_field[j, :nghost,:,:] = pot_field[j,-2*nghost:-nghost,:,:]
        pot_field[j,-nghost:,:,:] = pot_field[j, nghost: 2*nghost,:,:]
        pot_field[j,:, :nghost,:] = pot_field[j,:,-2*nghost:-nghost,:]
        pot_field[j,:,-nghost:,:] = pot_field[j,:, nghost: 2*nghost,:]
        pot_field[j,:,:, :nghost] = pot_field[j,:,:,-2*nghost:-nghost]
        pot_field[j,:,:,-nghost:] = pot_field[j,:,:, nghost: 2*nghost]
        rot_field[j, :nghost,:,:] = rot_field[j,-2*nghost:-nghost,:,:]
        rot_field[j,-nghost:,:,:] = rot_field[j, nghost: 2*nghost,:,:]
        rot_field[j,:, :nghost,:] = rot_field[j,:,-2*nghost:-nghost,:]
        rot_field[j,:,-nghost:,:] = rot_field[j,:, nghost: 2*nghost,:]
        rot_field[j,:,:, :nghost] = rot_field[j,:,:,-2*nghost:-nghost]
        rot_field[j,:,:,-nghost:] = rot_field[j,:,:, nghost: 2*nghost]
    if nonperi_bc:
        print('Please implement new nonperi_bc not yet implemented.\n',
              'Applying periodic boundary conditions for now.')
    #check div and curl approximate/equal zero
    print('Max {} and mean {} of abs(curl(pot field))'.format(np.abs(
          curl(pot_field, grid.dx, grid.dy, grid.dz)
              )[:,nghost:-nghost,nghost:-nghost,nghost:-nghost].max(), np.abs(
          curl(pot_field, grid.dx, grid.dy, grid.dz)
              )[:,nghost:-nghost,nghost:-nghost,nghost:-nghost].mean()))
    print('Max {} and mean {} of abs(div(rot field))'.format(np.abs(
          div(rot_field, grid.dx, grid.dy, grid.dz)
             )[nghost:-nghost,nghost:-nghost,nghost:-nghost].max(), np.abs(
          div(rot_field, grid.dx, grid.dy, grid.dz)
             )[nghost:-nghost,nghost:-nghost,nghost:-nghost].mean()))
    return rot_field, pot_field

