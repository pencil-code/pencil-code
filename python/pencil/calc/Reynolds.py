#!/usr/bin/env python
# coding: utf-8
#
# Reynolds.py
# Written by Fred Gent (fred.gent.ncl@gmail.com)
#
"""
   Calculate the Reynolds number fields from the momentum and induction
   equations using the ratio of advective to diffusive expressions.
"""

from pencil.math import dot, dot2, cross
from pencil.math.derivatives import div, curl, curl2, grad, del2, del6
import numpy as np

def fluid_reynolds(uu, param, grid, lnrho=list(), shock=list(), nghost=3,
                   lmix=True, quiet=True):
    """
    Computes the fluid Reynolds number from the advective and effective
    viscous expressions in the momentum equation.

    call signature:

    fluid_reynolds(uu, ivisc, grid, rho=None, shock=None, nghost=3)

    Keyword arguments:

     *uu*:
       The velocity field [3,mz,my,mx] from the simulation data

     *param*:
       The Param simulation object with viscosity data information

     *grid*:
       The Grid simulation object

     *lnrho*:
       The log density field if it is non-uniform

     *shock*:
       The shock variable if shock viscosity is applied

     *nghost*:
       The number of ghost zones appropriate to the order of accuracy

     *lmix*:
       Option not to include hyper values when Laplacian values present
    """
    #viscous forces
    th2 = 2./3
    th1 = 1./3
    fvisc = np.zeros_like(uu)
    #molecular viscosity contribution
    ldel2, lshock, lhyper3 = False, False, False
    for ivisc in param.ivisc:
        if not 'shock' in ivisc and not 'hyper' in ivisc\
                                and not '\n' in ivisc:
            ldel2 = True
        if 'shock' in ivisc:
            if not isinstance(shock, list):
                lshock = True
        if 'hyper3' in ivisc:
            lhyper3 = True

    if ldel2:
        if lhyper3:
            lhyper3 = lhyper3==lmix
        del2u = np.zeros_like(uu)
        for j in range(0,3):
            del2u[j] = del2(uu[j],grid.dx,grid.dy,grid.dz,x=grid.x,y=grid.y,
                            coordinate_system=param.coord_system)
            del2u[j,       : nghost, nghost:-nghost, nghost:-nghost] = del2u[j,-2*nghost: -nghost,   nghost: -nghost,   nghost: -nghost]
            del2u[j,-nghost:       , nghost:-nghost, nghost:-nghost] = del2u[j,   nghost:2*nghost,   nghost: -nghost,   nghost: -nghost]
            del2u[j, nghost:-nghost,       : nghost, nghost:-nghost] = del2u[j,   nghost: -nghost,-2*nghost: -nghost,   nghost: -nghost]
            del2u[j, nghost:-nghost,-nghost:       , nghost:-nghost] = del2u[j,   nghost: -nghost,   nghost:2*nghost,   nghost: -nghost]
            del2u[j, nghost:-nghost, nghost:-nghost,       : nghost] = del2u[j,   nghost: -nghost,   nghost: -nghost,-2*nghost: -nghost]
            del2u[j, nghost:-nghost, nghost:-nghost,-nghost:       ] = del2u[j,   nghost: -nghost,   nghost: -nghost,   nghost:2*nghost]
        for ivisc in param.ivisc:
            ivisc = str.strip(ivisc,'\n')
            if 'nu-const' not in ivisc and 'shock' not in ivisc\
                                   and 'hyper' not in ivisc and len(ivisc) > 0:
                print('fluid_reynolds WARNING: '+ivisc+' not implemented\n'+
                'terms may be missing from the standard rate of strain tensor')
        fvisc = fvisc + param.nu*del2u
        del(del2u)
    tmp0 = grad(uu[0],grid.dx,grid.dy,grid.dz,x=grid.x,y=grid.y,
                coordinate_system=param.coord_system)
    for j in range(0,3):
        tmp0[j, :nghost,nghost:-nghost,nghost:-nghost] = tmp0[j,-2*nghost:-nghost,nghost:-nghost,nghost:-nghost]
        tmp0[j,-nghost:,nghost:-nghost,nghost:-nghost] = tmp0[j, nghost: 2*nghost,nghost:-nghost,nghost:-nghost]
        tmp0[j,nghost:-nghost, :nghost,nghost:-nghost] = tmp0[j,nghost:-nghost,-2*nghost:-nghost,nghost:-nghost]
        tmp0[j,nghost:-nghost,-nghost:,nghost:-nghost] = tmp0[j,nghost:-nghost, nghost: 2*nghost,nghost:-nghost]
        tmp0[j,nghost:-nghost,nghost:-nghost, :nghost] = tmp0[j,nghost:-nghost,nghost:-nghost,-2*nghost:-nghost]
        tmp0[j,nghost:-nghost,nghost:-nghost,-nghost:] = tmp0[j,nghost:-nghost,nghost:-nghost, nghost: 2*nghost]
    tmp1 = grad(uu[1],grid.dx,grid.dy,grid.dz,x=grid.x,y=grid.y,
                coordinate_system=param.coord_system)
    for j in range(0,3):
        tmp1[j, :nghost,nghost:-nghost,nghost:-nghost] = tmp1[j,-2*nghost:-nghost,nghost:-nghost,nghost:-nghost]
        tmp1[j,-nghost:,nghost:-nghost,nghost:-nghost] = tmp1[j, nghost: 2*nghost,nghost:-nghost,nghost:-nghost]
        tmp1[j,nghost:-nghost, :nghost,nghost:-nghost] = tmp1[j,nghost:-nghost,-2*nghost:-nghost,nghost:-nghost]
        tmp1[j,nghost:-nghost,-nghost:,nghost:-nghost] = tmp1[j,nghost:-nghost, nghost: 2*nghost,nghost:-nghost]
        tmp1[j,nghost:-nghost,nghost:-nghost, :nghost] = tmp1[j,nghost:-nghost,nghost:-nghost,-2*nghost:-nghost]
        tmp1[j,nghost:-nghost,nghost:-nghost,-nghost:] = tmp1[j,nghost:-nghost,nghost:-nghost, nghost: 2*nghost]
    tmp2 = grad(uu[2],grid.dx,grid.dy,grid.dz,x=grid.x,y=grid.y,
                coordinate_system=param.coord_system)
    for j in range(0,3):
        tmp2[j, :nghost,nghost:-nghost,nghost:-nghost] = tmp2[j,-2*nghost:-nghost,nghost:-nghost,nghost:-nghost]
        tmp2[j,-nghost:,nghost:-nghost,nghost:-nghost] = tmp2[j, nghost: 2*nghost,nghost:-nghost,nghost:-nghost]
        tmp2[j,nghost:-nghost, :nghost,nghost:-nghost] = tmp2[j,nghost:-nghost,-2*nghost:-nghost,nghost:-nghost]
        tmp2[j,nghost:-nghost,-nghost:,nghost:-nghost] = tmp2[j,nghost:-nghost, nghost: 2*nghost,nghost:-nghost]
        tmp2[j,nghost:-nghost,nghost:-nghost, :nghost] = tmp2[j,nghost:-nghost,nghost:-nghost,-2*nghost:-nghost]
        tmp2[j,nghost:-nghost,nghost:-nghost,-nghost:] = tmp2[j,nghost:-nghost,nghost:-nghost, nghost: 2*nghost]
    #effect of compressibility
    if len(lnrho) > 0:
        divu = div(uu,grid.dx,grid.dy,grid.dz,x=grid.x,y=grid.y,
                   coordinate_system=param.coord_system)
        divu[ :nghost,nghost:-nghost,nghost:-nghost] = divu[-2*nghost:-nghost,nghost:-nghost,nghost:-nghost]
        divu[-nghost:,nghost:-nghost,nghost:-nghost] = divu[ nghost: 2*nghost,nghost:-nghost,nghost:-nghost]
        divu[nghost:-nghost, :nghost,nghost:-nghost] = divu[nghost:-nghost,-2*nghost:-nghost,nghost:-nghost]
        divu[nghost:-nghost,-nghost:,nghost:-nghost] = divu[nghost:-nghost, nghost: 2*nghost,nghost:-nghost]
        divu[nghost:-nghost,nghost:-nghost, :nghost] = divu[nghost:-nghost,nghost:-nghost,-2*nghost:-nghost]
        divu[nghost:-nghost,nghost:-nghost,-nghost:] = divu[nghost:-nghost,nghost:-nghost, nghost: 2*nghost]
        gradlnrho = grad(lnrho,grid.dx,grid.dy,grid.dz,x=grid.x,y=grid.y,
                         coordinate_system=param.coord_system)
        for j in range(0,3):
            gradlnrho[j, :nghost,nghost:-nghost,nghost:-nghost] = gradlnrho[j,-2*nghost:-nghost,nghost:-nghost,nghost:-nghost]
            gradlnrho[j,-nghost:,nghost:-nghost,nghost:-nghost] = gradlnrho[j, nghost: 2*nghost,nghost:-nghost,nghost:-nghost]
            gradlnrho[j,nghost:-nghost, :nghost,nghost:-nghost] = gradlnrho[j,nghost:-nghost,-2*nghost:-nghost,nghost:-nghost]
            gradlnrho[j,nghost:-nghost,-nghost:,nghost:-nghost] = gradlnrho[j,nghost:-nghost, nghost: 2*nghost,nghost:-nghost]
            gradlnrho[j,nghost:-nghost,nghost:-nghost, :nghost] = gradlnrho[j,nghost:-nghost,nghost:-nghost,-2*nghost:-nghost]
            gradlnrho[j,nghost:-nghost,nghost:-nghost,-nghost:] = gradlnrho[j,nghost:-nghost,nghost:-nghost, nghost: 2*nghost]
        Sglnrho = np.zeros_like(uu)
        Sglnrho[0] = dot(tmp0,gradlnrho) +\
                        (tmp0[0]+tmp1[0]+tmp2[0]-th2*divu)*gradlnrho[0]
        Sglnrho[1] = dot(tmp1,gradlnrho) +\
                        (tmp0[1]+tmp1[1]+tmp2[1]-th2*divu)*gradlnrho[1]
        Sglnrho[2] = dot(tmp2,gradlnrho) +\
                        (tmp0[2]+tmp1[2]+tmp2[2]-th2*divu)*gradlnrho[2]
        graddivu = grad(divu,grid.dx,grid.dy,grid.dz,x=grid.x,y=grid.y,
                        coordinate_system=param.coord_system)
        for j in range(0,3):
            graddivu[j, :nghost,nghost:-nghost,nghost:-nghost] = graddivu[j,-2*nghost:-nghost,nghost:-nghost,nghost:-nghost]
            graddivu[j,-nghost:,nghost:-nghost,nghost:-nghost] = graddivu[j, nghost: 2*nghost,nghost:-nghost,nghost:-nghost]
            graddivu[j,nghost:-nghost, :nghost,nghost:-nghost] = graddivu[j,nghost:-nghost,-2*nghost:-nghost,nghost:-nghost]
            graddivu[j,nghost:-nghost,-nghost:,nghost:-nghost] = graddivu[j,nghost:-nghost, nghost: 2*nghost,nghost:-nghost]
            graddivu[j,nghost:-nghost,nghost:-nghost, :nghost] = graddivu[j,nghost:-nghost,nghost:-nghost,-2*nghost:-nghost]
            graddivu[j,nghost:-nghost,nghost:-nghost,-nghost:] = graddivu[j,nghost:-nghost,nghost:-nghost, nghost: 2*nghost]
        fvisc = fvisc + param.nu*(th1*graddivu+Sglnrho)
        del(Sglnrho)
    elif param.ldensity:
        print('fluid_reynolds WARNING: no lnrho provided\n'+
              'rate of strain tensor likely incomplete')
    #shock contribution
    if lshock:
        if len(shock) == 0:
            print('fluid_reynolds WARNING: no shock provided\n'+
                  'rate of strain tensor likely incomplete')
        else:
            shock[ :nghost,nghost:-nghost,nghost:-nghost] = shock[-2*nghost:-nghost,nghost:-nghost,nghost:-nghost]
            shock[-nghost:,nghost:-nghost,nghost:-nghost] = shock[ nghost: 2*nghost,nghost:-nghost,nghost:-nghost]
            shock[nghost:-nghost, :nghost,nghost:-nghost] = shock[nghost:-nghost,-2*nghost:-nghost,nghost:-nghost]
            shock[nghost:-nghost,-nghost:,nghost:-nghost] = shock[nghost:-nghost, nghost: 2*nghost,nghost:-nghost]
            shock[nghost:-nghost,nghost:-nghost, :nghost] = shock[nghost:-nghost,nghost:-nghost,-2*nghost:-nghost]
            shock[nghost:-nghost,nghost:-nghost,-nghost:] = shock[nghost:-nghost,nghost:-nghost, nghost: 2*nghost]
            divugradlnrho = np.zeros_like(uu)
            gradshock = grad(shock,grid.dx,grid.dy,grid.dz,x=grid.x,y=grid.y,
                             coordinate_system=param.coord_system)
            for j in range(0,3):
                gradshock[j, :nghost,nghost:-nghost,nghost:-nghost] = gradshock[j,-2*nghost:-nghost,nghost:-nghost,nghost:-nghost]
                gradshock[j,-nghost:,nghost:-nghost,nghost:-nghost] = gradshock[j, nghost: 2*nghost,nghost:-nghost,nghost:-nghost]
                gradshock[j,nghost:-nghost, :nghost,nghost:-nghost] = gradshock[j,nghost:-nghost,-2*nghost:-nghost,nghost:-nghost]
                gradshock[j,nghost:-nghost,-nghost:,nghost:-nghost] = gradshock[j,nghost:-nghost, nghost: 2*nghost,nghost:-nghost]
                gradshock[j,nghost:-nghost,nghost:-nghost, :nghost] = gradshock[j,nghost:-nghost,nghost:-nghost,-2*nghost:-nghost]
                gradshock[j,nghost:-nghost,nghost:-nghost,-nghost:] = gradshock[j,nghost:-nghost,nghost:-nghost, nghost: 2*nghost]
            for j in range(0,3):
                divugradlnrho[j] = param.nu_shock*divu*gradshock[j] +\
                          param.nu_shock*shock*(divu*gradlnrho[j] + graddivu[j])
            del(divu,gradshock,gradlnrho,graddivu)
            fvisc = fvisc + divugradlnrho
            del(divugradlnrho)
    if lhyper3:
        #deluij5 = np.zeros_like([uu,uu,uu])
        #uij5glnrho to be included
        del6u = np.zeros_like(uu)
        for j in range(0,3):
            del6u[j] = del6(uu[j],grid.dx,grid.dy,grid.dz)
            del6u[j, :nghost,nghost:-nghost,nghost:-nghost] = del6u[j,-2*nghost:-nghost,nghost:-nghost,nghost:-nghost]
            del6u[j,-nghost:,nghost:-nghost,nghost:-nghost] = del6u[j, nghost: 2*nghost,nghost:-nghost,nghost:-nghost]
            del6u[j,nghost:-nghost, :nghost,nghost:-nghost] = del6u[j,nghost:-nghost,-2*nghost:-nghost,nghost:-nghost]
            del6u[j,nghost:-nghost,-nghost:,nghost:-nghost] = del6u[j,nghost:-nghost, nghost: 2*nghost,nghost:-nghost]
            del6u[j,nghost:-nghost,nghost:-nghost, :nghost] = del6u[j,nghost:-nghost,nghost:-nghost,-2*nghost:-nghost]
            del6u[j,nghost:-nghost,nghost:-nghost,-nghost:] = del6u[j,nghost:-nghost,nghost:-nghost, nghost: 2*nghost]
            #del6 for non-cartesian tba
            #del6u[j] = del6(uu[j],grid.dx,grid.dy,grid.dz,x=grid.x,y=grid.y,
            #                coordinate_system=param.coord_system)
            fvisc = fvisc + param.nu_hyper3*del6u
        del(del6u)
    fvisc2 = np.sqrt(dot2(fvisc))
    #advective forces
    advec = np.zeros_like(uu)
    advec[0] = dot(uu,tmp0)
    advec[1] = dot(uu,tmp1)
    advec[2] = dot(uu,tmp2)
    del(tmp0,tmp1,tmp2)
    advec2 = np.sqrt(dot2(advec))
    del(advec)
    #avoid division by zero
    if fvisc2.max() > 0:
        fvisc2[np.where(fvisc2==0)] = fvisc2[np.where(fvisc2>0)].min()
        Re = advec2/fvisc2
        #set minimum floor to exclude zero-valued Re
        Re[np.where(Re==0)] = Re[np.where(Re>0)].min()
    else:
        Re = advec2
        print('Re undefined')
    return Re

def magnetic_reynolds(uu, param, grid, aa=list(), bb=list(), jj=list(),
                      nghost=3, lmix=True, quiet=True):
    """
    Computes the magnetic Reynolds number from the advective and effective
    resistive expressions in the induction equation.

    call signature:

    magnetic_reynolds(uu, param, grid, aa=None, bb=None, jj=None, nghost=3):

    Keyword arguments:

     *uu*:
       The velocity field [3,mz,my,mx] from the simulation data

     *param*:
       The Param simulation object with resistivity data information

     *grid*:
       The Grid simulation object

     *aa*:
       The vector potential if bb is not present or hyper diffusion

     *bb*:
       The magnetic field

     *jj*:
       The current density field

     *nghost*:
       The number of ghost zones appropriate to the order of accuracy

     *lmix*:
       Option not to include hyper values when Laplacian values present
    """
    if len(bb) ==0 and len(aa) ==0 and len(jj) ==0:
        print('magnetic_reynolds WARNING: no aa, bb nor jj provided\n'+
              'aa or bb must be provided or aa for only hyper resistivity')
    #resistive force
    lres, lhyper3 = False, False
    for iresi in param.iresistivity:
        iresi = str.strip(iresi,'\n')
        if 'hyper' not in iresi and len(iresi) > 0:
            lres = True
        if 'hyper3' in iresi:
            lhyper3 = True
    fresi = np.zeros_like(uu)
    if lres:
        if lhyper3:
            lhyper3 = lhyper3==lmix
        if len(jj) == 0:
            if len(aa) == 0:
                print('magnetic_reynolds WARNING: calculating jj without aa\n',
                      'provide aa or jj directly for accurate boundary values')
                jj = curl(bb,grid.dx,grid.dy,grid.dz,x=grid.x,y=grid.y,
                            coordinate_system=param.coord_system)
            else:
                jj = curl2(aa,grid.dx,grid.dy,grid.dz,x=grid.x,y=grid.y,
                            coordinate_system=param.coord_system)
            for j in range(0,3):
                jj[j, :nghost,:,:] = jj[j,-2*nghost:-nghost,:,:]
                jj[j,-nghost:,:,:] = jj[j, nghost: 2*nghost,:,:]
                jj[j,:, :nghost,:] = jj[j,:,-2*nghost:-nghost,:]
                jj[j,:,-nghost:,:] = jj[j,:, nghost: 2*nghost,:]
                jj[j,:,:, :nghost] = jj[j,:,:,-2*nghost:-nghost]
                jj[j,:,:,-nghost:] = jj[j,:,:, nghost: 2*nghost]
        fresi = fresi + param.eta*param.mu0*jj
        for iresi in param.iresistivity:
            iresi = str.strip(iresi,'\n')
            if 'eta-const' not in iresi and 'hyper' not in iresi\
                                        and len(iresi) > 0:
                print('magnetic_reynolds WARNING: '+iresi+' not implemented\n'+
                      'terms may be missing from the standard resistive forces')
    if lhyper3:
        if len(aa) == 0:
            print('magnetic_reynolds WARNING: no aa provided\n'+
                  'aa must be provided for hyper resistivity')
            return 1
        else:
            del6a = np.zeros_like(aa)
            for j in range(0,3):
                del6a[j] = del6(aa[j],grid.dx,grid.dy,grid.dz)
                del6a[j, :nghost,:,:] = del6a[j,-2*nghost:-nghost,:,:]
                del6a[j,-nghost:,:,:] = del6a[j, nghost: 2*nghost,:,:]
                del6a[j,:, :nghost,:] = del6a[j,:,-2*nghost:-nghost,:]
                del6a[j,:,-nghost:,:] = del6a[j,:, nghost: 2*nghost,:]
                del6a[j,:,:, :nghost] = del6a[j,:,:,-2*nghost:-nghost]
                del6a[j,:,:,-nghost:] = del6a[j,:,:, nghost: 2*nghost]
                #del6 for non-cartesian tba
                #del6a[j] = del6(aa[j],grid.dx,grid.dy,grid.dz,x=grid.x,y=grid.y,
                #                coordinate_system=param.coord_system)
            #effective at l > 5 grid.dx?
            fresi = fresi + param.eta_hyper3*del6a
            del(del6a)
    fresi2 = np.sqrt(dot2(fresi))
    del(fresi)
    #advective force
    if len(bb) == 0:
        if len(aa) == 0:
            print('magnetic_reynolds WARNING: calculating uu x bb without bb\n',
                  'provide aa or bb directly to proceed')
            return 1
        else:
            bb = curl(aa,grid.dx,grid.dy,grid.dz,x=grid.x,y=grid.y,
                      coordinate_system=param.coord_system)
            for j in range(0,3):
                bb[j, :nghost,:,:] = bb[j,-2*nghost:-nghost,:,:]
                bb[j,-nghost:,:,:] = bb[j, nghost: 2*nghost,:,:]
                bb[j,:, :nghost,:] = bb[j,:,-2*nghost:-nghost,:]
                bb[j,:,-nghost:,:] = bb[j,:, nghost: 2*nghost,:]
                bb[j,:,:, :nghost] = bb[j,:,:,-2*nghost:-nghost]
                bb[j,:,:,-nghost:] = bb[j,:,:, nghost: 2*nghost]
    advec = cross(uu,bb)
    advec2 = np.sqrt(dot2(advec))
    del(advec)
    #avoid division by zero
    if fresi2.max() > 0:
        fresi2[np.where(fresi2==0)] = fresi2[np.where(fresi2>0)].min()
        Rm = advec2/fresi2
        #set minimum floor to exclude zero-valued Rm
        if Rm.max() > 0:
            Rm[np.where(Rm==0)] = Rm[np.where(Rm>0)].min()
        else:
            print('Rm undefined')
    else:
        Rm = advec2
        print('Rm undefined')
    return Rm
