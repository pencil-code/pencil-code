"""
   shocktube.py


   Author: fg (fred.gent.ncl@gmail.com)
   Date:   10-Apr-2018

   Description
   Adapted from IDL script shocktube.pro
   Analytical solution of the shocktube problem

   Initial state (pressure jump)

   ------------------+------------------
           pl        |        pr
   ------------------+------------------
                     0

   Evolved state: membrane has snapped, four domains

   -----------------------------------------
         1 (l)   |   2  |   3   |  4 | 5 (r)
   --------------+------+-------+----+------
                 x      x       x    x
                 1      2       3    4

   Velocity
                         *************
                       *
                     *
                   *
   -**************--------------------***---

   Nomenclature
   pl, pr - pressure far left/right from the shock structure
   p2     - pressure in the expansion wave
   p3=p4  - no pressure jump across the contact discontinuity
   p: pressure, T: temperature, rho: density,
   u: flow velocity, cs: sound speed

   gamma (adiabatic index) is assumed to be constant

   The borders between the different domains are x1, x2, x3, x4
   and move at velocities ul-csl, u4-c4,  u2, u3, u4, respectively.

   Warning: This works so far only in the case ul=ur=0
"""
import numpy as np
import matplotlib.pyplot as plt


def calc_shocktube(xarr, time, par=list(), lreference=False, 
              DEBUG=False, lplot=False
             ):
    """
      X          : coordinate vector (the initial discontinuity is always
                   at x=0, so you may want to call this routine in the form
                     shocktube, x-x0, t, p, rho, parl, parr, gamma
      T          : time after membrane snapped
      p, u, rho     : pressure and density at positions X
      PARL, PARR : parameters [u,p,rho] to left and right of membrane
      GAMMA      : adiabatic index

     Default parameters are for Sod's reference problem
    """
    if isinstance(par,list):
        par = pc.read.param() 
    if lreference:
        print("lreference is True: running Sod's reference problem")
        par.uu_right = [0.,]
        par.uu_left = [0.,]
        par.rho_right = [0.125,]
        par.rho_left = [1.,]
        par.gamma = 1.4
        sr = np.log(0.1)/par.gamma - np.log(0.125) 
        sl = np.log(3.0)/par.gamma - np.log(1.0) 
        par.ss_right = [sr,] #pressure=0.1 np.exp(gamma*(ss+lnrho)))
        par.ss_left = [sl,] #pressure=3.0
    if not isinstance(par.ss_right, list):
        par.ss_right = [par.ss_right,]
    if not isinstance(par.uu_left, list):
        par.ss_left = [par.ss_left,]
    par.__setattr__('pp_left',[np.exp(par.gamma*(par.ss_left[0]+
                     np.log(par.rho_left[0]))),])
    par.__setattr__('pp_right',[np.exp(par.gamma*(par.ss_right[0]+
                     np.log(par.rho_right[0]))),])
    if not isinstance(par.uu_right, list):
        par.uu_right = [par.uu_right,]
    if not isinstance(par.uu_left, list):
        par.uu_left = [par.uu_left,]
    ## Warn about imperfections:
    if not par.uu_left[0] == 0. or not par.uu_right[0] == 0.:
        print("Case initially not at rest not yet implemented"+
              " -- results not valid")

    gamm1=par.gamma-1.

    csl=np.sqrt(par.gamma*par.pp_left[0]/par.rho_left[0])       # left sound speed

    ##   declare fields:
    uu  = xarr*0.
    pp  = xarr*0.
    rho = xarr*0.

    #    iteratively find p3/pl

    p3=par.pp_left[0]*(par.pp_right[0]/par.pp_left[0])**0.2             # initial guess
    for i in range(1,21):
        u3 = csl*2/gamm1*(1-(p3/par.pp_left[0])**(gamm1/2/par.gamma))
        p3 = par.pp_right[0] + (u3-par.uu_right[0])*np.sqrt(
             par.rho_right[0]/2*((par.gamma+1)*p3+gamm1*par.pp_right[0]))

        if DEBUG:
            print("p3/pl {}, u3 {}".format(p3/par.pp_left[0], u3))

    rho3 = par.rho_left[0]*(p3/par.pp_left[0])**(1./par.gamma)
    cs3  = np.sqrt(par.gamma*p3/rho3)

    p4 = p3
    u4 = u3
    us = par.uu_right[0] + (par.pp_right[0]-p4)/(par.uu_right[0]-u4)/par.rho_right[0]  # velocity of shock front
    rho4 = -(par.pp_right[0]-p4)/(par.uu_right[0]-u4)/(u4-us)
    cs4 = np.sqrt(par.gamma*p4/rho4)

    ##   positions of separating faces
    x1 = (par.uu_left[0]-csl)*time
    x2 = (u3-cs3)*time
    x3 = u4*time
    x4 = us*time

    ##   calculate profiles
    left  = np.where(xarr <= x1)[0]
    if len(left)>0:
        reg2  = np.where(xarr[left[-1]+1:] < x2)[0]+left[-1]+1 # expansion region
    else:
        reg2  = np.where(xarr < x2)[0]
    if len(reg2)>0:
        reg3  = np.where(xarr[reg2[-1]+1:] < x3)[0]+reg2[-1]+1
    else:
        if len(left)>0:
            reg3  = np.where(xarr[left[-1]+1:] < x3)[0]+left[-1]+1 # expansion region
        else:
            reg3  = np.where(xarr < x3)[0]
    if len(reg3)>0:
        reg4  = np.where(xarr[reg3[-1]+1:] < x4)[0]+reg3[-1]+1
    else:
        if len(reg2)>0:
            reg4  = np.where(xarr[reg2[-1]+1:] < x4)[0]+reg2[-1]+1
        else:
            if len(left)>0:
                reg4  = np.where(xarr[left[-1]+1:] < x4)[0]+left[-1]+1 # expansion region
            else:
                reg4  = np.where(xarr < x4)[0]
    right = np.where(xarr > x4)

    if len(left)>0:
        uu[left] = par.uu_left[0]
        pp[left] = par.pp_left[0]
        rho[left] = par.rho_left[0]

    if len(reg2)>0:                                    # expansion region
        uu[reg2]   = 2/(par.gamma+1)*(csl+xarr[reg2]/time+gamm1/2*par.uu_right[0])
        pp[reg2]   = par.pp_left[0]*(1-gamm1/2*uu[reg2]/csl)**(2*par.gamma/gamm1)
        rho[reg2] = par.rho_left[0]*(1-gamm1/2*uu[reg2]/csl)**(2/gamm1)

    if len(reg3)>0:
        uu[reg3] = u3
        pp[reg3] = p3
        rho[reg3] = rho3

    if len(reg4)>0:                            # expansion region
        uu[reg4] = u4
        pp[reg4] = p4
        rho[reg4] = rho4

    if len(right)>0:
        uu[right] = par.uu_right[0]
        pp[right] = par.pp_right[0]
        rho[right] = par.rho_right[0]

    if lplot:
        plt.figure()
        fig,(ax1,ax2,ax3) = plt.subplots(3,1,sharex=True)
        ax1.semilogy( xarr, pp)
        ax1.set(ylabel=r'$p$')
        ax2.semilogy( xarr, rho)
        ax2.set(ylabel=r'$\rho$')
        ax3.plot( xarr, uu)
        ax3.set(ylabel=r'$u$')
        plt.show()

    if DEBUG:
        print( 'u3={}, u4 ={}'.format( u3, u4))
        print( 'p3={}, p4 ={}'.format( p3, p4))
        print( 'rho4 ={}'.format( rho4))
        print( 'rho3 ={}'.format( rho3))
        print( 'V1 ={}'.format( par.uu_left[0]-csl))
        print( 'V2 ={}'.format( u4-cs3))
        print( 'V3 ={}'.format( u4))
        print( 'V4 ={}'.format( us))

    return pp, uu, rho
    ## End of file shocktube.py
