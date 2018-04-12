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


def shocktube(x, t, parl=[], parr=[], gamma=1.4,
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
    if len(parl) <= 0:
        print("No parameters specified, running Sod's reference problem")
        ul,ur    = 0.,0.
        pl,pr    = 3.,0.1
        rhol,rhor= 1.,0.125
    else:
        ul  ,ur  = parl[0], parr[0]
        pl  ,pr  = parl[1], parr[1]
        rhol,rhor= parl[2], parr[2]
    ##   Warn about imperfections:
    if (ul != 0) or (ur != 0):
        print("Case ur=0 or ul=0 is not OK yet -- results won't make sense")

    gamm1=gamma-1.

    csl=np.sqrt(gamma*pl/rhol)       # left sound speed

    ##   declare fields:
    u   = x*0.
    p   = x*0.
    rho = x*0.

    #    iteratively find p3/pl

    p3=pl*(pr/pl)**0.2             # initial guess
    for i in range(1,21):
        u3 = csl*2/gamm1*(1-(p3/pl)**(gamm1/2/gamma))
        p3 = pr + (u3-ur)*np.sqrt(rhor/2*((gamma+1)*p3+gamm1*pr))

        if DEBUG:
            print("p3/pl {}, u3 {}".format(p3/pl, u3))

    rho3 = rhol*(p3/pl)**(1./gamma)
    cs3  = np.sqrt(gamma*p3/rho3)

    p4 = p3
    u4 = u3
    us = ur + (pr-p4)/(ur-u4)/rhor  # velocity of shock front
    rho4 = -(pr-p4)/(ur-u4)/(u4-us)
    cs4 = np.sqrt(gamma*p4/rho4)

    ##   positions of separating faces
    x1 = (ul-csl)*t
    x2 = (u3-cs3)*t
    x3 = u4*t
    x4 = us*t

    ##   calculate profiles
    left  = np.where(x <= x1)[0]
    if len(left)>0:
        reg2  = np.where(x[left[-1]+1:] < x2)[0]+left[-1]+1 # expansion region
    else:
        reg2  = np.where(x < x2)[0]
    if len(reg2)>0:
        reg3  = np.where(x[reg2[-1]+1:] < x3)[0]+reg2[-1]+1
    else:
        if len(left)>0:
            reg3  = np.where(x[left[-1]+1:] < x3)[0]+left[-1]+1 # expansion region
        else:
            reg3  = np.where(x < x3)[0]
    if len(reg3)>0:
        reg4  = np.where(x[reg3[-1]+1:] < x4)[0]+reg3[-1]+1
    else:
        if len(reg2)>0:
            reg4  = np.where(x[reg2[-1]+1:] < x4)[0]+reg2[-1]+1
        else:
            if len(left)>0:
                reg4  = np.where(x[left[-1]+1:] < x4)[0]+left[-1]+1 # expansion region
            else:
                reg4  = np.where(x < x4)[0]
    right = np.where(x > x4)

    if len(left)>0:
        u[left] = ul
        p[left] = pl
        rho[left] = rhol

    if len(reg2)>0:                                    # expansion region
        u[reg2]   = 2/(gamma+1)*(csl+x[reg2]/t+gamm1/2*ur)
        p[reg2]   = pl*(1-gamm1/2*u[reg2]/csl)**(2*gamma/gamm1)
        rho[reg2] = rhol*(1-gamm1/2*u[reg2]/csl)**(2/gamm1)

    if len(reg3)>0:
        u[reg3] = u3
        p[reg3] = p3
        rho[reg3] = rho3

    if len(reg4)>0:                            # expansion region
        u[reg4] = u4
        p[reg4] = p4
        rho[reg4] = rho4

    if len(right)>0:
        u[right] = ur
        p[right] = pr
        rho[right] = rhor

    if lplot:
        plt.figure()
        fig,(ax1,ax2,ax3) = plt.subplots(3,1,sharex=True)
        ax1.semilogy( x, p)
        ax1.set(ylabel=r'$p$')
        ax2.semilogy( x, rho)
        ax2.set(ylabel=r'$\rho$')
        ax3.plot( x, u)
        ax3.set(ylabel=r'$u$')
        plt.show()

    if DEBUG:
        print( 'u3={}, u4 ={}'.format( u3, u4))
        print( 'p3={}, p4 ={}'.format( p3, p4))
        print( 'rho4 ={}'.format( rho4))
        print( 'rho3 ={}'.format( rho3))
        print( 'V1 ={}'.format( ul-csl))
        print( 'V2 ={}'.format( u4-cs3))
        print( 'V3 ={}'.format( u4))
        print( 'V4 ={}'.format( us))

    return p, u, rho
    ## End of file shocktube.py
