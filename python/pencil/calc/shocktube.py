# shocktube.py
#
# Construct Sod solution to shock tube test from simulation
"""
Contains the classes containing the SOD analytic data
"""

import numpy as np
import matplotlib.pyplot as plt
from pencil.read import param

def sod(*args, **kwargs):
    """
    Solve Sod shock tube for given parameters and resolution 

    Signature:

    calc_shocktube(xarr, time, par=list(), lreference=False, DEBUG=False,
                   lplot=False, itplot=0, magic=['ee','tt','Ms'])


    Parameters
    ----------
    *xarr*: Coordinate vector (the initial discontinuity is always at x=0)

    *time*: Time after membrane snapped.
  
    *par*: Param object 
  
    *lreference*: Use default parameters for Sod's reference problem.

    *DEBUG*: Flag to switch on output

    *lplot*: Plot first snapshot profiles

    *itplot*: Iteration index of snaphot to plot

    *magic*: Optional profiles to include in Sod object.
 
    Returns
    -------
    Class containing coordinates and shock profiles

    Notes
    ----- 
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
    Examples
    --------
    >>> shock = pc.calc.sod()
    >>> shock.keys()
    t
    x
    rho
    ux
    pp
    ee
    tt
    Ms
    """

    sod_tmp = SodShock()
    sod_tmp.calc_shocktube(*args, **kwargs)
    return sod_tmp

class SodShock(object):
    """
    SodShock -- holds shock profiles and coordinate arrays
    """

    def __init__(self):
       """
       Fill members with default values
       """

       self.t = np.array([])

    def keys(self):
        for i in self.__dict__.keys():
            print(i)

    def calc_shocktube(self, xarr, time, par=list(), lreference=False,
                       DEBUG=False, lplot=False, itplot=0, 
                       magic=['ee','tt','Ms']
                      ):
        """
        *xarr*: Coordinate vector (the initial discontinuity is always at x=0)

        *time*: Time after membrane snapped.

        *par*: Param object 
  
        *lreference*: Use default parameters for Sod's reference problem.

        *DEBUG*: Flag to switch on output

        *lplot*: Plot first snapshot profiles

        *itplot*: Iteration index of snaphot to plot

        *magic*: Optional profiles to include in Sod object.
 
        """

        #Apply, update and append parameters from simulation directory
        if isinstance(par, list):
            par = param()
        if lreference:
            print("lreference is True: running Sod's reference problem")
            par.uu_right = 0.
            par.uu_left = 0.
            par.rho_right = [0.125,]
            par.rho_left = [1.,]
            par.gamma = 1.4
            par.ss_right = np.log(0.1)/par.gamma - np.log(0.125) #pressure=0.1
            par.ss_left = np.log(3.0)/par.gamma - np.log(1.0) #pressure=3.0
        cv1 = par.gamma/par.cp
        cp1 = 1./par.cp
        cv = par.cp/par.gamma
        gamma_m1 = par.gamma-1.
        gamma1 = 1./par.gamma
        cpmcv1 = 1./(par.cp-cv)
        if par.gamma == 1.0:
            lnTT0 = np.log(par.cs0**2/par.cp)
        else:
            lnTT0 = np.log(par.cs0**2/(par.cp*gamma_m1))
        lnrho0 = np.log(par.rho0)
        lnTT_l = lnTT0+cv1*par.ss_left +gamma_m1*(np.log(par.rho_left[0] )-lnrho0)
        lnTT_r = lnTT0+cv1*par.ss_right+gamma_m1*(np.log(par.rho_right[0])-lnrho0)
        par.__setattr__('pp_left',
                         (par.cp-cv)*np.exp(lnTT_l+np.log(par.rho_left[0] )))
        par.__setattr__('pp_right',
                         (par.cp-cv)*np.exp(lnTT_r+np.log(par.rho_right[0])))
        ## Warn about imperfections:
        if not par.uu_left == 0. or not par.uu_right == 0.:
            print("Case initially not at rest not yet implemented"+
                  " -- results not valid")
        if DEBUG:
            for key in ['pp_left','pp_right']:
                print(key, par.__getattribute__(key))

        csl=np.sqrt(par.gamma*par.pp_left/par.rho_left[0]) # left sound speed
        #    iteratively find p3/pl
        p3=par.pp_left*(par.pp_right/par.pp_left)**0.2     # initial guess
        for i in range(1,21):
            u3 = csl*2/gamma_m1*(1-(p3/par.pp_left)**(gamma_m1/2/par.gamma))
            p3 = par.pp_right + (u3-par.uu_right)*np.sqrt(
                 par.rho_right[0]/2*((par.gamma+1)*p3+gamma_m1*par.pp_right))
            if DEBUG:
                print("p3/pl {}, u3 {}".format(p3/par.pp_left, u3))
    
        rho3 = par.rho_left[0]*(p3/par.pp_left)**(1./par.gamma)
        cs3  = np.sqrt(par.gamma*p3/rho3)
    
        p4 = p3
        u4 = u3
        # velocity of shock front
        us = par.uu_right + (par.pp_right-p4)/(
                             par.uu_right-u4)/par.rho_right[0]
        rho4 = -(par.pp_right-p4)/(par.uu_right-u4)/(u4-us)
        cs4 = np.sqrt(par.gamma*p4/rho4)

        if not isinstance(time, np.ndarray):
            time = np.array(time)
        if not isinstance(xarr, np.ndarray):
            xarr = np.array(xarr)
        print('xarr shape: {}'.format(xarr.shape))
        # declare fields
        ux = np.empty([time.size, xarr.size])
        pp = ux.copy()
        rh = ux.copy()
        for mag in magic:
           if 'tt' == mag:
              tt = ux.copy()
           elif 'ee' == mag:
              ee = ux.copy()
           elif 'Ms' == mag:
              Ms = ux.copy()
           else:
              print('Please implement calculation of {}'.format(mag))
              magic.remove(mag)

        #iterate solution for each time
        for it in range(time.size):
            ##   positions of separating faces
            x1 = (par.uu_left-csl)*time[it]
            if DEBUG:
                print("x1 {}, csl {}, par.uu_left {}".format(
                       x1,    csl,    par.uu_left))
            x2 = (u3-cs3)*time[it]
            x3 = u4*time[it]
            x4 = us*time[it]
    
            ##   calculate profiles
            left  = np.where(xarr <= x1)[0]
            if len(left) > 0:
                # expansion region
                reg2 = np.where(xarr[left[-1]+1:] < x2)[0]+left[-1]+1
            else:
                reg2 = np.where(xarr < x2)[0]
            if len(reg2) > 0:
                reg3 = np.where(xarr[reg2[-1]+1:] < x3)[0]+reg2[-1]+1
            else:
                if len(left) > 0:
                    reg3 = np.where(xarr[left[-1]+1:] < x3)[0]+left[-1]+1
                else:
                    reg3 = np.where(xarr < x3)[0]
            if len(reg3) > 0:
                reg4 = np.where(xarr[reg3[-1]+1:] < x4)[0]+reg3[-1]+1
            else:
                if len(reg2) > 0:
                    reg4 = np.where(xarr[reg2[-1]+1:] < x4)[0]+reg2[-1]+1
                else:
                    if len(left) > 0:
                        reg4 = np.where(xarr[left[-1]+1:] < x4)[0]+left[-1]+1
                    else:
                        reg4 = np.where(xarr < x4)[0]
            right = np.where(xarr > x4)
    
            if len(left) > 0:
                ux[it,left] = par.uu_left
                pp[it,left] = par.pp_left
                rh[it,left] = par.rho_left[0]
    
            if len(reg2) > 0:
                ux[it,reg2] = 2/(par.gamma+1)*(csl+xarr[reg2]/time[it]+\
                              gamma_m1/2*par.uu_right)
                pp[it,reg2] = par.pp_left*(1-gamma_m1/2*ux[it,reg2]/csl)**(
                              2*par.gamma/gamma_m1)
                rh[it,reg2] = par.rho_left[0]*(1-gamma_m1/2*ux[it,reg2]/
                               csl)**(2/gamma_m1)

            if len(reg3) > 0:
                ux[it,reg3] = u3
                pp[it,reg3] = p3
                rh[it,reg3] = rho3
    
            if len(reg4) > 0:
                ux[it,reg4] = u4
                pp[it,reg4] = p4
                rh[it,reg4] = rho4
    
            if len(right) > 0:
                ux[it,right] = par.uu_right
                pp[it,right] = par.pp_right
                rh[it,right] = par.rho_right[0]
            if 'tt' in magic:
                tt[it] = pp[it]*cpmcv1/rh[it]
            if 'ee' in magic:
                if not 'tt' in magic:
                    ee[it] = pp[it]/(gamma_m1*rh[it])
                else:
                    ee[it] = cv*tt[it]
            if 'Ms' in magic:
                if not 'tt' in magic:
                    Ms[it] = ux[it]*np.sqrt(gamma1*rh[it]/pp[it])
                else:
                    Ms[it] = ux[it]/np.sqrt(par.cp*gamma_m1*tt[it])
        setattr(self, 't', time)
        setattr(self, 'x', xarr)
        setattr(self, 'rho', rh)
        setattr(self, 'ux', ux)
        setattr(self, 'pp', pp)
        if 'ee' in magic:
            setattr(self, 'ee', ee)
        if 'tt' in magic:
            setattr(self, 'tt', tt)
        if 'Ms' in magic:
            setattr(self, 'Ms', Ms)


        if lplot:
            fig,(ax1,ax2,ax3) = plt.subplots(3,1,sharex=True)
            ax1.semilogy( xarr, pp[itplot])
            ax1.set(ylabel=r'$p$')
            ax2.semilogy( xarr, rh[itplot])
            ax2.set(ylabel=r'$\rho$')
            ax3.plot( xarr, ux[itplot])
            ax3.set(ylabel=r'$u$')
            ax1.set_title(r'$t={:.3e}$'.format(time[itplot]))
            plt.tight_layout()
            plt.show()
            if len(magic)>0:
                fig, ax = plt.subplots(3,1,sharex=True)
                for ip in range(len(magic)):
                    if magic[ip]=='tt':
                        ylabel = r'$T$'
                    elif magic[ip]=='ee':
                        ylabel = r'$e$'
                    else:
                        ylabel = r'${}$'.format(magic[ip])
                    ax[ip].plot(xarr, self.__getattribute__(magic[ip])[itplot])
                    ax[ip].set(ylabel=ylabel)
                ax[0].set_title(r'$t={:.3e}$'.format(time[itplot]))
                plt.tight_layout()
                plt.show()
    
        if DEBUG:
            print( 'u3={}, u4 ={}'.format( u3, u4))
            print( 'p3={}, p4 ={}'.format( p3, p4))
            print( 'rho4 ={}'.format( rho4))
            print( 'rho3 ={}'.format( rho3))
            print( 'V1 ={}'.format( par.uu_left-csl))
            print( 'V2 ={}'.format( u4-cs3))
            print( 'V3 ={}'.format( u4))
            print( 'V4 ={}'.format( us))
    
        ## End of file shocktube.py
