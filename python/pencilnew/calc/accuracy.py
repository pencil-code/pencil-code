# accuracy.py
#
# Module containing routines used for assessment of formal accuracy of runs
# using the ogrid module. Formal order of accuracy computed by L2-norms
#
# Can be extended to also be used on regular grids, if needed
#
# Author:
# J. Aarnes (jorgenaarnes@gmail.com)
#
def twonorm_accuracy(simulations, field='ux', strip=0, varfile='ogvar.dat',direction='x',noerr=True):
    """
    Assessment of accuracy of simulation:
    Computes the two-norm error of all available simulation, where the simulation
    with the maximum amount of grid points is used as the correct/reference solution.
    E.g., for runs with grid points assessment of accuracy of x-component of velocity
    along the y-direction, for runs with grid points nxgrid = n, 2n, 4n, 8n, compute

    || u_n - u_0 || = dy \sum\limits_{n=0}^n (sqrt(u_n(x_n)-u_0(x_8n)))
    
    for all runs (except for the 8n case, used as reference).

    Requires that the runs have matching grids, that is, grid refined by a factor of 2m,
    and grids adjusted so that the grid point overlap (needs ofset if periodic BC is used).

    call signature:
      twonorm_accuracy(simulations)

    Keyword arguments:

    *simulations*
      array of simulation names to be included in the computations

    *field*
      variable used in accuracy assessment

    *strip*:
      index for strip along coordinate 

    *varfile*:
      name of varfile to read from each sim

    *direction*:
      compute two-norm along 'x' or 'y' direction

    *noerr*:
      set to false if you want to return an array of maximum error along strip, in 
      addition to the two-norm

    Returns
      array of two-norms where the larges array is used as base 
    """

    import numpy as np
    import os as os
    from pencilnew import read
    from pencilnew import sim

    # Find directories to include in accuracy assessment
    sims = []
    for simulation in simulations:
        sims.append(sim.get(simulation,quiet=True))

    # Sort runs by size of grid spacing
    if(direction=='x' or direction=='r'):
        sims = sim.sort(sims,'dx',reverse=True)
    elif(direction=='y' or direction=='th'):
        sims = sim.sort(sims,'dy',reverse=True)

    # From this point we only need the varfile to compute accuracy
    # Need to update dim to ogdim for reading of ogvar to be done
    # correctly from simulation object
    for i,thissim in enumerate(sims):
        sims[i].dim = read.ogdim(datadir=thissim.datadir)
        sims[i] = read.ogvar(sim=thissim,trimall=True)

    # Check that increase in size is correct for use in two-norm calculation
    nsims = len(sims)
    nx_min = sims[0].r.size
    ny_min = sims[0].th.size
    for thissim in sims:
        if((thissim.r.size-1)%(nx_min-1) != 0):
            print('ERROR: Incorrect size in r-dir')
            print('sims.r',thissim.r.size)
            return False

        if(thissim.th.size%ny_min != 0):
            print('ERROR: Incorrect size in th-dir')
            return False

    # Now we are sure that first coordinate of r and th are the same for all runs
    # Can compute the two-norms for increasing sizes. Use largest size as normalization of error
    twonorms = np.zeros(nsims-1)
    maxerrs = np.zeros(nsims-1)
    if(direction=='x' or direction=='r'):
        dh = sims[0].dx
        n2_factor = int(dh/sims[-1].dx)
    elif(direction=='y' or direction=='th'):
        dh = sims[0].dy
        n2_factor = int(dh/sims[-1].dy)

    attribute = getattr(sims[-1],field)
    if(field=='ux' or field=='uy'):
        u2 = attribute[0::n2_factor,0::n2_factor]
    else:
        u2 = attribute[0,0::n2_factor,0::n2_factor]

    strip=int(strip)
    i = 0
    for thissim in sims[:-1]:
        if(direction=='x' or direction=='r'):
            n2_factor = int(dh/thissim.dx)
        elif(direction=='y' or direction=='th'):
            n2_factor = int(dh/thissim.dy)
        attribute = getattr(thissim,field)
        if(field=='ux' or field=='uy'):
            u1 = attribute[0::n2_factor,0::n2_factor]
        else:
            u1 = attribute[0,0::n2_factor,0::n2_factor]
        if(direction=='x' or direction=='r'):
           twonorms[i] = (twonorm(u1[:,strip],u2[:,strip],dh))
           maxerrs[i] = (maxerror(u1[:,strip],u2[:,strip]))
        elif(direction=='y' or direction=='th'):
           twonorms[i] = (twonorm(u1[strip,:],u2[strip,:],dh))
           maxerrs[i] = (maxerror(u1[strip,:],u2[strip,:]))
        i = i + 1

    print('Two-norm computed for field:',field,', along strip:',strip)
    if(direction==1):
        print('Along r-direction')
    else:
        print('Along th-direction')

    if not noerr:
        return twonorms, maxerrs
    else:
        return twonorms

def order_accuracy(simulations=[], nstrips=0, twonorm_arr= [], field='ux', varfile='ogvar.dat',direction='x'):
    """
    Compute an estimate of the order of accuracy, using two-norms where 
    the finest grid is used as reference solution u_0.
    Return array of orders of accuracy, where each order p is for computation along one strip.

    p = log(||u(\Delta x) - u_0|| / ||u(\Delta x/2) - u_0||)/log(2)

    Keyword arguments:

    *simulations*
      array of simulation names to be included in the computations

    *nstrips*
      number of strips to include in twonorm_array
      twonorm is computed along each strip, and the strips must coencide on the each
      grid used in such a grid refinement study

    *twonorm_arr*
      array of computed two-normes, if these are available
      if empty, routine will compute these

    *field*
      variable used in accuracy assessment

    *varfile*:
      name of varfile to read from each sim

    *direction*:
      compute two-norm along 'x' or 'y' direction
    """
    import numpy as np

    if(twonorm_arr == []):
        twonorm_arr = twonorm_array(simulations=simulations,nstrips=nstrips,field=field,
                                    varfile=varfile,direction=direction)

    if(twonorm_arr[0,0]==0):
        print('Remove twonorms for strip=0, as these are zero for surface point velocity')
        twonorm_arr = twonorm_arr[1:,:]

    nstrips, ntwonorms = twonorm_arr.shape
    p_arr = np.empty([nstrips,ntwonorms-1])
    
    for i in range(ntwonorms-1):
        p_arr[:,i] = np.log(twonorm_arr[:,i]/twonorm_arr[:,i+1])/np.log(2)

    return p_arr

    
    return p_arr


def twonorm_array(simulations, nstrips, field='ux', varfile='ogvar.dat',direction='x'):
    """
    Compute twonorm_accuracy along the selected direction, for a number of 'strips'.
    See twonorm_accuracy for details.

    call signature:
        twonorm_array(simulations,nstrips)

    Keyword arguments:

    *simulations*
      array of simulation names to be included in the computations

    *nstrips*
      number of strips to include in twonorm_array
      twonorm is computed along each strip, and the strips must coencide on the each
      grid used in such a grid refinement study

    *field*
      variable used in accuracy assessment

    *varfile*:
      name of varfile to read from each sim

    *direction*:
      compute two-norm along 'x' or 'y' direction
    """
    import numpy as np

    tn_arr = np.empty([nstrips,len(simulations)-1])
    for i in range(nstrips):
        tn_arr[i,:] = (twonorm_accuracy(simulations=simulations,field=field,
                                        strip=i,varfile=varfile,direction=direction))

    return tn_arr



def twonorm(u1,u2,dx):
    """
    Compute the two-norm error of two spatial vectors u1 and u2.
    The distance between grid points used in calculation is dx.
    """
    from numpy import sqrt

    twonorm = sqrt(dx*sum((u1-u2)**2))

    return twonorm

def maxerror(u1,u2):
    """
    Find the larges error between values from u1 and u2 (for the 
    same indices).
    """
    from numpy import subtract

    diff = subtract(u1,u2)
    diff = abs(diff)

    return max(diff)
