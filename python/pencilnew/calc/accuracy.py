def accuracy(directory='./', exception=[], field='ux', strip=0, varfile='ogvar.dat',direction=1):
    """
    Assessment of accuracy of simulation by comparison
    Computes the two-norm of velocity fields along strip or in the entire
    2D domain

    call signature:
      accuracy(arrays, strip=0)

    Keyword arguments:

    *directory*
      where to loop over different runs

    *exeptions*
      list of runs or other directories to skip

    *field*
      variable used in accuracy assessment

    *strip*:
      index for strip along coordinate 

    *varfile*:
      name of varfile to read from each sim

    Returns
      array of two-norms where the larges array is used as base 
    """

    import numpy as np
    import os as os
    from pencilnew import read

    # Find directories to include in accuracy assessment
    dirs = []
    for item in os.listdir(directory):
        if os.path.isdir(os.path.join(directory, item)) and item not in exception:
            dirs.append(item)

    # Collect flow data
    sims = []    
    #sims = np.empty([len(dirs)])
    for runs in dirs:
        os.chdir(runs)
        print('Read varfile from:',runs)
        if(varfile=='ogvar.dat'):
            sims.append(read.ogvar())
        else:
            sims.append(read.var(varfile))
        os.chdir('..')

    tmpsims=sims
    # Sort runs by size of r-direction
    sims= sorted(tmpsims, key=lambda x: x.r.size)

    # Check that increase in size is correct for use in two-norm calculation
    nsim = len(sims)
    nx_min = sims[0].r.size
    ny_min = sims[0].th.size
    for i in range(nsim):
        if((sims[i].r.size-1)%(nx_min-1) != 0):
            print('ERROR: Incorrect size in r-dir')
            print('sim[i].r',sims[i].r.size)
            return 0

        if(sims[i].th.size%ny_min != 0):
            print('ERROR: Incorrect size in th-dir')
            return 0

    # Now we are sure that first coordinate of r and th are the same for all runs
    # Can compute the two-norms for increasing sizes. Use largest size as normalization of error
    twonorm = []
    maxerr = []
    dx = sims[0].dx
    n2_factor = int(dx/sims[-1].dx)
    attribute = getattr(sims[-1],field)
    if(field=='ux' or field=='uy'):
        u2 = attribute[0::n2_factor,0::n2_factor]
    else:
        u2 = attribute[0,0::n2_factor,0::n2_factor]

    #sims[-1].rho
    #r64 = sims[-1].r[0::n2_factor]
    #th64 = sims[-1].th[0::n2_factor]
    strip=int(strip)
    for i in range(nsim-1):
        n2_factor = int(dx/sims[i].dx)
        #r = sims[i].r[0::n2_factor]
        #th = sims[i].th[0::n2_factor]
        attribute = getattr(sims[i],field)
        if(field=='ux' or field=='uy'):
            u1 = attribute[0::n2_factor,0::n2_factor]
        else:
            u1 = attribute[0,0::n2_factor,0::n2_factor]

        #u1 = attribute[0,0::n2_factor,0::n2_factor]
        #u1 = sims[i].rho[0,0::n2_factor,0::n2_factor]
        if(direction==1):
            twonorm.append(twonorm_err(u1[:,strip],u2[:,strip],dx))
            maxerr.append(maxerror(u1[:,strip],u2[:,strip]))
        elif(direction==2):
            twonorm.append(twonorm_err(u1[strip,:],u2[strip,:],dx))
            maxerr.append(maxerror(u1[strip,:],u2[strip,:]))

    print('Two-norm computed for field:',field,', along strip:',strip)
    if(direction==1):
        print('Along r-direction')
    else:
        print('Along th-direction')
    return twonorm, maxerr


def twonorm_err(u1,u2,dx):
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
