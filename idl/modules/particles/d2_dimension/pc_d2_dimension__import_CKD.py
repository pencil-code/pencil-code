# This piece of code is after work of Lars Mattsson (Nordita).
#
# The number of features here is reduced compared to the IDL
# version of the code, but the purpose of this piece of software
# is to provide a basis for comparison with the IDL software.
#
# This routine is not used as it is, but is called from
# the IDL tool pc_d2_dimension when the keyword ALLPYTHON is set.
#
# Written by: Lars Mattsson, adapted for IDL by Christer Sandin
#
import numpy as np
from scipy.spatial import cKDTree
from scipy.optimize import curve_fit

def dim_fit(xmil, nbins=10): 
    obs = np.zeros(nbins)
    bin_edge = np.zeros(nbins + 1)
    obs, bin_edge = np.histogram(xmil, bins=nbins, density=True)
    obs_max = np.amax(obs)

    ip_max = np.where(obs == obs_max)
    ip_max = np.mean(ip_max[0])
    ip_max = int(0.33 * ip_max)

    obs_cut = obs[0 : ip_max]
    len_obs_cut = ip_max + 1

    nbins_fit = int(100.0 / len_obs_cut * nbins)

    obs, bin_edge = np.histogram(xmil, bins=nbins_fit, density=True)
    obs_max = np.amax(obs)

    ip_max = np.where(obs == obs_max)
    ip_max = np.mean(ip_max[0])
    ip_max = int(0.33 * ip_max)

    xp_mid = np.array([(bin_edge[j] + bin_edge[j + 1]) / 2 for j in range(len(obs))])
    xp_fit = xp_mid[0 : ip_max]
    obs_fit = obs[0 : ip_max]

    # Curve fit with Lambda operator to power-law tail; free parameters: a, D:
    f_D = lambda xp_fit, a, D: a * D * xp_fit**(D - 1) * np.exp(-1.0 * a * xp_fit**D)
    popt, pcov = curve_fit(f_D, xp_fit, obs_fit)

    return popt[1]


def pd2_import_CKD(i, xp, yp, zp, ngrsizes, xdes, ap0j):
    Nsamp = len(xp)
    pxp = np.array(xp)
    pyp = np.array(yp)
    pzp = np.array(zp)

    # Create KD-tree:
    print('  Tree: create')
    tree = cKDTree(np.stack((pxp,pyp,pzp)).T)
    print('  Tree: created')

    # Loop over each element in the position arrays (inefficient):
    xmil = np.zeros(Nsamp)
    for rps in range(0, Nsamp):
        rps3d = np.stack((pxp[rps], pyp[rps], pzp[rps]))
        xmin_kdt, index = tree.query(rps3d, k=2)
        xmil[rps] = xmin_kdt[1] / xdes

    # Compute correlation dimension:
    nbins = int(np.ceil(Nsamp ** 0.5))   # Estimate
    d2 = dim_fit(xmil, nbins=nbins)

    return str(i) + '\t' + str(ap0j) + \
                    '\t' + str(d2) + \
                    '\t' + str(np.mean(xmil)) + \
                    '\t' + str(np.var(xmil))
