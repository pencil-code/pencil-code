import numpy as np
from scipy.spatial import cKDTree
from scipy.optimize import curve_fit

def dim_fit(xmin_list, nbins=10): 
    obs, bin_edge = np.histogram(xmin_list, bins=nbins, density=True)
    obs_max = np.amax(obs)
    ip_max = np.where(obs == obs_max)
    ip_max = np.mean(ip_max[0])
    ip_max = int(0.33 * ip_max)
    obs_cut = obs[0:ip_max]
    nbins_fit = int(150.0/len(obs_cut)*nbins)
    print("dim_fit: nbins,nbins_fit=", nbins,nbins_fit)

    obs, bin_edge = np.histogram(xmin_list, bins=nbins_fit, density=True)
    xp_mid = np.array([(bin_edge[j] + bin_edge[j+1])/2 for j in range(len(obs))])
    obs_max = np.amax(obs)
    ip_max = np.where(obs == obs_max)
    ip_max = np.mean(ip_max[0])
    ip_max = int(0.33*ip_max)
    xp_fit = xp_mid[0:ip_max]
    obs_fit = obs[0:ip_max]

    # Curve fit with Lambda operator to power-law tail:
    f_D = lambda xp_fit, a, D: a*D*xp_fit**(D - 1)*np.exp(-1.0*a*xp_fit**D)
    popt, pcov = curve_fit(f_D, xp_fit, obs_fit)

    return popt[1]

def pd2_import_CKD(ipos, ipos3d, xp, yp, zp, xmin_list, ngrsizes):
    tree = cKDTree(ipos3d.T)
    Nsamp=shape(ipos)
    print(Nsamp)

    ipos0=ipos[0]
    for isamp in range(0,Nsamp):
      rps=ipos0[isamp]
      rps3d=np.stack((xp,yp,zp))
      xmin_kdt, index = tree.query(rps3d.T,k=2)
      xmin_list.append(xmin_kdt[1]/xdes)

    # Estimate number of bins needed for histogram:
    nbins=int(np.ceil(Nsamp**0.5))

    # Compute NND stats:
    xmin_mean = np.mean(xmin_list)
    xmin_var  = np.var(xmin_list)
    xmin_max  = np.max(xmin_list)
    xmin_min  = np.min(xmin_list)

    # Compute correlation dimension:
    d2=dim_fit(xmin_list, nbins=ngrsizes)

    return d2
