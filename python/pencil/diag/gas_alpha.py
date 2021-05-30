def gas_alpha(sim=False, t_range=[0, -1], OVERWRITE=False):
    """ This calculates the global alpha-value from ts object via
    		alphah=(2./3.)*(ts3.uxuym)/(ts3.rhom*ts3.csm**2)
        for the lastNentries of the time series.
        Input:
            t_range         use this time range of the timeseries
            OVERWRITE       overwrite alpha file in sim.pc_datadir/alpha_<lastNentries>

        return:
            alpha dictionary
    """

    from pencil import get_sim
    from pencil.sim import sim
    from pencil import io
    from os.path import exists, join
    import numpy as np

    def empirical_std_deviation(x):
        """
        (Geschaetzte) Streuung der Messwerte x um den (unbekannten) wahren Wert
        (estimated) stray of the measurments x around the (unknown) true value

        s(x) = SQRT( 1./(M-1) * SUM( (x-<x>)**2 ) )"""
        import numpy as np
        x = np.array(x)
        M = np.size(x)
        xm = np.mean(x)

        #return np.sqrt(1./(M-1.)*np.sum((x-xm)**2))
        return np.sqrt( M/(M-1.) * ( (1./M*np.sum(x**2)) - xm**2 ) )


    def std_deviation_of_mean_value(x):
        """
        Messunsicherheit des Mittelwertes <x>
        Empirical standarddeviation of the arithmetical mean value

        s(<x>) = s(x)/SQRT( M ) = SQRT( 1./(M-1) * SUM( (x-<x>)**2 ) ) / SQRT( M )"""

        import numpy as np
        x = np.array(x)
        M = np.size(x)

        if M == 1: return 0

        return empirical_std_deviation(x)/np.sqrt(M)


    if type(sim) == sim.__Simulation__:
        SIM = sim
    else:
        SIM = get_sim()

    filename = 'alpha_'+str(t_range[0])+'_'+str(t_range[1])

    ## skip if nothing is new
    if not OVERWRITE and io.exists(name=filename, sim=SIM):
        print('~ Alpha for "'+SIM.name+'" already exists. Loading file...')
        return io.load(name=filename, sim=SIM)

    print('~ Calculating alpha for "'+SIM.name+'" in "'+SIM.path+'"')

    ## import time series object
    try:
      print('~ reading time series..')
      ts = SIM.get_ts()
    except:
      print('! ERROR: Couldnt read time series!')
      return False

    ## correct if csm quantity is not correctly exported
    try:
      csm = ts.csm
      if csm[0] == 0: csm = 1.
    except:
      print('? WARNING: Couldnt find >csm< in time series, set it to csm = 1. This may be incorrect!')
      csm = 1.

    if t_range[1] == -1: t_range[1] = ts.t[-1]
    id_min = np.argmin(np.abs(ts.t-t_range[0]))
    id_max = np.argmin(np.abs(ts.t-t_range[1]))

    alpha_dict = {}

    ## calculate alpha
    print('~ calculating alpha, its mean value and standard deviation for '+str(t_range))
    alpha_dict['alpha'] = -(2./3.)*(ts.uxuym)/(ts.rhom*csm**2)
    alpha_dict['alpha_mean'] = np.mean(alpha_dict['alpha'][id_min:id_max])
    alpha_dict['alpha_stddev'] = std_deviation_of_mean_value(alpha_dict['alpha'][id_min:id_max])
    print('~ alpha_mean = '+str(alpha_dict['alpha_mean'])+' and alpha_stddev = '+str(alpha_dict['alpha_stddev']))

    ## calculate alpha minus mean_flux
    print('~ calculating alpha minus mean_flux (alpha_mmf), its mean value and standard deviation')
    alpha_dict['alpha_mmf'] = -(2./3.)*(ts.uxuym-ts.uxm*ts.uym)/(ts.rhom*csm**2)
    alpha_dict['alpha_mmf_mean'] = np.mean(alpha_dict['alpha_mmf'][id_min:id_max])
    alpha_dict['alpha_mmf_stddev'] = std_deviation_of_mean_value(alpha_dict['alpha_mmf'][id_min:id_max])
    print('~ alpha_mmf_mean = '+str(alpha_dict['alpha_mmf_mean'])+' and alpha_mmf_stddev = '+str(alpha_dict['alpha_mmf_stddev']))

    import math
    for v in alpha_dict.values():
        if type(v) == type(1.1) and math.isnan(v):
            io.debug_breakpoint()

    ## save alpha in plk
    print('~ saving alpha values in '+SIM.pc_datadir+'/'+filename)
    io.save(alpha_dict, filename, folder=SIM.pc_datadir)

    return alpha_dict
