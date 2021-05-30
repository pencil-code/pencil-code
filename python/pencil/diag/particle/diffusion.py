
def diffusion(directions=['x'], trange=[0, -1], sim='.', OVERWRITE=False, quiet=True, jump_distance=0.5,
              use_existing_pstalk_sav=False):
    """
    Calculate particle diffusion via stalked particles: PSTALK
    Therefore, it read PSTALK files from Pencil Code using IDL by
    using the IDL<->Python Bridge. It may be activated manually!

    Generated DiffusionData object will be stored in:
        sim.pc_datadir/particle/diffusion/<direction>_<t_start>:<t_end>
    and be reloaded and returned if already existing and OVERWRITE=False!

    Args:
        - sim:         simulation object, but also simulations object or
                       list of simulations or path of simulation,
                        default ='.'

        - directions   list or directions to calculate diffusion for, e.g. ['x','y','z'],
                        default =['x']

        - trange       range in code time for which diffusion shall be calc.
                        default =[0, -1]

        - jump_distance expressed in L_(x,y,z), default: half the domain size

        - OVERWRITE    previously calculated diffusions shall be overwritten
                        default =False

        - quiet        verbosity, default =False

        - use_existing_pstalk_sav
                        use existing <sim.datadir>/data/pc/tmp/pstalk.sav for speed up

    Returns:
        diffusion object for simulation,
        if function is called for a set of simulations, a dictionary with all
        diffusion objects and simulation-name as key is returned
    """

    from pencil.sim import simulations
    import numpy as np

    sims = simulations(sim).sims

    if np.size(sims) == 1:
        if np.size(directions) == 1:
            return DiffusionData(direction=directions[0],
                                trange=trange,
                                sim=sims[0],
                                OVERWRITE=OVERWRITE,
                                quiet=quiet,
                                use_existing_pstalk_sav=use_existing_pstalk_sav)
        else:
            for ii, direction in enumerate(directions):
                directions[ii] = DiffusionData(direction=direction,
                                                trange=trange,
                                                sim=sims[0],
                                                OVERWRITE=OVERWRITE,
                                                quiet=quiet,
                                                use_existing_pstalk_sav=use_existing_pstalk_sav)
            return directions
    else:
        DiffusionData_Dict = {}
        for sim in sims:
            if sim.hidden: print('?? Skipping hidden simulation: '+sim.name); continue
            DiffusionData_Dict[sim.name] = DiffusionData(directions=directions,
                                             trange=trange,
                                             sim=sim,
                                             OVERWRITE=OVERWRITE,
                                             quiet=quiet,
                                             use_existing_pstalk_sav=use_existing_pstalk_sav)
        return DiffusionData_Dict
    return False


class DiffusionData(object):
    """
    DiffusionData -- holds diffusion data for particles.
    """

    def __init__(self, direction='x', trange=[0, -1], sim='.',
                 OVERWRITE=False, quiet=True, jump_distance=0.5, use_existing_pstalk_sav=False):

        from os.path import join
        from os.path import exists as path_exists
        import numpy as np
        from scipy.stats import linregress

        from pencil.io import load, save, exists
        from pencil.read import pstalk as read_pstalk
        from pencil.math.derivatives import simple_centered as derivatives_centered
        from pencil.backpack import printProgressBar

        # .. no pstalk file is found
        if not path_exists(join(sim.datadir, 'proc0', 'particles_stalker.dat')):
          print('?? WARNING: No particles_stalker.dat found for simulation '+sim.name+'!! Skipping this run!')
          return False

        # check if already existing. if, then load it
        out_path = join(sim.pc_datadir, 'particle', 'diffusion')
        out_name = direction+'_'+str(trange[0])+'_'+str(trange[1])

        # skip if diffusion already exists
        if not OVERWRITE and exists(name=out_name, folder=out_path):
            self_tmp = load(out_name, folder=out_path)
            for key in [a for a in dir(self_tmp) if not a.startswith('__')]:
                setattr(self, key, getattr(self_tmp, key))

        else:
            #### start calculations ####
            print('##')
            print('## Calculating particle diffusion for "'+sim.name+'" in "'+sim.path+'"')

            print('## reading particle stalker file..')
            pstalk = read_pstalk(sim=sim, use_existing_pstalk_sav=use_existing_pstalk_sav, tmin=trange[0], tmax=trange[1])
            grid = sim.grid
            dim = sim.dim

            # get time range as index for pstalk dataset
            argmin = np.abs(pstalk.t-trange[0]).argmin()
            if trange[1] < 0.:
              argmax = pstalk.t.argmax()
            else:
              argmax = np.abs(pstalk.t-trange[1]).argmin()
            time_range = pstalk.t[argmin:argmax+1]

            # modify time range with time_range[0] == 0 by substracting first time entry from time list
            #time_offset = pstalk.t[argmin:argmax][0]
            #time_range = pstalk.t[argmin:argmax]-time_offset

            print('\n## doing the calculation for '+direction)
            L = getattr(grid, 'L'+direction)                      # domain size in direction
            N = getattr(dim, 'n'+direction)                       # number grid cells in direction
            pos_series = getattr(pstalk,direction+'p').T[argmin:argmax+1]     # get position timeseries series for direction for all particles
            N_dt = pos_series.shape[0]                            # get number of time steps available
            N_par = pos_series.shape[1]                           # get number of stalked particles

            travel_distance = 0. * pos_series       # prepare 'travel distance' array for all particles
            mean_travel_dist = np.array(0);         # prepare mean, var and sigma arrays
            variance = np.array(0);
            sigma = np.array(0)

            ## calulate travel distances for each particle with corrections for jumps at the boundary
            pbar = False; Nt = np.size(pos_series)
            for i_t,pos in enumerate(pos_series):
                if i_t == 0: continue                                 # skip first time_step
                pbar = printProgressBar(i_t, Nt, pbar=pbar)

                # calculate the distance dx made in dt for all particles at once
                dx = pos - pos_series[i_t-1]
                travel_distance[i_t] = travel_distance[i_t-1] + dx

                # correct for negative jumps
                jumps = np.where(dx > jump_distance*L)
                travel_distance[i_t][jumps] = travel_distance[i_t][jumps] - L

                # correct for positive jumps
                jumps = np.where(dx < -jump_distance*L)
                travel_distance[i_t][jumps] = travel_distance[i_t][jumps] + L

            # calculate mean, variance and sigma for at 'time' i
            mean_travel_dist = np.mean(travel_distance, axis=1)

            # 1. estimate variance as time series
            variance = np.var(travel_distance, axis=1)
            sigma = np.sqrt(variance)

            # 2. estimate diffusion by using centered derivatives method  !!!
            diff_dvar = 0.5*derivatives_centered(time_range, variance)
            diffusion_mean = np.mean(diff_dvar)
            diffusion_std = np.std(diff_dvar)

            # create diffusion object
            self.diffusion = diffusion_mean
            self.diffusion_error = diffusion_std
            self.travel_distance = travel_distance
            self.mean_travel_distance = mean_travel_dist
            self.timerange = time_range
            self.variance = variance
            self.sigma = sigma
            self.direction = direction

            print('~ diffusion = '+str(self.diffusion))
            print('~ diffusion std = '+str(self.diffusion_error))
            print('~ direction = '+str(self.direction))

            try:
                print('## saving results in' + join(out_path,out_name))
                save(obj=self, name=out_name, folder=out_path)
            except:
                print("!! Unexpected error:", sys.exc_info()[0])
                print("!! Check if you have writing rights.")
                raise
