
def dispersion_and_drift(sim=False, OVERWRITE=False, GLOBAL=True, LOCAL=True, use_IDL=False, recalculate_gas_velo_at_particle_pos=False):
    """This calculates the dispersion (sigma) and drift (zeta) locally and globally
    by using the gas_velo_at_particle_pos script and dataset for all particles.

    With sigma = sqrt( 1/N_par * sum_i^N_par( (v_par(i) - <v_par>)^2 ) ) and
    zeta = sqrt( 1/N_par * sum_i^N_par( (v_par(i) - u(xp_i))^2 ) )

    Arg:
	  OVERWRITE:		set to True to overwrite already calculated results
      GLOBAL:         Calculate drift and dispersion globally, i.e. whole simulation domain
      LOCAL:          Calculate drift and dispersion locally, i.e. grid cell wise
      recalculate_gas_velo_at_particle_pos:
                    if the dataset shall be recalcualted
      use_IDL:  use backup solution of IDL script and sav files

      returns True if successfull
    """

    from pencil import get_sim
    from pencil import io
    from pencil import read
    from pencil.diag.particle import gas_velo_at_particle_pos
    from scipy.io import readsav
    from os import listdir
    from os.path import exists, join, dirname
    import numpy as np

    if sim == False:
        sim = get_sim()
        if sim == False:
            print('! ERROR: Specify simulation object!')
            return False
    SIM = sim

    ## calculate gas speed at particle position dataset
    gas_velo_at_particle_pos(OVERWRITE=recalculate_gas_velo_at_particle_pos, use_IDL=use_IDL, sim=sim)

    print('\n##################### Starting the whole calculation process of DISPERSON and DRIFT for '+SIM.name+' #####################')

    ## default and setup
    GASVELO_DESTINATION = 'gas_velo_at_particle_pos'
    GASVELO_DIR = join(SIM.pc_datadir, GASVELO_DESTINATION)


    ## get list of available files
    if use_IDL:
        file_filetype = '.sav'
    else:
        file_filetype = '.pkl'
    files = []
    if exists(GASVELO_DIR):
        files = [i for i in listdir(GASVELO_DIR) if i.startswith(GASVELO_DESTINATION) and i.endswith(file_filetype)]
    if files == []: print('!! ERROR: No calc_gas_speed_at_particle_position-files found for '+SIM.name+'! Use idl script to produce them first!')
    if use_IDL:
        USE_PKL_FILES = False
        files = [i.split('_')[-1].split(file_filetype)[0] for i in files]
        scheme = ''
    else:
        USE_PKL_FILES = True
        files = [i.split('_')[-1].split(file_filetype)[0] for i in files]
        scheme = '_tsc'

    ## calculate global dispersion for all snapshot for which gas_velo_at_particle_pos files are found
    for file_no in files:
      print('## Starting the calculation for DISPERSON and DRIFT for  ### VAR'+str(file_no)+' ###')

      # check if files already exist
      if (not OVERWRITE) and io.exists('sigma_'+file_no, folder=SIM.pc_datadir) and io.exists('zeta_'+file_no, folder=SIM.pc_datadir) and io.exists('sigma_l_'+file_no, folder=join(SIM.pc_datadir, 'sigma_l')) and io.exists('zeta_l_'+file_no, folder=join(SIM.pc_datadir, 'zeta_l')):
        print('## Skipping calculations')
        continue

      ## read sav and var file
      print('## reading gas_velo_at_particle_pos file and VAR')
      if USE_PKL_FILES:
          sav_file	= io.load(join(GASVELO_DIR, GASVELO_DESTINATION+scheme+'_'+file_no+'.pkl'))[GASVELO_DESTINATION]
      else:
          sav_file	= readsav(join(GASVELO_DIR, GASVELO_DESTINATION+scheme+'_'+file_no+'.sav'))[GASVELO_DESTINATION]
      var_file	= read.var(varfile='VAR'+file_no, quiet=True, trimall=True, datadir=SIM.datadir)

      ## get everything ready
      dim = SIM.dim
      pdim = read.pdim(sim=SIM)
      npar = pdim.npar
      npar1 = 1./npar

      ## get first quantities and setup the DATA_SET
      if use_IDL:
          time = sav_file['time'][0]
          DATA_SET = np.core.records.fromarrays(
                            [range(0,npar),
                            sav_file['par_idx'][0][0].astype('int'), 	sav_file['par_idx'][0][1].astype('int'), 	sav_file['par_idx'][0][2].astype('int'),
                            sav_file['par_pos'][0][0], 			sav_file['par_pos'][0][1], 			sav_file['par_pos'][0][2],
                            sav_file['par_velo'][0][0], 			sav_file['par_velo'][0][1], 			sav_file['par_velo'][0][2],
                            sav_file['npar'][0],
                            sav_file['gas_velo'][0][0],			sav_file['gas_velo'][0][1],			sav_file['gas_velo'][0][2],
                            var_file.rho[sav_file['par_idx'][0][2].astype('int'), sav_file['par_idx'][0][1].astype('int'), sav_file['par_idx'][0][0].astype('int')],
                            var_file.rhop[sav_file['par_idx'][0][2].astype('int'), sav_file['par_idx'][0][1].astype('int'), sav_file['par_idx'][0][0].astype('int')]
                            ],
                        names = 'parid, idx, idy, idz, posx, posy, posz, vx, vy, vz, npar, gasv_x, gasv_y, gasv_z, rho, rhop',
                        formats =  'int, int,int,int, float,float,float, float,float,float, int, float,float,float, float, float')
      else:
          time = sav_file['time']
          DATA_SET = np.core.records.fromarrays(
                            [range(0,npar),
                            sav_file['par_idx'][0].astype('int'), 	sav_file['par_idx'][1].astype('int'), 	sav_file['par_idx'][2].astype('int'),
                            sav_file['par_pos'][0], 			sav_file['par_pos'][1], 			sav_file['par_pos'][2],
                            sav_file['par_velo'][0], 			sav_file['par_velo'][1], 			sav_file['par_velo'][2],
                            sav_file['npar'],
                            sav_file['gas_velo'][0],			sav_file['gas_velo'][1],			sav_file['gas_velo'][2],
                            var_file.rho[sav_file['par_idx'][2].astype('int'), sav_file['par_idx'][1].astype('int'), sav_file['par_idx'][0].astype('int')],
                            var_file.rhop[sav_file['par_idx'][2].astype('int'), sav_file['par_idx'][1].astype('int'), sav_file['par_idx'][0].astype('int')]
                            ],
                        names = 'parid, idx, idy, idz, posx, posy, posz, vx, vy, vz, npar, gasv_x, gasv_y, gasv_z, rho, rhop',
                        formats =  'int, int,int,int, float,float,float, float,float,float, int, float,float,float, float, float')

      DATA_SET = np.sort(DATA_SET, order=['idx', 'idy', 'idz'])


      # calculate GLOBAL DISPERSION in x, y and z direction, also the absolute magnitude
      if GLOBAL:
        print('## Calculating GLOBAL DISPERSION values in x,y,z direction and abs value')
        mean_vx = np.mean(DATA_SET['vx']); mean_vy = np.mean(DATA_SET['vy']); mean_vz = np.mean(DATA_SET['vz'])
        SIGMA = {'SIGMA_o_x': np.sqrt(npar1 * np.sum( (DATA_SET['vx'] - mean_vx)**2 ) ),
                'SIGMA_o_y': np.sqrt(npar1 * np.sum( (DATA_SET['vy'] - mean_vy)**2 ) ),
                'SIGMA_o_z': np.sqrt(npar1 * np.sum( (DATA_SET['vz'] - mean_vz)**2 ) ),
                'SIGMA_o': np.sqrt(npar1 * np.sum( (DATA_SET['vx'] - mean_vx)**2 + (DATA_SET['vy'] - mean_vy)**2 + (DATA_SET['vz'] - mean_vz)**2 ) )}


        # calculate GLOBAL DRIFT in x, y and z direction, also the absolute magnitude
        print('## Calculating GLOBAL DRIFT values in x,y,z direction and abs value')
        ZETA = {'ZETA_o_x': np.sqrt(npar1 * np.sum((DATA_SET['vx'] - DATA_SET['gasv_x'])**2)),
                'ZETA_o_y': np.sqrt(npar1 * np.sum((DATA_SET['vy'] - DATA_SET['gasv_y'])**2)),
                'ZETA_o_z': np.sqrt(npar1 * np.sum((DATA_SET['vz'] - DATA_SET['gasv_z'])**2)),
                'ZETA_o': np.sqrt(npar1 * np.sum((DATA_SET['vx'] - DATA_SET['gasv_x'])**2 + (DATA_SET['vy'] - DATA_SET['gasv_y'])**2 + (DATA_SET['vz'] - DATA_SET['gasv_z'])**2))}

        print('## saving calculated GLOBAL DISPERSION and DRIFT')
        io.save(SIGMA, 'sigma_'+file_no, folder=SIM.pc_datadir)
        io.save(ZETA, 'zeta_'+file_no, folder=SIM.pc_datadir)


      # calculate LOCAL DISPERSION and DRIFT
      if LOCAL:
        print('## Calculating LOCAL DISPERSION and DRIFT values')
        tmp_id = [DATA_SET[0]['idx'], DATA_SET[0]['idy'], DATA_SET[0]['idz']];
        sigma_l = np.zeros((dim.nx,dim.ny,dim.nz))
        zeta_l = np.zeros((dim.nx,dim.ny,dim.nz))
        particles_l = []
        for particle in DATA_SET:
            if [particle['idx'],particle['idy'],particle['idz']] == tmp_id:
                particles_l.append(particle)
            else:
              np_l = np.size(particles_l)
              if np_l != 0:
                  np1_l = 1./np_l
              mean_vx_l = 0; mean_vy_l = 0; mean_vz_l = 0
              sum_s_l = 0; sum_z_l = 0

              for entry in particles_l:
                  mean_vx_l = mean_vx_l + np1_l*entry['vx']
                  mean_vy_l = mean_vy_l + np1_l*entry['vy']
                  mean_vz_l = mean_vz_l + np1_l*entry['vz']

              for entry in particles_l:
                  sum_s_l = sum_s_l + (entry['vx'] - mean_vx_l)**2 + (entry['vy'] - mean_vy_l)**2 + (entry['vz'] - mean_vz_l)**2
                  sum_z_l = sum_z_l + (entry['vx'] - entry['gasv_x'])**2 + (entry['vy'] - entry['gasv_y'])**2 + (entry['vz'] - entry['gasv_z'])**2

              sigma_l[tmp_id[0]][tmp_id[1]][tmp_id[2]] = np.sqrt(np1_l * sum_s_l)
              zeta_l[tmp_id[0]][tmp_id[1]][tmp_id[2]]  = np.sqrt(np1_l * sum_z_l)

              # reset all local variables and lists to the new state (ie towards the newst particle)
              tmp_id = [particle['idx'],particle['idy'],particle['idz']]
              particles_l = []
              particles_l.append(particle)

        # do this local calculations a last time for the last grid cell
        np_l = np.size(particles_l)
        np1_l = 1./np_l
        mean_vx_l = 0; mean_vy_l = 0; mean_vz_l = 0
        sum_s_l = 0; sum_z_l = 0

        for entry in particles_l:
            mean_vx_l = mean_vx_l + np1_l*entry['vx']
            mean_vy_l = mean_vy_l + np1_l*entry['vy']
            mean_vz_l = mean_vz_l + np1_l*entry['vz']

        for entry in particles_l:
            sum_s_l = sum_s_l + (entry['vx'] - mean_vx_l)**2 + (entry['vy'] - mean_vy_l)**2 + (entry['vz'] - mean_vz_l)**2
            sum_z_l = sum_z_l + (entry['vx'] - entry['gasv_x'])**2 + (entry['vy'] - entry['gasv_y'])**2 + (entry['vz'] - entry['gasv_z'])**2

        sigma_l[tmp_id[0]][tmp_id[1]][tmp_id[2]] = np.sqrt(np1_l * sum_s_l)
        zeta_l[tmp_id[0]][tmp_id[1]][tmp_id[2]]  = np.sqrt(np1_l * sum_z_l)

        # save sigma, zeta locally and globally to SIM.pc_datadir
        print('## saving calculated LOCAL DISPERSION and DRIFT')
        io.save(sigma_l, 'sigma_l_'+file_no, folder=join(SIM.pc_datadir, 'sigma_l'))
        io.save(zeta_l, 'zeta_l_'+file_no, folder=join(SIM.pc_datadir, 'zeta_l'))

      ## Please keep this lines as a reminder on how to add columns to an record array!
      # add colums to DATA_SET fror local zeta and sigma for the individuel particle
      #DATA_SET = add_column_to_record_array(DATA_SET, 'zeta_l', zeta_l[DATA_SET['idx'],DATA_SET['idy'],DATA_SET['idz']], dtypes='float', usemask=False, asrecarray=True)
      #DATA_SET = add_column_to_record_array(DATA_SET, 'sigma_l', sigma_l[DATA_SET['idx'],DATA_SET['idy'],DATA_SET['idz']], dtypes='float', usemask=False, asrecarray=True)



    return True
