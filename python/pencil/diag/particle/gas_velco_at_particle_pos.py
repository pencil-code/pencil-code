def gas_velo_at_particle_pos(varfiles='last4', sim=False, scheme='tsc', use_IDL=False, OVERWRITE=False):
  """This script calulates the gas velocity at the particle position and stores this together
  with particle position, containing grid cell idicies, particle velocities, and particle index
  in a gas_velo_at_particle_pos file.

  Args:
    varfiles:       specifiy varfiles for calculation, e.g. 'last', 'first',
                        'all', 'VAR###', 'last4', 'first3'
    scheme:         possible are:
                        - ngp: nearest grid point
                        - cic: cloud in cell
                        - tsc: triangular shaped cloud
    OVERWRITE:		set to True to overwrite already calculated results
  """

  import os
  from pencil import io
  from pencil import read
  from pencil import get_sim
  from os.path import exists
  import numpy as np

  GAS_VELO_TAG = 'gas_velo_at_particle_pos'

  if sim == False:
      sim = get_sim()
      if sim == False:
          print('! ERROR: Specify simulation object!')
          return False
  SIM = sim

  if use_IDL:
      print('? WARNING: IDL VERSION OF THIS SCRIPT BY JOHANSEN, not recommended for 2D data')
      from pencil.backpack import pidly
      print('## starting IDL engine..')
      IDL = pidly.IDL(long_delay=0.05)		# start IDL engine

      ## skip if nothing is new
      if (not OVERWRITE) and (exists(os.path.join(SIM.pc_datadir, 'sigma.pkl'))) and (exists(os.path.join(SIM.pc_datadir, 'zeta.pkl'))):
          print('~ '+SIM.name+' is already calculated and up-to-date! -> skipping it!')

      else:
          ## start calculations
          print('~ Calculating gas_velo_at_particle_pos for "'+SIM.name+'" in "'+SIM.path+'"')
          IDL.pro('gas_velo_at_particle_pos', datadir=SIM.datadir, destination=GAS_VELO_TAG, doforthelastNvar=varfiles[4:])
          files = [i.split('_')[-1].split('.sav')[0] for i in os.listdir(os.path.join(SIM.pendatadir,GAS_VELO_TAG)) if i.startswith(GAS_VELO_TAG) and i.endswith('.sav') or i.endswith('.pkl')]
          if files == []: print('!! ERROR: No calc_gas_speed_at_particle_position-files found for '+SIM.name+'! Use idl script to produce them first!')

      IDL.close()
      return True

  else:
      print('~ Calculating gas_velo_at_particle_pos for "'+SIM.name+'" in "'+SIM.path+'"')
      save_destination = os.path.join(SIM.pc_datadir, GAS_VELO_TAG); io.mkdir(save_destination)
      varlist = SIM.get_varlist(pos=varfiles, particle=False); pvarlist = SIM.get_varlist(pos=varfiles, particle=True)

      for f, p in zip(varlist, pvarlist):
          save_filename = GAS_VELO_TAG+'_'+scheme+'_'+f[3:]
          if not OVERWRITE and exists(save_filename, folder=save_destination): continue

          print('## Reading '+f+' ...')
          ff = read.var(datadir=SIM.datadir, varfile=f, quiet=True, trimall=False)
          pp = read.pvar(datadir=SIM.datadir, varfile=p)

          ## remove ghost zones from grid, call the reduced grid the "real grid"
          realgridx = ff.x[ff.l1:ff.l2]; realgridy = ff.y[ff.m1:ff.m2]; realgridz = ff.z[ff.n1:ff.n2]
          nx = ff.l2-ff.l1; ny = ff.m2-ff.m1; nz = ff.n2-ff.n1

          ## prepare list for all quantities
          l_ipars = pp.ipars                                # particle number   KNOWN
          l_px = pp.xp;  l_py = pp.yp;  l_pz = pp.zp        # particle absolut position KNOWN
          l_vx = pp.vpx; l_vy = pp.vpy; l_vz = pp.vpz       # particle velocity KNOWN
          l_rix = [];    l_riy = [];    l_riz = []          # particle untrimmed realgrid index (grid index = l/m/n + readgrid index ???)
          l_ix = [];     l_iy = [];     l_iz = []           # particle grid index (in untrimmed grid)
          l_ux = [];     l_uy = [];     l_uz = []           # underlying gas velocity at position of particle

          ## get index of realgrid cell for each particle
          for i in range(len(l_ipars)):
              l_rix.append(np.abs(realgridx-l_px[i]).argmin())
              l_riy.append(np.abs(realgridy-l_py[i]).argmin())
              l_riz.append(np.abs(realgridz-l_pz[i]).argmin())

          ## convert into untrimmed grid
          l_ix = np.array(l_rix)+ff.l1; l_iy = np.array(l_riy)+ff.m1; l_iz = np.array(l_riz)+ff.n1

          ## NGP
          if scheme == 'ngp' or scheme == 'NGP':
              print('## Calculating gas velocities via '+scheme)
              l_ux = ff.ux[l_iz, l_iy, l_ix]
              l_uy = ff.uy[l_iz, l_iy, l_ix]
              l_uz = ff.uz[l_iz, l_iy, l_ix]

          ## CIC
          if scheme == 'cic' or scheme == 'CIC':
              print('## Calculating gas velocities via '+scheme)
              for ix0, iy0, iz0, px, py, pz in zip(l_ix, l_iy, l_iz, l_px, l_py, l_pz):     # for each particle
                  if ff.x[ix0] > px: ix0 = ix0-1            # ix0 must be left to particle
                  if ff.y[iy0] > py: iy0 = iy0-1            # iy0 must be below the particle
                  if ff.z[iz0] > pz: iz0 = iz0-1            # iz0 must be under particle

                  ix1=ix0; iy1=iy0; iz1=iz0                 # if a dim. is zero, this is default, else:
                  if nx > 1: ix1 = ix0+1; dx_1 = 1./ff.dx   # if a dim is non-zero, ajust ix1 to right cell
                  if ny > 1: iy1 = iy0+1; dy_1 = 1./ff.dy   # if a dim is non-zero, ajust iy1 to above cell
                  if nz > 1: iz1 = iz0+1; dz_1 = 1./ff.dz   # if a dim is non-zero, ajust iz1 to above cell

                  ux = 0.; uy = 0.; uz = 0.
                  for ix in [ix0, ix1]:
                      for iy in [iy0, iy1]:
                          for iz in [iz0, iz1]:
                              weight = 1.
                              if nx > 1: weight = weight * ( 1. - abs(px-ff.x[ix])*dx_1)
                              if ny > 1: weight = weight * ( 1. - abs(py-ff.y[iy])*dy_1)
                              if nz > 1: weight = weight * ( 1. - abs(pz-ff.z[iz])*dz_1)

                              ux = ux + weight*ff.ux[iz, iy, ix]
                              uy = uy + weight*ff.uy[iz, iy, ix]
                              uz = uz + weight*ff.uz[iz, iy, ix]

                              if iz0 == iz1: break      # beware of degeneracy:
                          if iy0 == iy1: break      # beware of degeneracy:
                      if ix0 == ix1: break      # beware of degeneracy:

                  l_ux.append(ux); l_uy.append(uy); l_uz.append(uz)



          ## TSC
          if scheme == 'tsc' or scheme == 'TSC':
              for ix0, iy0, iz0, px, py, pz in zip(l_ix, l_iy, l_iz, l_px, l_py, l_pz):  # for each particle
                ixx0 = ix0; ixx1 = ix0      # beware of degeneracy
                iyy0 = iy0; iyy1 = iy0
                izz0 = iz0; izz1 = iz0

                if nx > 1: ixx0 = ix0-1; ixx1 = ix0+1; dx_1 = 1./ff.dx; dx_2 = 1./ff.dx**2
                if ny > 1: iyy0 = iy0-1; iyy1 = iy0+1; dy_1 = 1./ff.dy; dy_2 = 1./ff.dy**2
                if nz > 1: izz0 = iz0-1; izz1 = iz0+1; dz_1 = 1./ff.dz; dz_2 = 1./ff.dz**2

                ux = 0.; uy = 0.; uz = 0.
                for ix in [ix0, ixx0, ixx1]:
                    weight_x = 0.
                    if ix-ix0 == -1 or ix-ix0 == 1:
                        weight_x = 1.125 - 1.5*abs(px-ff.x[ix])*dx_1 + 0.5*abs(px-ff.x[ix])**2*dx_2
                    elif nx != 1:
                        weight_x = 0.75  - (px-ff.x[ix])**2*dx_2

                    for iy in [iy0, iyy0, iyy1]:
                        weight_y = 0.
                        if iy-iy0 == -1 or iy-iy0 == 1:
                            weight_y = 1.125 - 1.5*abs(py-ff.y[iy])*dy_1 + 0.5*abs(py-ff.y[iy])**2*dy_2
                        elif ny != 1:
                            weight_y = 0.75  - (py-ff.y[iy])**2*dy_2

                        for iz in [iz0, izz0, izz1]:
                            weight_z = 0.
                            if iz-iz0 == -1 or iz-iz0 == 1:
                                weight_z = 1.125 - 1.5*abs(pz-ff.z[iz])*dz_1 + 0.5*abs(pz-ff.z[iz])**2*dz_2
                            elif nz != 1:
                              weight_z = 0.75  - (pz-ff.z[iz])**2*dz_2

                            weight = 1.
                            if nx > 1: weight = weight * weight_x
                            if ny > 1: weight = weight * weight_y
                            if nz > 1: weight = weight * weight_z

                            ux = ux + weight*ff.ux[iz, iy, ix]
                            uy = uy + weight*ff.uy[iz, iy, ix]
                            uz = uz + weight*ff.uz[iz, iy, ix]

                            if izz0 == izz1: break      # beware of degeneracy:
                        if iyy0 == iyy1: break      # beware of degeneracy:
                    if ixx0 == ixx1: break      # beware of degeneracy:

                l_ux.append(ux); l_uy.append(uy); l_uz.append(uz)


          ## Convert all information into a single record array
          data_set = np.core.records.fromarrays(
                        [l_ipars.astype('int'),
                        l_px, l_py, l_pz,
                        l_vx, l_vy, l_vz,
                        l_rix, l_riy, l_riz,
                        l_ix, l_iy, l_iz,
                        l_ux, l_uy, l_uz],
                        names = 'ipar, ipx, ipy, ipz, vx, vy, vz, rix, riy, riz, ix, iy, iz, ux, uy, uz',
                        formats = 'int, float, float, float, float, float, float, int, int, int, int, int, int, float, float, float'
                        )
          gas_velo_at_particle_pos = np.sort(data_set, order=['ix', 'iy', 'iz'])

          Nix = int(gas_velo_at_particle_pos['rix'].max()+1)
          Niy = int(gas_velo_at_particle_pos['riy'].max()+1)
          Niz = int(gas_velo_at_particle_pos['riz'].max()+1)

          Npar_arr = np.array([gas_velo_at_particle_pos['rix'], gas_velo_at_particle_pos['riy'], gas_velo_at_particle_pos['riz']])
          #rgrid_edges = (grid.x[1:]-(grid.x[1:]-grid.x[:-1])/2)[2:-2]
          xrange=np.arange(0,float(gas_velo_at_particle_pos['rix'].max())+2); xrange=xrange-0.5
          yrange=np.arange(0,float(gas_velo_at_particle_pos['riy'].max())+2)
          zrange=np.arange(0,float(gas_velo_at_particle_pos['riz'].max())+2)

          Npar_hist, edges = np.histogramdd(Npar_arr.T, bins=(xrange, yrange, zrange))
          Npar_hist, edges = np.histogramdd(Npar_arr.T, bins=(Nix, Niy, Niz))


          gas_velo_at_particle_pos = {
                    'time': ff.t,
                    'par_pos': np.array([gas_velo_at_particle_pos['ipx'],
                                         gas_velo_at_particle_pos['ipy'],
                                         gas_velo_at_particle_pos['ipz']]),
                    'par_velo': np.array([gas_velo_at_particle_pos['vx'],
                                          gas_velo_at_particle_pos['vy'],
                                          gas_velo_at_particle_pos['vz']]),
                    'par_idx': np.array([gas_velo_at_particle_pos['rix'],
                                         gas_velo_at_particle_pos['riy'],
                                         gas_velo_at_particle_pos['riz']]),
                    'npar': np.array(Npar_hist[gas_velo_at_particle_pos['rix'],
                                               gas_velo_at_particle_pos['riy'],
                                               gas_velo_at_particle_pos['riz']]),
                    'gas_velo': np.array([gas_velo_at_particle_pos['ux'],
                                         gas_velo_at_particle_pos['uy'],
                                         gas_velo_at_particle_pos['uz']])
                    }


          print('## Saving dataset into '+save_destination+'...')
          io.pkl_save({'gas_velo_at_particle_pos': gas_velo_at_particle_pos, 't': ff.t}, save_filename, folder=save_destination)
      print('## Done!')
