
from .. import get_sims
from ..math import is_float, is_int, is_number
import sys
import os
from os.path import join
import numpy as np

print('### PENCIL - FIX DIFFUSIVITIES STARTED ###')

######################
# defaults

## typical diffusivity values:

### for lower vaule use:
typ_L = 0.001           # simulation domain size
typ_N = 128             # number of grid cells
typ_diff = 1.7e-29      # corrosponding diffusivity value

### from Anders Johansen Sedimentation Sample:
#typ_L = 0.2           # simulation domain size
#typ_N = 128             # number of grid cells
#typ_diff = 4e-18      # corrosponding diffusivity value

### for higher value use:
#typ_L = 0.1           # simulation domain size
#typ_N = 125             # number of grid cells
#typ_diff = 1.7e-18      # corrosponding diffusivity value

## sets of diffusivity quantities
diff_set_runin = ['diffrho_hyper3', 'nu_hyper3']
print('## Fixing the following values in run.in: '+str(diff_set_runin))

########################
## main script
##
## doing:  nu_new = nu_old * (dx_new / dx_old)**5, maximum in x, y, z
SIMs = get_sims()

if SIMs == False: sys.exit('!! ERROR: No simulation found! No SIMDICT found!')

for SIM in SIMs:
  print('\n## Doing fixing for Simulation: '+SIM.name)
  ## paths
  p_cparam = join(SIM.path, 'src/cparam.local')
  p_runin = join(SIM.path, 'run.in')
  p_startin = join(SIM.path, 'start.in')

  ## get data
  Nx = SIM.get_value('nxgrid')
  Ny = SIM.get_value('nygrid')
  Nz = SIM.get_value('nzgrid')

  num_arr = SIM.get_value('Lxyz')
  Lx = num_arr[0]
  Ly = num_arr[1]
  Lz = num_arr[2]

  dx_sim = np.min([Lx/Nx,Ly/Ny,Lz/Nz])
  dx_old = typ_L/typ_N

  new_diff = typ_diff * (dx_sim / dx_old)**5
  print('## Old diffusivity value is: '+str(SIM.get_value(diff_set_runin[0])))
  print('## New diffusivity value is: '+str(new_diff))

  ## change value in file
  for q in diff_set_runin:
    SIM.change_value_in_file(p_runin, q, newValue=new_diff)


print('## DONE!\n')
