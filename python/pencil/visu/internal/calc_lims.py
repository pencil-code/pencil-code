def calc_lims(val_min=False, val_max=False, offset=False):
  """
  Gives back an optimal limit for an plot upper or lower value. Therefor, it uses an inbuild offset from the PROFILE.
  """
  from pen.profiles import get
  import numpy as np
  import math
  if not offset: offset = get()['ax_xylim_moffset']


  if not isinstance(val_min, (bool)):
    if val_min == 0.: return 0.

    val_min = float(val_min - offset * np.abs(val_min))

    e = int(math.log10(abs(val_min)))
    if e < 0: e = e-1
    M = val_min/(10**e)
    val_min = np.floor(M*10.)/10.*10**e

    return val_min


  if not isinstance(val_max, (bool)):
    if val_max == 0.: return 0.

    val_max = float(val_max*1.0 + offset * np.abs(val_max))

    e = int(math.log10(abs(val_max)))
    if e < 0: e = e-1
    M = val_max/(10**e)
    val_max = np.ceil(M*10.)/10.*10**e

    return val_max

  return False
