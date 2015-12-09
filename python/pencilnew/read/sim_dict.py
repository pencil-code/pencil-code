# $Id$
#
# read simulation dictionary
#
# Author: A. Schreiber (aschreiber@mpia.de)
#
#

import pen.io.load_obj as load_obj


def sim_dict(*args, **kwargs):
  """
  Read simulation setup
  """
  
  return sim_dict(*args, **kwargs)

class sim_dict(object):
  
  def __init__(self, datadir='./data'):
    """
    Read simulation data and store it in a dictionaly like object.
    """
    print '## Loading simulation parameters..'
    simuDict = io.load_obj(simuDict_name)
    if not simuDict:
      simuDict = {}
      grid = pen.read.grid(quiet=True,trim=True)
      simuDict['L_x'] = grid.Lx
      simuDict['L_y'] = grid.Ly
      simuDict['L_z'] = grid.Lz
      simuDict['d_x'] = grid.dx
      simuDict['d_y'] = grid.dy
      simuDict['d_z'] = grid.dz
      simuDict['grid_x'] = grid.x
      simuDict['grid_y'] = grid.y
      simuDict['grid_z'] = grid.z
      if len(grid.x)==1:
	simuDict['2d']='yz'
      elif len(grid.y)==1:
	simuDict['2d']='xz'
      elif len(grid.z)==1:
	simuDict['2d']='xy'
      else:
	simuDict['2d']=False
      io.save_obj(simuDict, simuDict_name)