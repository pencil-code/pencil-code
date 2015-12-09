#################################
#
# The __init__ file is used not only to import the sub-modules, but also to set everything up properly. 
#
#################################

############
# read external modules
import os

############
# read internal sub-modules
from .. import io		# input und output functions, like save data or call IDL scripts
from .. import diag		# diagnostic scripts and functions
from .. import visu		# visualisation routines
from .. import calc		# math functions and further calculations
from .. import vtk		# 
from .. import math
from .. import read		# read data and parameters from pencil code directory
from .. import tool_kit		# all nice workarounds get stored here (e.g., resubmit script)
from .. import export		# exporter (e.g., vtk, xml)


# do a check on simulation data dir and read your simulation setup in simuDict
if os.path.isdir("./data"):
  print("## Pencil Code Simulation found!")
  # -> produce simulation setup dictionary here 
else:
  print(("?? WARNING: No './data' directory found! Current dir is '"+os.getcwd()+"'. Please give your data directory manually to all pen functions."))
 
  
# -> load existing or create new diagnostics defaults