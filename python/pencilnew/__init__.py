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
import io		# input und output functions, like save data or call IDL scripts
import diag		# diagnostic scripts and functions
import visu		# visualisation routines
import calc		# math functions and further calculations
import vtk		# 
import math
import read		# read data and parameters from pencil code directory
import tool_kit		# all nice workarounds get stored here (e.g., resubmit script)
import export		# exporter (e.g., vtk, xml)


# do a check on simulation data dir and read your simulation setup in simuDict
if os.path.isdir("./data"):
  print "## Pencil Code Simulation found!"
  # -> produce simulation setup dictionary here 
else:
  print "?? WARNING: No './data' directory found! Current dir is '"+os.getcwd()+"'. Please give your data directory manually to all pen functions."
 
  
# -> load existing or create new diagnostics defaults