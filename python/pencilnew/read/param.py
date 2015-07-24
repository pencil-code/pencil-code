# $Id$

def param(datadir='data/', param2=False, quiet=False):
  """
  param.py -- reads parameters from the f90 namelists
  Author: J. Oishi (joishi@amnh.org)

  REQUIRES: nl2python perl script (based on Wolfgang Dobler's nl2idl script)

  todo: think about single/double precision; make things into numpy arrays?
  """
  import numpy as N
  import os
  
  datadir = os.path.expanduser(datadir)

  if (param2):
    filen = os.path.join(datadir,'param2.nml')	# read &INIT_PARS
  else:
    filen = os.path.join(datadir,'param.nml')	# read &RUN_PARS

  # will need to read dim.dat to deal with precision, should that be necessary
  #dim = read_dim(datadir)

  # execute output of nl2python script
  if (not os.path.exists(filen)):
      print "read_param: no such file",filen
      raise ValueError
  
  cmd = 'nl2python '+filen
  script = os.popen(cmd).read()
  if (not quiet): print script;
  if (len(script) != 0):
      exec(script)
  else:
      print "read_param: nl2python returned nothing! is $PENCIL_HOME/bin in the path?"
      return -1
  
  ret = Params() # par is the name of the class

  return ret

