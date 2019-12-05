
def walklevel(some_dir, level=1):
  """Walklevel walks through higher level dirs to a certain depth.
  Is used as an interator in loops, e.g.:
	  sim_paths = []
	  for path,dirs,files in walklevel('.', depth):
		if 'start.in' in files and 'run.in' in files:
		  print '## Found Simulation in '+path
		  sim_paths.append(path)

  Args:
	some_dir:		at which dir to start the walklevel algorithm
	level:		how deep this algorithm shall walk through dirs, default is 1

  Yields:
	yield root, dirs, files
  """
  import os

  some_dir = some_dir.rstrip(os.path.sep)
  assert os.path.isdir(some_dir)
  num_sep = some_dir.count(os.path.sep)
  for root, dirs, files in os.walk(some_dir):
	  yield root, dirs
	  num_sep_this = root.count(os.path.sep)
	  if num_sep + level <= num_sep_this:
		  del dirs[:]
