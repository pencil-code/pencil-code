
def dill_exists(name, folder=False):
  """This scripts checks if a certain dill-file already exists.

  Args:
	name:		Name of dill file  (<name>.dill)
	folder:		Folder containing dill file
  """

  if (not name.endswith('.dill')):	name = name+'.dill'

  # if folder is not defined try to find file at typical places
  if not folder:
      if __exists__(__join__('.pc', name)):
          folder = '.pc'
      elif __exists__(__join__('data/.pc', name)):
          folder = 'data/.pc'
      else:
          return False

  file = __join__(folder, name)
  try:							# check on existance
    if not __exists__(file):
      return False
    return True

  except:						# if anything goes wrong
	print '!! ERROR: Something went wrong when checking for the dill file!'
	return False
