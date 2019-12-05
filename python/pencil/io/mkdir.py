#####################################################################################
##
##	class: IO
##
##	This scripts checks the existance of a folder an creates it if necessary 
##
#####################################################################################


def mkdir(destination):
  ##
  ## external modules
  import os
  
  ## 
  ## create dir if not existing
  if not os.path.exists(destination):
    os.makedirs(destination)
    return True
  return False