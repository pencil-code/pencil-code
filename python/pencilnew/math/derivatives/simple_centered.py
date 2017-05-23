
def simple_centered(x, y):
  """
  calculate dy by center differencing using array slices
  """
  import numpy as np

  dy = np.zeros(y.shape,np.float) #we know it will be this size
  dy[1:-1] = (y[2:] - y[0:-2])/(x[2:] - x[0:-2])
  dy[0] = (y[1]-y[0])/(x[1]-x[0])
  dy[-1] = (y[-1] - y[-2])/(x[-1] - x[-2])

  return dy
