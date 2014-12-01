#!/usr/bin/env python

class zaverage:
  """A dummy class for explicitly adding attributes.

  Something like collections.namedtuple may be a better solution.
  """
  pass


class ZAverage():

  def __init__(self,filename='zaverages.dat',datadir='data',infile='zaver.in',
                    record_length=4,read_deltay=False):

    from re import compile
    from dim import read_dim
    from numpy import float32,float64,array

    p = compile('mxy$')
    fid = open(infile,'rt')
    varnames = [p.sub('',line.strip()) for line in fid.readlines()]
    fid.close()

    while True:
      try:
        varnames.remove('')
      except:
        break

    self.__varnames = varnames

    dim = read_dim(datadir)

    self.__dtype = {'S': lambda: float32, 'D': lambda: float64}[dim.precision]()
    self.__shape = array((dim.nx,dim.ny,len(varnames)))
    self.__count = self.__shape.prod()

    self.__record_length = record_length

    self.__file = open(datadir+'/'+filename,'rb')

    self.__read_deltay = read_deltay


  def next(self):
    from numpy import fromfile

    file = self.__file
    dtype = self.__dtype
    shape = self.__shape
    count = self.__count
    record_length = self.__record_length
    varnames = self.__varnames
    read_deltay = self.__read_deltay

    if not file.read(record_length): raise StopIteration
    self.t = fromfile(file, dtype=dtype, count=1)[0]
    file.read(record_length)

    file.read(record_length)
    data = fromfile(file, dtype=dtype, count=count).reshape(shape, order='F')
    file.read(record_length)

    if read_deltay:
      file.read(record_length)
      self.deltay = fromfile(file, dtype=dtype, count=1)[0]
      file.read(record_length)

    zaver = zaverage()
    zaver.t = self.t
    for i in range(len(varnames)):
      setattr(zaver, varnames[i], data[:,:,i])

    return zaver

  def __iter__(self):

    return self

  def __del__(self):

    self.__file.close()
