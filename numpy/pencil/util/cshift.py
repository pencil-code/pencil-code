# $Id: cshift.py,v 1.1 2008-02-05 15:29:32 theine Exp $
#
# Implemented fortran style cshift.
# (Apparently this is not provided by numpy.)
#
# Author: T. Heinemann (T.Heinemann@damtp.cam.ac.uk)
#

__version__ = "$Id: cshift.py,v 1.1 2008-02-05 15:29:32 theine Exp $"

import numpy

def cshift(f,shift,axis=0):

  """Provides Fortran 90 style cshift.

  f     -- Input array of arbitrary dimension.
  shift -- Shift by this much.
  axis  -- The dimension along which to shift.
           (optional; defaults to first dimension)

  """

  if axis >= f.ndim:
    print 'Error: axis keyword must not exceed number of dimensions.'
    return False

  colons = [':' for i in range(f.ndim)]

  a = colons; a[axis] = '%d:' % shift; a = ','.join(a)
  b = colons; b[axis] = ':%d' % shift; b = ','.join(b)

  exec('g = numpy.concatenate([f['+a+'],f['+b+']],axis=%d)') % axis

  return g
