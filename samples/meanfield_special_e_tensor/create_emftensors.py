import argparse
import gc
import os
import re

import numpy as np
from scipy.special import erf
import h5py

import pencil as pc

def labelDataset(dataset):

    labels  = ['t','z','y','x']

    coefdims    = len(dataset.shape) - len(labels)

    for coefdim in np.arange(0,coefdims):
        labels = ['coef%d' % coefdim] + labels

    for dim,label in zip(np.arange(0,len(labels)), labels):
        dataset.dims[dim].label = label

def printH5Structure(h5file):
    print 'File has following structure:' 
    datanames = [] 
    h5file.visit(datanames.append)
    for dataname in datanames:
        print dataname 

def getDatasetNames(h5file):
    datanames = []
    def addDataname(dataname):
        if isinstance(h5file[dataname],h5py.Dataset):
            datanames.append('/{0}'.format(dataname))
    h5file.visit(addDataname)
    return datanames

class FieldCreator:

  def __init__(self, pencilfolder):

    self.datadir    = os.path.join(pencilfolder, 'data')
    try:
        self.grid       = pc.read_grid(self.datadir, trim=True, quiet=True)
    except Exception,e:
        print "Cannot read grid. Have you run 'pc_run start' already?"
    self.datatype   = self.grid.x.dtype
    self.h5file     = os.path.join(self.datadir, 'emftensors.h5')

  def createField(self, **fieldargs):
    
    print 'Adding dataset emftensor/{field}/{dataset} to the tensor file'.format(**fieldargs)
    extradims_dict =  { 
                      'alpha' : [3,3],
                      'beta'  : [3,3],
                      'gamma' : [3],
                      'delta' : [3],
                      'kappa' : [3,3,3],
                      'utensor' : [3],
                      'omega' : [3]
                      }

    field = fieldargs['field']
    dataset   = fieldargs['dataset']
    fullname  = '/emftensor/{field}/{dataset}'.format(**fieldargs)
    with h5py.File(self.h5file,'a') as f:
      
      datasets          = getDatasetNames(f)
      extradims = extradims_dict.get(field,[])

      datashape   = [fieldargs['tlen'], len(self.grid.z), len(self.grid.y), len(self.grid.x)]
      fieldshape  = extradims + datashape
      maxshape    = [fieldargs['tlen'], len(self.grid.z), len(self.grid.y), len(self.grid.x)]
      maxshape    = extradims + maxshape
      chunks      = tuple([ 1 for x in extradims ] + [1, len(self.grid.z), len(self.grid.y), len(self.grid.x)])

      data        = np.zeros(fieldshape, dtype=self.datatype)

      t, z, y, x = np.meshgrid(np.arange(fieldargs['tlen']), self.grid.z, self.grid.y, self.grid.x)
  
      x = np.reshape(x, datashape)
      y = np.reshape(y, datashape)
      z = np.reshape(z, datashape)
      t = np.reshape(t, datashape)

      values = np.asarray(fieldargs['values'], dtype='O')
      for i in xrange(len(values.shape)-1):
        values = np.transpose(values)
      if values.shape == tuple(extradims) and len(extradims) > 0:
        for index in np.ndindex(values.shape):
          #print index, values[index]
          if callable(values[index]):
            data[index] = values[index](t,z,y,x)
          else:
            data[index] = values[index]
      elif len(values) == 1:
        pass
      else:
        raise Exception('Got {0} values but do not know how to fill {1} indices.'.format(values.shape, extradims))
  
      if 'emftensor' not in f:
        emftensor_grp = f.create_group('emftensor')
      else:
        emftensor_grp = f['/emftensor']
      if field not in emftensor_grp:
        field_grp = emftensor_grp.create_group(field)
      else:
        field_grp = emftensor_grp[field]
      
      if dataset in field_grp:
        del field_grp[dataset]
      ds    = field_grp.create_dataset(dataset, fieldshape,  \
                                              maxshape=maxshape, \
                                              chunks=chunks, \
                                              dtype=self.datatype,
                                              data=data)
      labelDataset(ds)

emffields = [ 
              {
              'field'     : 'alpha',
              'dataset'   : 'diag',
              'values'    : [ \
                              [ 1, 0, 0 ], \
                              [ 0, 1, 0 ], \
                              [ 0, 0, 1 ]
                            ],
              'tlen'      : 1
              },
              {
              'field'     : 'alpha',
              'dataset'   : 'uniform',
              'values'    : [ \
                              [ 1, 1, 1 ], \
                              [ 1, 1, 1 ], \
                              [ 1, 1, 1 ]
                            ],
              'tlen'      : 1
              },
              {
              'field'     : 'alpha',
              'dataset'   : 'ns-asymmetric-cos',
              'values'    : [ \
                              [ lambda t,z,y,x: -np.cos(y), 0, 0 ], \
                              [ 0, lambda t,z,y,x: -np.cos(y), 0 ], \
                              [ 0, 0, lambda t,z,y,x: -np.cos(y) ], \
                            ],
              'tlen'      : 1
              },
              {
              'field'     : 'alpha',
              'dataset'   : 'ns-asymmetric-erf',
              'values'    : [ \
                              [ lambda t,z,y,x: erf((y-np.pi/2)/0.05), 0, 0 ], \
                              [ 0, lambda t,z,y,x: erf((y-np.pi/2)/0.05), 0 ], \
                              [ 0, 0, lambda t,z,y,x: erf((y-np.pi/2)/0.05) ], \
                            ],
              'tlen'      : 1
              },
              {
              'field'     : 'alpha',
              'dataset'   : 'Steenbeck-Krause-1969-model1',
              'values'    : [ \
                              [ 0, 0, 0 ], \
                              [ 0, 0, 0 ], \
                              [ 0, 0, lambda t,z,y,x: 0.5*(1+erf((x-0.9)/0.075))*np.cos(y) ], \
                            ],
              'tlen'      : 1
              },
              {
              'field'     : 'alpha',
              'dataset'   : 'Jouve-2008-benchmark',
              'values'    : [ \
                              [ lambda t,z,y,x: (3.0*np.sqrt(3.0)/4.0)*np.power(np.sin(y),2.0)*np.cos(y)*(1.0+erf((x-0.7)/0.02)), 0, 0 ], \
                              [ 0, lambda t,z,y,x: (3.0*np.sqrt(3.0)/4.0)*np.power(np.sin(y),2.0)*np.cos(y)*(1.0+erf((x-0.7)/0.02)), 0 ], \
                              [ 0, 0, lambda t,z,y,x: (3.0*np.sqrt(3.0)/4.0)*np.power(np.sin(y),2.0)*np.cos(y)*(1.0+erf((x-0.7)/0.02)) ], \
                            ],
              'tlen'      : 1
              },
              {
              'field'     : 'beta',
              'dataset'   : 'diag',
              'values'    : [ \
                              [ 1, 0, 0 ], \
                              [ 0, 1, 0 ], \
                              [ 0, 0, 1 ]
                            ],
              'tlen'      : 1
              },
              {
              'field'     : 'beta',
              'dataset'   : 'Jouve-2008-benchmark',
              'values'    : [ \
                              [ lambda t,z,y,x: 0.01 + 0.5*(1-0.01)*(1.0+erf((x-0.7)/0.02)), 0, 0 ], \
                              [ 0, lambda t,z,y,x: 0.01 + 0.5*(1-0.01)*(1.0+erf((x-0.7)/0.02)), 0 ], \
                              [ 0, 0, lambda t,z,y,x: 0.01 + 0.5*(1-0.01)*(1.0+erf((x-0.7)/0.02)) ], \
                            ],
              'tlen'      : 1
              },
              {
              'field'     : 'gamma',
              'dataset'   : 'uniform',
              'values'    : [ 1, 1, 1 ],
              'tlen'      : 1
              },
              {
              'field'     : 'utensor',
              'dataset'   : 'uniform',
              'values'    : [ 1, 1, 1 ],
              'tlen'      : 1
              },
              {
              'field'     : 'utensor',
              'dataset'   : 'constant-omega',
              'values'    : [ 
                              0, \
                              0, \
                              lambda t,z,y,x: x*np.sin(y)  \
                            ],
              'tlen'      : 1
              },
              {
              'field'     : 'utensor',
              'dataset'   : 'simple-diffrot',
              'values'    : [ 
                              0, \
                              0, \
                              lambda t,z,y,x: -x*np.sin(y)*(x-x.min())/(x.max()-x.min()) 
                            ],
              'tlen'      : 1
              },
              {
              'field'     : 'utensor',
              'dataset'   : 'Steenbeck-Krause-1969-model1',
              'values'    : [ 
                              0, \
                              0, \
                              lambda t,z,y,x: 0.5*(1-erf((x-0.7)/0.075))*x*np.sin(y)  \
                            ],
              'tlen'      : 1
              },
              {
              'field'     : 'utensor',
              'dataset'   : 'Joure-2008-benchmark',
              'values'    : [
                              0, \
                              0, \
                              lambda t,z,y,x: x*np.sin(y)*(0.92 + 0.5*(1.0 + erf((x-0.7)/0.02))*(1.0-0.92-0.2*np.power(np.cos(y),2.0))) \
                            ],
              'tlen'      : 1
              },
              {
              'field'     : 'utensor',
              'dataset'   : 'Jouve-2008-benchmark-noavg',
              'values'    : [
                              0, \
                              0, \
                              lambda t,z,y,x: x*np.sin(y)*(-0.011125 + 0.5*(1.0 + erf((x-0.7)/0.02))*(1.0-0.92-0.2*np.power(np.cos(y),2.0))) \
                            ],
              'tlen'      : 1
              }
             ]

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Create artificial emf mean fields')
    parser.add_argument('pencilfolder', nargs='1', help='Folder to analyze')
    #parser.add_argument('folder',nargs=1, help='Pencil-code folder to process')
    args = parser.parse_args()
    pencilfolder    = args.pencilfolder[0]
    fieldcreator    = FieldCreator(pencilfolder)
    for fieldargs in emffields:
      fieldcreator.createField(**fieldargs)
    with h5py.File(fieldcreator.h5file) as f:
        printH5Structure(f)
        #print 'File has the following datasets:'
        #for ds in datasets:
        #    print ds
