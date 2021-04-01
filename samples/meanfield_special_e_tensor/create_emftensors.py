from __future__ import print_function
import argparse
import os

import numpy as np
from scipy.special import erf
import h5py

import pencil_old as pc


def labelDataset(dataset):
    """Labels all dimensions of a dataset.

    Parameters
    ----------
    dataset : h5py.Dataset
        Dataset to label.
    """

    labels = ['z', 'y', 'x', 't']

    coefdims = len(dataset.shape) - len(labels)

    for coefdim in np.arange(0, coefdims):
        labels = ['coef%d' % coefdim] + labels

    for dim, label in zip(np.arange(0, len(labels)), labels):
        dataset.dims[dim].label = label


def printH5Structure(h5file):
    """Prints all groups and datasets in a file.

    Parameters
    ----------
    h5file : h5py.File
        File that will be printed.
    """

    print('\n' + 80*'-' + '\n')
    print('File has the following structure:\n')
    datanames = []
    h5file.visit(datanames.append)
    for dataname in datanames:
        ngroups = dataname.count('/')
        endstr = dataname.split('/')[ngroups]
        spaces = len(dataname) - len(endstr) - ngroups
        print('  ' + spaces*' ' + endstr)
    print('\n' + 80*'-')

class CoefficientCreator:
    """CoefficientCreator creates mean-field coeffients.
    """

    def __init__(self, pencilfolder):
        """Initializes CoefficientCreator based on a Pencil-code folder.

        Parameters
        ----------
        pencilfolder : str
            Pencil-code that provides the grid and will be used
            to store the coefficients. pc_start needs to be run
            on the folder beforehand.
        """

        self.datadir = os.path.join(pencilfolder, 'data')
        try:
            self.grid = pc.read_grid(self.datadir, trim=True, quiet=True)
        except Exception as e:
            print("Cannot read grid. Have you run 'pc_run start' already?")
            raise e
        self.datatype = self.grid.x.dtype
        self.h5file = os.path.join(self.datadir, 'emftensors.h5')

    def createCoefficient(self,
                          coefficient='alpha',
                          dataset='isotropic',
                          values=np.diag([1.0, 1.0, 1.0]),
                          tvals=np.array([1.0])):
        """Creates mean-field coefficient based on given configuration parameters.

        Parameters
        ----------
        coefficient : str
            Name of the coeffient: alpha, beta, gamma, delta,
            kappa or utensor.
        dataset : str
            Name of the dataset to be stored in the coefficient
            group.
        values : array-like, float or callable
            A multi-dimensional array-like structure with floating point
            numbers or callable functions. Entries in array should be 
            either floating point numbers or functions that - when given
            z, y, x and t coordinate meshes - will produce a mesh with
            the values of the dataset at each point. Each index of the
            values array will then correspond with the dimension of the
            coefficient. E.g. alpha requires 3x3 array-like while
            kappa requires 3x3x3 array-like.
        tvals : array-like, float
            Time values to evaluate. For time-dependent coefficients.
        """

        print('Creating dataset emftensor/{coefficient}/{dataset}'.format(**coef_conf))

        # Define the number of dimensions in coefficients
        extradims_dict = {
                         'alpha': [3, 3],
                         'beta': [3, 3],
                         'gamma': [3],
                         'delta': [3],
                         'kappa': [3, 3, 3],
                         'utensor': [3],
                         }

        # Open data/emftensors.h5 for writing
        with h5py.File(self.h5file, 'a') as f:

            extradims = extradims_dict.get(coefficient, [])

            # Define the dataset shape and chunk shape
            datashape = [len(self.grid.z), len(self.grid.y), len(self.grid.x), len(tvals)]
            coefficientshape = extradims + datashape
            maxshape = [len(self.grid.z), len(self.grid.y), len(self.grid.x), len(tvals)]
            maxshape = extradims + maxshape
            chunks = tuple(
                    [1 for x in extradims] +
                    [len(self.grid.z), len(self.grid.y), len(self.grid.x), 1])

            data = np.zeros(coefficientshape, dtype=self.datatype)

            # Create meshgrid for evaluating values
            t, z, y, x = np.meshgrid(tvals, self.grid.z, self.grid.y, self.grid.x)

            x = np.moveaxis(x, 1, -1)
            y = np.moveaxis(y, 1, -1)
            z = np.moveaxis(z, 1, -1)
            t = np.moveaxis(t, 1, -1)

            # For index in values, initialize the dataset
            values = np.asarray(values, dtype='O')
            for i in range(len(values.shape)-1):
                values = np.transpose(values)
            if values.shape == tuple(extradims) and len(extradims) > 0:
                for index in np.ndindex(values.shape):
                    # If value is a function, evaluate it
                    if callable(values[index]):
                        data[index] = values[index](z, y, x, t)
                    else:
                        data[index] = values[index]
            elif len(values) == 1:
                pass
            else:
                raise Exception(
                        ('Got {0} values '
                         'but do not know how to '
                         'fill {1} indices.').format(values.shape, extradims))

            # Create groups for dataset
            if 'emftensor' not in f:
                emftensor_grp = f.create_group('emftensor')
            else:
                emftensor_grp = f['/emftensor']
            if coefficient not in emftensor_grp:
                coefficient_grp = emftensor_grp.create_group(coefficient)
            else:
                coefficient_grp = emftensor_grp[coefficient]

            if dataset in coefficient_grp:
                del coefficient_grp[dataset]

            # Save dataset
            ds = coefficient_grp.create_dataset(
                    dataset,
                    coefficientshape,
                    maxshape=maxshape,
                    chunks=chunks,
                    dtype=self.datatype,
                    data=data)

            # Label dataset
            labelDataset(ds)

"""
This list contains example coefficient configurations

alpha/isotropic - isotropic alpha
alpha/Steenbeck-Krause-1969-model1 - alpha used in model 1 in paper by Steenbeck and Krause
alpha/Jouve-2008-benchmark - alpha used in benchmarks A and B in paper by Jouve et al.
beta/isotropic - isotropic beta
beta/Jouve-2008-benchmark - beta used in benchmark B in paper by Jouve et al.
utensor/Steenbeck-Krause-1969-model1 - velocity field used in model 1 in paper by Steenbeck and Krause
utensor/Jouve-2008-benchmark - velocity field used in benchmarks A and B in paper by Jouve et al.
"""

coef_configs = [
    {
        'coefficient': 'alpha',
        'dataset': 'isotropic',
        'values': [
                  [1, 0, 0],
                  [0, 1, 0],
                  [0, 0, 1]
                  ],
        'tvals': [1]
    },
    {
        'coefficient': 'alpha',
        'dataset': 'Steenbeck-Krause-1969-model1',
        'values': [
                  [0, 0, 0],
                  [0, 0, 0],
                  [0, 0, lambda z, y, x, t: 0.5*(1+erf((x-0.9)/0.075))*np.cos(y)],
                  ],
        'tvals': [1]
    },
    {
        'coefficient': 'alpha',
        'dataset': 'Jouve-2008-benchmark',
        'values': [
                  [0, 0, 0],
                  [0, 0, 0],
                  [0, 0, lambda z, y, x, t: (3.0*np.sqrt(3.0)/4.0)*np.power(np.sin(y), 2.0)*np.cos(y)*(1.0+erf((x-0.7)/0.02))],
                  ],
        'tvals': [1]
    },
    {
        'coefficient': 'beta',
        'dataset': 'isotropic',
        'values': [
                  [1, 0, 0],
                  [0, 1, 0],
                  [0, 0, 1]
                  ],
        'tvals': [1]
    },
    {
        'coefficient': 'beta',
        'dataset': 'Jouve-2008-benchmark',
        'values': [
                  [lambda z, y, x, t: 0.01 + 0.5*(1-0.01)*(1.0+erf((x-0.7)/0.02)), 0, 0],
                  [0, lambda z, y, x, t: 0.01 + 0.5*(1-0.01)*(1.0+erf((x-0.7)/0.02)), 0],
                  [0, 0, lambda z, y, x, t: 0.01 + 0.5*(1-0.01)*(1.0+erf((x-0.7)/0.02))],
                  ],
        'tvals': [1]
    },
    {
        'coefficient': 'utensor',
        'dataset': 'Steenbeck-Krause-1969-model1',
        'values': [
                  0,
                  0,
                  lambda z, y, x, t: 0.5*(1-erf((x-0.7)/0.075))*x*np.sin(y)
                  ],
        'tvals': [1]
    },
    {
        'coefficient': 'utensor',
        'dataset': 'Jouve-2008-benchmark',
        'values': [
                  0,
                  0,
                  # lambda z, y, x, t: x*np.sin(y)*(-0.011125 + 0.5*(1.0 + erf((x-0.7)/0.02))*(1.0-0.92-0.2*np.power(np.cos(y), 2.0)))
                  lambda z, y, x, t: x*np.sin(y)*(0.5*(1.0 + erf((x-0.7)/0.02))*(1.0-0.92-0.2*np.power(np.cos(y), 2.0)))
                  ],
        'tvals': [1]
    }
]

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
            description='Create mean field coefficients for meanfield_e_tensor'
            )
    parser.add_argument(
            'pencilfolder',
            nargs='?',
            default=[os.getcwd()],
            help='Folder to analyze')
    args = parser.parse_args()
    pencilfolder = args.pencilfolder[0]
    coefficientcreator = CoefficientCreator(pencilfolder)

    print(80*'-' + '\n\nCreating coefficients:\n')
    for coef_conf in coef_configs:
        coefficientcreator.createCoefficient(**coef_conf)
    with h5py.File(coefficientcreator.h5file) as f:
        printH5Structure(f)
