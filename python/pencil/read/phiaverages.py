# phiaverages.py
#
# Read phi-average files.
"""
Contains the classes and methods to read phi-averaged files.
"""

import sys
from pencil.util import copy_docstring, PencilArray
import pathlib

from .averages import _Plane, Averages

class PhiAverages(Averages):
    """
    Old class to handle reading phi averages.
    """

    def read(
        self,
        datadir="data",
        avg_file=None,
        var_index=-1,
        iter_list=None,
        precision="f",
        param=None,
    ):
        """
         read(datadir="data", avg_file="1", var_index=-1, iter_list=None, precision="f", param=None)

         Read Pencil Code phi-averaged data.


         Keyword arguments:

        datadir : string
            Directory where the data is stored.

        avg_file : int
            Number of average file to be read.

        var_index : int
            Index of single variable taken from among the 'phi' averages.
            Takes an integer value < len(phiaver.in).

        iter_list : list of int
            Iteration indices for which to sample the slices.

        precision : string
            Float (f), double (d) or half (half).

        param     : Param object
            Param object belonging to the simulation.

        Notes
        -----
        Each returned array (e.g. self.phiavg.rho) is a
        pencil.util.PencilArray with axis_order ('t', 'z', 'r'); its
        .reorder(...) method returns a view with axes permuted by label.
        """

        import os
        import numpy as np
        from pencil import read

        datadir = pathlib.Path(datadir)

        if param is None:
            param = read.param(datadir=datadir, quiet=True)

        lh5 = (param.io_strategy == "HDF5")

        # Determine which average files to read.
        in_file_name_list = []
        aver_file_name_list = []
        if lh5:
            raise NotImplementedError("read.ogrid: not implemented for hdf5")
        else:
            if os.path.exists(datadir/"averages/phiavg.list"):
                naverages = sum(1 for line in open(datadir/"averages/phiavg.list"))
                in_file_name_list.append(datadir/"averages/phiavg.list")
                if avg_file:
                    aver_file_name_list.append(
                    os.path.join("averages", "PHIAVG" + avg_file)
                    )
                else:
                    nphiavg = sum(1 for line in open(datadir/"averages/phiavg.files"))
                    in_file_name_list = []
                    for phi_index in range(1,nphiavg+1):
                        in_file_name_list.append(datadir/"averages/phiavg.list")
                        aver_file_name_list.append(
                        os.path.join("averages", "PHIAVG" + str(phi_index))
                        )

        if len(in_file_name_list) == 0:
            raise ValueError("could not find phiavg.list")

        glob_dim = read.dim(datadir)
        nu = glob_dim.nx / 2
        nv = glob_dim.nz
        t_list = np.zeros(shape=(len(aver_file_name_list)))
        raw_data_list = PencilArray(
            np.zeros(shape=(len(aver_file_name_list),naverages,int(nv),int(nu))),
            axis_order=("t", "var", "z", "r"),
        )

        ii=0
        for in_file_name, aver_file_name in zip(in_file_name_list, aver_file_name_list):
            # This one will store the data.
            ext_object = _Plane()

            # Get the averaged quantities.
            file_id = open(os.path.join(os.path.dirname(datadir), in_file_name))
            variables = file_id.readlines()
            file_id.close()

            #remove newlines
            while True:
                try:
                    variables.remove("\n")
                except ValueError:
                    break

            n_vars = len(variables)

            t, r_cyl, z_cyl, raw_data = self._read_phiaver(
                datadir,
                variables,
                aver_file_name,
                n_vars,
                var_index,
                iter_list,
                precision=precision,
                lh5=lh5,
            )

            t_list[ii] = t
            raw_data_list[ii,:] = raw_data
            ii=ii+1

        # Add the raw data to self.
        var_idx = 0
        for var in variables:
            if var_index >= 0:
                if var_idx == var_index:
                    setattr(ext_object, var.strip(), raw_data_list[:,:, ...])
            else:
                setattr(ext_object, var.strip(), raw_data_list[:,var_idx, :, ...])
            var_idx += 1

        self.t = t_list
        self.r_cyl = r_cyl
        self.z_cyl = z_cyl
        self.phiavg = ext_object

@copy_docstring(PhiAverages.read)
def phiaver(*args, **kwargs):
    """
    Wrapper for :py:meth:`PhiAverages.read`
    """
    averages_tmp = PhiAverages()
    averages_tmp.read(*args, **kwargs)
    return averages_tmp
