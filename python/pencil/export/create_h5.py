# create_h5.py
"""
Created on Fri Apr 22 11:14:28 2016

@author: fagent

Create structure for hdf5 storage of 2D averages time series, motivated
to save emf tensors for use with mean-field pencil code module.
To generalise hdf5 for all data structures
"""


def create_aver_sph(
    filename, dataset, fields, nzyxt, zyxt, hdf5dir="data", dgroup="emftensor"
):
    """
    call signature:

    create_aver_sph(filename, dataset, fields, nzyxt, zyxt,
                    hdf5dir='data', dgroup='emftensor')

    Keyword arguments:

      *filename*:
        Output filename.

      *dataset*:
        HDF5 dataset.

      *fields*:
        Names of the fields to be written.

      *nzyxt*:
        Array containing [nz, ny, nx, nt].

      *zyxt*:
        Array containing the arrays z, y, x and t.

      *hdf5dir*:
        Output data directory.

      *dgroup*:
        Data group name.
    """

    import numpy as np
    import os
    import h5py

    # Unpack the tuples.
    nz, ny, nx, nt = nzyxt
    z, y, x, t = zyxt

    # Prepare hdf5 directory.
    if not os.path.exists(hdf5dir):
        os.makedirs(hdf5dir)
    print(hdf5dir, os.path.exists(hdf5dir), "constructing hdf5 file")
    if os.path.exists(filename):
        print(filename + " already exists, skipping file creation.")
    with h5py.File(filename, "a") as hf:
        # Existing grid will not be rewritten by subsequent datasets.
        if "grid" not in hf.keys():
            hf.create_group("grid")
        if "t" not in hf["grid"].keys():
            hf.create_dataset("grid/t", (t.size,), data=np.asfarray(t))
        if "x" not in hf["grid"].keys():
            hf.create_dataset("grid/x", (nx,), data=np.asfarray(x))
        if "y" not in hf["grid"].keys():
            hf.create_dataset("grid/y", (ny,), data=np.asfarray(y))
        if "z" not in hf["grid"].keys():
            hf.create_dataset("grid/z", (nz,), data=np.asfarray(z))
        if dgroup not in hf.keys():
            hf.create_group(dgroup)
        for field, comp in fields:
            shape = comp + (nz, ny, nx, nt)
            if field not in hf[dgroup].keys():
                hf.create_group("{0}/{1}".format(dgroup, field))
            if dataset not in hf["{0}/{1}".format(dgroup, field)].keys():
                hf.create_dataset(
                    "{0}/{1}/{2}".format(dgroup, field, dataset), shape, "f"
                )
            if shape != hf["{0}/{1}/{2}".format(dgroup, field, dataset)].shape:
                del hf["{0}/{1}/{2}".format(dgroup, field, dataset)]
                hf.create_dataset(
                    "{0}/{1}/{2}".format(dgroup, field, dataset), shape, "f"
                )


fvars = [
    ("utensor", (3,)),
    (
        "alpha",
        (
            3,
            3,
        ),
    ),
    (
        "beta",
        (
            3,
            3,
        ),
    ),
    ("gamma", (3,)),
    ("delta", (3,)),
    (
        "kappa",
        (
            3,
            3,
            3,
        ),
    ),
    (
        "acoef",
        (
            3,
            3,
        ),
    ),
    (
        "bcoef",
        (
            3,
            3,
            3,
        ),
    ),
]
