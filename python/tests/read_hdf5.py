import numpy as np
import os
import pencil as pc

datadir_novar = os.path.realpath(
    os.path.join(
        __file__, os.path.pardir, os.path.pardir, os.path.pardir, "tests", "input", "hdf5-novar"
        )
    ) #a directory without var.h5
datadir_nogrid = os.path.realpath(
    os.path.join(
        __file__, os.path.pardir, os.path.pardir, os.path.pardir, "tests", "input", "hdf5-nogrid"
        )
    ) #a directory without grid.h5

def test_read_dim_novar():
    """
    Check if dim info can be read without using var.h5 (i.e. from grid.h5)
    """
    dim = pc.read.dim(datadir=datadir_novar)

    assert dim.nx == 1152
    assert dim.nz == 288
    assert dim.precision == 'D'

def test_read_dim_nogrid():
    """
    Check if dim info can be read from var.h5 in the absence of grid.h5
    """
    dim = pc.read.dim(datadir=datadir_novar)

    assert dim.nx == 1152
    assert dim.nz == 288
    assert dim.precision == 'D'

def test_read_grid_novar():
    """
    Check if grid info can be read without using var.h5 (i.e. from grid.h5)
    """
    grid = pc.read.grid(datadir=datadir_novar, trim=True, quiet=True)

    assert grid.Lx == 16
    assert len(grid.x) == 1152
    assert len(grid.y) == 288
    assert len(grid.z) == 288
    assert grid.y[10] == -1.8541666
    assert grid.dz_tilde[30] == 0

def test_read_grid_nogrid():
    """
    Check if grid info can be read from var.h5 in the absence of grid.h5
    """
    grid = pc.read.grid(datadir=datadir_nogrid, trim=True, quiet=True)

    assert grid.Lx == 16
    assert len(grid.x) == 1152
    assert len(grid.y) == 288
    assert len(grid.z) == 288
    assert grid.y[10] == -1.8541666
    assert grid.dz_tilde[30] == 0
