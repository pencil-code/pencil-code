import numpy as np
import os
import pencil as pc
import pytest

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

@pytest.mark.integration
def test_read_var(datadir_helical_MHDTurb_HDF5):
    var = pc.read.var(
        datadir=datadir_helical_MHDTurb_HDF5,
        trimall=True,
        lpersist=True,
        magic=["bb"],
        )

    assert len(var.x) == 32
    assert len(var.y) == 32
    assert len(var.z) == 32
    assert np.isclose(var.f[0,3,5,6], 0.02334117174211011)
    assert np.isclose(var.uy[0,6,2], -0.05910974500656841)
    assert np.isclose(var.uz[5,10,27], -0.02635602018447831)
    assert np.isclose(var.persist.forcing_tsforce, 0.3999999999999999)

@pytest.mark.integration
def test_read_var_irangex(datadir_helical_MHDTurb_HDF5):
    var = pc.read.var(datadir=datadir_helical_MHDTurb_HDF5, trimall=False, irange_x=[15,30])

    assert len(var.x) == 22
    assert len(var.y) == 38
    assert len(var.z) == 38
    assert np.isclose(var.uz[8,13,18], -0.02635602018447831)

@pytest.mark.integration
def test_read_var_irangex_trim(datadir_helical_MHDTurb_HDF5):
    var = pc.read.var(datadir=datadir_helical_MHDTurb_HDF5, trimall=True, irange_x=[15,30])

    assert len(var.x) == 16
    assert len(var.y) == 32
    assert len(var.z) == 32
    assert np.isclose(var.uz[5,10,15], -0.02635602018447831)

@pytest.mark.integration
def test_read_var_rangex_trim(datadir_helical_MHDTurb_HDF5):
    var = pc.read.var(datadir=datadir_helical_MHDTurb_HDF5, trimall=True, range_x=[-0.69,2.26])

    assert len(var.x) == 16
    assert len(var.y) == 32
    assert len(var.z) == 32
    assert np.isclose(var.uz[5,10,15], -0.02635602018447831)

@pytest.mark.integration
def test_read_var_irangexslice(datadir_helical_MHDTurb_HDF5):
    var = pc.read.var(datadir=datadir_helical_MHDTurb_HDF5, trimall=False, irange_x=slice(15,30))

    assert len(var.x) == 22
    assert len(var.y) == 38
    assert len(var.z) == 38
    assert np.isclose(var.uz[8,13,18], -0.02635602018447831)
