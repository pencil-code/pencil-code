import numpy as np
import os
import pencil as pc
import pytest

from test_utils import require_sample

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

@require_sample("samples/helical-MHDturb_HDF5")
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

@require_sample("samples/helical-MHDturb_HDF5")
def test_read_var_selective(datadir_helical_MHDTurb_HDF5):
    var = pc.read.var(
        datadir=datadir_helical_MHDTurb_HDF5,
        trimall=True,
        lpersist=True,
        var_list=['uy', 'uz'],
        )

    assert len(var.x) == 32
    assert len(var.y) == 32
    assert len(var.z) == 32
    assert np.isclose(var.uy[0,6,2], -0.05910974500656841)
    assert np.isclose(var.uz[5,10,27], -0.02635602018447831)
    assert np.isclose(var.persist.forcing_tsforce, 0.3999999999999999)

    assert not hasattr(var, 'lnrho')
    assert not hasattr(var, 'ss')
    assert not hasattr(var, 'ux')
    assert var.f.shape[0] == 2

@require_sample("samples/helical-MHDturb_HDF5")
def test_read_var_irangex(datadir_helical_MHDTurb_HDF5):
    var = pc.read.var(datadir=datadir_helical_MHDTurb_HDF5, trimall=False, irange_x=[15,30])

    assert len(var.x) == 22
    assert len(var.y) == 38
    assert len(var.z) == 38
    assert np.isclose(var.uz[8,13,18], -0.02635602018447831)

@require_sample("samples/helical-MHDturb_HDF5")
def test_read_var_irangex_trim(datadir_helical_MHDTurb_HDF5):
    var = pc.read.var(datadir=datadir_helical_MHDTurb_HDF5, trimall=True, irange_x=[15,30])

    assert len(var.x) == 16
    assert len(var.y) == 32
    assert len(var.z) == 32
    assert np.isclose(var.uz[5,10,15], -0.02635602018447831)

@require_sample("samples/helical-MHDturb_HDF5")
def test_read_var_rangex_trim(datadir_helical_MHDTurb_HDF5):
    var = pc.read.var(datadir=datadir_helical_MHDTurb_HDF5, trimall=True, range_x=[-0.69,2.26])

    assert len(var.x) == 16
    assert len(var.y) == 32
    assert len(var.z) == 32
    assert np.isclose(var.uz[5,10,15], -0.02635602018447831)

@require_sample("samples/helical-MHDturb_HDF5")
def test_read_var_irangexslice(datadir_helical_MHDTurb_HDF5):
    var = pc.read.var(datadir=datadir_helical_MHDTurb_HDF5, trimall=False, irange_x=slice(15,30))

    assert len(var.x) == 22
    assert len(var.y) == 38
    assert len(var.z) == 38
    assert np.isclose(var.uz[8,13,18], -0.02635602018447831)

@require_sample("samples/power_xy/complex_lpowerxyhdf5/nprocx_2")
def test_read_power_lazy(datadir_pxyhdf5_complex_npx_2):
    """
    Reference values here are the same as in
    samples/power_xy/complex_lpowerxyhdf5/nprocx_2/tests/read_powerxy.py
    where lazy=False
    """
    p = pc.read.power(datadir=datadir_pxyhdf5_complex_npx_2, lazy=True)

    assert np.all(np.isclose(
        p.t[:5],
        [0.0543473, 0.101901, 0.156246, 0.203799, 0.251352],
        atol=1e-2,
        rtol=0,
        ))

    assert p.uz_xy.ndim == 4
    assert p.uz_xy.shape == (5,32,32,32)
    assert not np.any(np.isnan(p.uz_xy[()]))

    def pow_is_close(val, ref):
        return np.all(np.isclose(val, ref, rtol=1.5e-2, atol=0))

    assert pow_is_close(
        np.real(p.uz_xy[0,16,1,:3]),
        [-1.2e-05, -2.36e-05, 9.05e-06],
        )
    assert pow_is_close(
        np.real(p.uz_xy[3,7:15,1,1]),
        [-0.000991, -0.00108, -0.00111, -0.00118, -0.00119, -0.00118, -0.00131, -0.00138],
        )
    assert pow_is_close(
        np.real(p.uz_xy[3,16:20,1,1]),
        [-5.63e-04, -1.90e-04, -3.67e-05, -2.66e-06],
        )
    assert pow_is_close(
        np.abs(p.uz_xy[3,9,26:,1]),
        [8.73e-06, 5.53e-06, 6.00e-05, 5.03e-04, 9.30e-04, 1.81e-03],
        )
    assert pow_is_close(
        np.imag(p.uz_xy[3,9,29,-5:]),
        [2.44e-06, -9.59e-06, -6.65e-06, -4.27e-05, 4.65e-04],
        )
