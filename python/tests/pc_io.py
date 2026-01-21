#!/usr/bin/python3
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

"""
Test functions in pc.io
"""

import os
import shutil
import pytest
import numpy as np

from pencil.io import (
    get_value_from_file,
    write_h5_snapshot,
    write_snapshot,
    )
from pencil import read
from test_utils import (
    make_test,
    assert_equal,
    get_rundir,
    require_sample,
    )


def test_get_value_from_file() -> None:
    """Reading values from {start,run}.in"""
    rundir = get_rundir("samples/conv-slab")

    def get_s(quantity):
        return get_value_from_file("start.in", quantity, filepath=rundir)

    def get_r(quantity):
        return get_value_from_file("run.in", quantity, filepath=rundir)

    assert_equal(get_s("gravz"), -1)
    assert_equal(get_s("gravz_profile"), "const")
    assert_equal(get_s('bcz'), ['s','s','a','a2','a2:cT'])
    assert_equal(get_r("wcool"), 0.2)
    assert_equal(get_r("cool"), 15)
    assert_equal(get_r("gravz"), -1)
    assert_equal(get_r("gravz_profile"), "const")
    assert_equal(get_r("lupw_lnrho"), True)

@require_sample("samples/helical-MHDturb")
def test_write_snapshot(datadir_helical_MHDTurb, temp_datadir):
    src_datadir = datadir_helical_MHDTurb
    dest_datadir = shutil.copytree(
        src_datadir,
        temp_datadir/"data",
        symlinks=True,
        )

    for f in dest_datadir.glob("proc*/var.dat"):
        f.unlink()

    src_var = read.var(datadir=src_datadir, trimall=True)
    src_dim = read.dim(datadir=src_datadir)

    if src_dim.precision == "D":
        precision = "d"
    else:
        precision = "f"

    dest_dim = read.dim(datadir=dest_datadir)
    write_snapshot(
        src_var.f,
        file_name = "var.dat",
        datadir = dest_datadir,
        nprocx = dest_dim.nprocx,
        nprocy = dest_dim.nprocy,
        nprocz = dest_dim.nprocz,
        precision = precision,
        nghost=3,
        )

    dest_var = read.var(datadir=dest_datadir, trimall=True)
    assert np.all(src_var.f == dest_var.f)

@require_sample("samples/helical-MHDturb_HDF5")
def test_write_h5_snapshot(datadir_helical_MHDTurb_HDF5, temp_datadir):
    src_datadir = datadir_helical_MHDTurb_HDF5
    dest_datadir = shutil.copytree(
        src_datadir,
        temp_datadir/"data",
        symlinks=True,
        )

    dest_varfile = dest_datadir/"allprocs/var.h5"
    dest_varfile.unlink()

    src_var = read.var(datadir=src_datadir, trimall=True)
    src_dim = read.dim(datadir=src_datadir)

    if src_dim.precision == "D":
        precision = "d"
    else:
        precision = "f"

    write_h5_snapshot(
        src_var.f,
        file_name = "var.h5",
        precision = precision,
        nghost=3,
        lghosts=False,
        sim_datadir=dest_datadir,
        datadir = dest_datadir/"allprocs",
        )

    dest_var = read.var(datadir=dest_datadir, trimall=True)
    assert np.all(src_var.f == dest_var.f)

@require_sample("samples/cartesian-convection-kramers-chimax_HDF5")
def test_write_h5_snapshot_nc(datadir_car_conv_kram_chimax_h5, temp_datadir):
    """
    Read snapshot from non-cubical grid (HDF5 case)
    """
    src_datadir = datadir_car_conv_kram_chimax_h5
    dest_datadir = shutil.copytree(
        src_datadir,
        temp_datadir/"data",
        symlinks=True,
        )

    dest_varfile = dest_datadir/"allprocs/var.h5"
    dest_varfile.unlink()

    src_var = read.var(datadir=src_datadir, trimall=True)
    src_dim = read.dim(datadir=src_datadir)

    if src_dim.precision == "D":
        precision = "d"
    else:
        precision = "f"

    write_h5_snapshot(
        src_var.f,
        file_name = "var.h5",
        precision = precision,
        nghost=3,
        lghosts=False,
        sim_datadir=dest_datadir,
        datadir = dest_datadir/"allprocs",
        )

    dest_var = read.var(datadir=dest_datadir, trimall=True)
    assert np.all(src_var.f == dest_var.f)
