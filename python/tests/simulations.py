"""
Tests for the `Simulations` object
"""

import pathlib
import pytest
from pencil.sim import (
    Simulations,
    get_sims,
    )
from test_utils import get_rundir

samples = [
    "samples/helical-MHDturb_HDF5",
    "samples/conv-slab_cp_2",
    ]

simdirs = [get_rundir(sample) for sample in samples]

def test_Sims_from_tuple():
    sims = Simulations(*simdirs)
    assert len(sims.sims) == len(simdirs)

def test_Sims_from_list():
    sims = Simulations(simdirs)
    assert len(sims.sims) == len(simdirs)

def test_Sims_add_Sims():
    sims = Simulations(simdirs)
    sims_2 = Simulations(simdirs[0])
    sims.add(sims_2)
    assert len(sims.sims) == len(simdirs) + 1

def test_Sims_is_iterable():
    sims = Simulations(*simdirs)
    
    assert len(sims) == len(simdirs)
    
    for simdir in sims:
        pass

def test_Sims_kwarg():
    sims = Simulations(*simdirs, custom_property="something")
    assert len(sims.sims) == len(simdirs)
    assert sims.custom_property is "something"

def test_Sims_filter():
    sims = Simulations(*simdirs)
    name = pathlib.Path(simdirs[0]).name

    for sim in sims:
        if sim.name == name:
            #Later, we check if this property is preserved
            sim.custom_property = Ellipsis

    sims_f = sims.filter(lambda sim: sim.name == name)
    assert len(sims_f) == 1
    assert sims_f[0].custom_property is Ellipsis

@pytest.mark.xfail(reason="documented, but not implemented")
def test_Sims_sort():
    sims = Simulations(*simdirs)
    sims.sort()
    #TODO: would be good to actually check the results once this is implemented

def test_get_sims():
    loc = pathlib.Path(get_rundir("samples/1d-tests/conduction"))
    sims = get_sims(str(loc))
    assert len(sims) == len(list(loc.glob("*")))

@pytest.mark.xfail(reason="pathlib.Path not supported yet")
def test_get_sims_pathlib():
    loc = pathlib.Path(get_rundir("samples/1d-tests/conduction"))
    sims = get_sims(loc)
    assert len(sims) == len(list(loc.glob("*")))
