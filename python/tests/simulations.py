"""
Tests for the `Simulations` object
"""

import pytest
from pencil.sim import Simulations
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
    sims = Simulations(*simdirs, quiet=True)
    assert len(sims.sims) == len(simdirs)
