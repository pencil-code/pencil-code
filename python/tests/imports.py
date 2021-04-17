#!/usr/bin/python3
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

"""Test Python imports.

Note: Given that Python tries to avoid importing the same module twice,
any repeated imports listed here (e.g. 'from pencil.sim import
simulation' followed by 'from pencil.sim.get', where get.py imports
simulation.py) may actually test a lot less than one would expect.
  It would be better to start a new Python process for each import, but
currently the latency of importing pencil is about 4 seconds, so these
tests would take a significant amount of time.

Note 2: The tests here are used to ensure that functionality is not quietly
removed when refactoring the Python modules.
  However, being listed here does not automatically declare an import
meaningful and important. Any import listed here that does not make sense,
can be removed in the code, provided these tests get adapted.

"""

try:
    from proboscis import test
except ImportError:
    from proboscis_dummy import test


@test
def import_stuff() -> None:
    """Import a number of modules.

    Ideally, this list will grow each time we refactor something.

    """
    from pencil.sim import get
    from pencil.diag.particle import dispersion_and_drift
    from pencil.diag.particle import diffusion
    from pencil.diag.particle import gas_velo_at_particle_pos
    from pencil.io import remove_files
    from pencil.io import fort2h5
    from pencil.io import snapshot
    from pencil.sim import remesh
    from pencil.sim import simulation
    from pencil import sim

    from pencil.util import is_sim_dir
    from pencil.calc import Reynolds
    from pencil import calc
    from pencil.calc import accuracy
    from pencil.calc import draglift
    from pencil.calc import tensors
    from pencil.calc import __init__
    from pencil.calc import shocktube
    from pencil.export import pc2vtk
    from pencil.io import get_value_from_file
    from pencil.io import mkdir
    from pencil.io import npfile

    from pencil.io import pc_hdf5
    from pencil.ism_dyn import derived_h5
    from pencil.ism_dyn import get_masks
    from pencil.ism_dyn import get_stats
    from pencil.ism_dyn import ism_cooling
    from pencil.ism_dyn import ism_sedov_taylor
    from pencil.math import Helmholtz
    from pencil.math.derivatives import der_6th_order_w_ghosts
    from pencil.math.derivatives import div_grad_curl

    from pencil.math import primes
    from pencil.read import allslices
    from pencil.read import averages
    from pencil.read import indices
    from pencil.read import params
    from pencil.read import power
    from pencil.read import pvarfile
    from pencil.read import varfile
    from pencil.sim import group
    from pencil.calc import streamlines
    from pencil.backpack import pidly

    # from pencil.calc import example_shocktube
    # from pencil import util
    pass


@test
def import_visu_stuff() -> None:
    """Import modules (or symbols?) related to visu. [commented out]

    These tests may fail, because visu is often not installed.

    """
    # from pencil.visu.internal import export_fig
    # from pencil.visu.internal import prepare_fig
    # from pencil.visu.lic_internal import lic
    # from pencil.visu import animate_multislices
    # from pencil.visu import animate_slices
    # from pencil.visu import animate_slices_maketomovie
    # from pencil.visu.internal import MinorSymLogLocator
    # from pencil.visu.internal import calc_lims
    # from pencil.visu.lic_internal import lic_demo
    # from pencil.visu import rvid_box
    pass


@test
def import_pencil() -> None:
    """Import the pencil module currently [commented out].

    This will in turn import all other modules, including visu.

    """
    # import pencil
    pass
