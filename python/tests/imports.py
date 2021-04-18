#!/usr/bin/python3
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

"""Test Python imports.

Note: Given that Python tries to avoid importing the same module twice,
any repeated imports listed here (e.g. 'from pencil.sim import
simulation' followed by 'from pencil.sim.get', where get.py imports
simulation.py) may actually test a lot less than one would expect.
  It would be better to start a new Python process for each import, but
currently the latency of importing pencil is many seconds, so these
tests would take a significant amount of time.

Note 2: The tests here are used to ensure that functionality is not quietly
removed when refactoring the Python modules.
  However, being listed here does not automatically declare an import
meaningful and important. Any import listed here that does not make sense,
can be removed in the code, provided these tests get adapted.

"""

from test_utils import test, standalone_test


@test
def import_stuff() -> None:
    """Import a number of modules.

    Ideally, this list will grow each time we refactor something.

    """
    from pencil.sim import get  # noqa
    from pencil.diag.particle import dispersion_and_drift  # noqa
    from pencil.diag.particle import diffusion # noqa
    from pencil.diag.particle import gas_velo_at_particle_pos # noqa
    from pencil.io import remove_files # noqa
    from pencil.io import fort2h5 # noqa
    from pencil.io import snapshot # noqa
    from pencil.sim import remesh # noqa
    from pencil.sim import simulation # noqa
    from pencil import sim # noqa

    from pencil.util import is_sim_dir # noqa
    from pencil.calc import Reynolds # noqa
    from pencil import calc # noqa
    from pencil.calc import accuracy # noqa
    from pencil.calc import draglift # noqa
    from pencil.calc import tensors # noqa
    from pencil.calc import __init__ # noqa
    from pencil.calc import shocktube # noqa
    from pencil.export import pc2vtk # noqa
    from pencil.io import get_value_from_file # noqa
    from pencil.io import mkdir # noqa
    from pencil.io import npfile # noqa

    from pencil.io import pc_hdf5 # noqa
    from pencil.ism_dyn import derived_h5 # noqa
    from pencil.ism_dyn import get_masks # noqa
    from pencil.ism_dyn import get_stats # noqa
    from pencil.ism_dyn import ism_cooling # noqa
    from pencil.ism_dyn import ism_sedov_taylor # noqa
    from pencil.math import Helmholtz # noqa
    from pencil.math.derivatives import der_6th_order_w_ghosts # noqa
    from pencil.math.derivatives import div_grad_curl # noqa

    from pencil.math import primes # noqa
    from pencil.read import allslices # noqa
    from pencil.read import averages # noqa
    from pencil.read import indices # noqa
    from pencil.read import params # noqa
    from pencil.read import power # noqa
    from pencil.read import pvarfile # noqa
    from pencil.read import varfile # noqa
    from pencil.sim import group # noqa
    from pencil.calc import streamlines # noqa
    from pencil.backpack import pidly # noqa

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
    """Import the pencil module in a separate Python process.

    This will in turn import all other modules, including visu.

    """
    standalone_test(["import pencil"])
