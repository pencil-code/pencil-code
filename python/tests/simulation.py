"""
Tests for the Simulation object.
Note that many aspects of the Simulation object (e.g. sim.bash, sim.run) are
already used while setting up the other tests.
"""

from test_utils import require_sample, get_rundir
from pencil.sim import Simulation

@require_sample("samples/helical-MHDturb")
def test_sim_param(datadir_helical_MHDTurb):
    simdir = get_rundir("samples/helical-MHDturb")
    sim = Simulation(path=simdir, hard=True, quiet=True)

    assert sim.param.lmagnetic is True
    assert sim.param['initaa'] == sim.param.initaa
    assert sim.param['magnetic']['initaa'] == sim.param.initaa
