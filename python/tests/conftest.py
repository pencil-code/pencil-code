import os
import pathlib
import pytest
import pencil as pc
from test_utils import get_rundir

os.environ["PENCIL_HOME"] = str(pathlib.Path(__file__).parent/"../..")

#The following seems dirty, but is needed for subprocesses to be initialized with the correct Python path
os.environ["PYTHONPATH"] = str(pathlib.Path(__file__).parent/"..")


#To populate data directories that are needed for integration tests.

@pytest.fixture(scope='session')
def datadir_conv_slab_noequi():
    rundir = get_rundir("samples/conv-slab-noequi")
    sim = pc.sim.get(rundir, quiet=True)

    if sim is False:
        raise RuntimeError(f"Could not get simulation in {rundir}")

    sim.compile(bashrc=False)
    sim.run(bashrc=False, cleardata=True)

    return sim.datadir

@pytest.fixture(scope='session')
def datadir_conv_slab():
    rundir = get_rundir("samples/conv-slab")
    sim = pc.sim.get(rundir, quiet=True)

    if sim is False:
        raise RuntimeError(f"Could not get simulation in {rundir}")

    sim.compile(bashrc=False)
    sim.run(bashrc=False)

    sim.bash("pc_build -t read_all_videofiles", bashrc=False, verbose=False)
    sim.bash("src/read_all_videofiles.x", bashrc=False, verbose=False)

    return sim.datadir

