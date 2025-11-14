import os
import pathlib
import pytest
import shutil
import tempfile

import pencil as pc
from test_utils import (
    get_rundir,
    _require_sample_markers,
    )

pencil_home = pathlib.Path(__file__).parent.parent.parent
os.environ["PENCIL_HOME"] = str(pencil_home)

#The following seems dirty, but is needed for subprocesses to be initialized with the correct Python path
os.environ["PYTHONPATH"] = str(pencil_home/"python")


#To populate data directories that are needed for integration tests.

@pytest.fixture(scope='session')
def datadir_conv_slab_noequi():
    rundir = get_rundir("samples/conv-slab-noequi")
    sim = pc.sim.get(rundir, quiet=True)

    if sim is False:
        raise RuntimeError(f"Could not get simulation in {rundir}")

    sim.compile(
        bashrc=False,
        cleanall=False,
        autoclean=True,
        raise_errors=True,
        )
    sim.run(bashrc=False, cleardata=True, raise_errors=True)

    return sim.datadir

@pytest.fixture(scope='session')
def datadir_conv_slab():
    rundir = get_rundir("samples/conv-slab")
    sim = pc.sim.get(rundir, quiet=True)

    if sim is False:
        raise RuntimeError(f"Could not get simulation in {rundir}")

    sim.compile(
        bashrc=False,
        cleanall=False,
        autoclean=True,
        raise_errors=True,
        )
    sim.run(bashrc=False, cleardata=True, raise_errors=True)

    sim.bash(
        "pc_build -t read_all_videofiles",
        bashrc=False,
        verbose=False,
        raise_errors=True,
        )
    sim.bash(
        "src/read_all_videofiles.x",
        bashrc=False,
        verbose=False,
        raise_errors=True,
        )

    return sim.datadir

@pytest.fixture(scope='session')
def datadir_conv_slab_cp_2():
    rundir = get_rundir("samples/conv-slab_cp_2")
    sim = pc.sim.get(rundir, quiet=True)

    if sim is False:
        raise RuntimeError(f"Could not get simulation in {rundir}")

    sim.compile(
        bashrc=False,
        cleanall=False,
        autoclean=True,
        raise_errors=True,
        )
    sim.run(bashrc=False, cleardata=True, raise_errors=True)

    return sim.datadir

@pytest.fixture(scope='session')
def datadir_helical_MHDTurb_HDF5():
    rundir = get_rundir("samples/helical-MHDturb_HDF5")
    sim = pc.sim.get(rundir, quiet=True)

    if sim is False:
        raise RuntimeError(f"Could not get simulation in {rundir}")

    sim.compile(
        bashrc=False,
        cleanall=False,
        autoclean=True,
        raise_errors=True,
        )
    sim.run(bashrc=False, cleardata=True, raise_errors=True)

    return sim.datadir

@pytest.fixture(scope='session')
def datadir_helical_MHDTurb():
    rundir = get_rundir("samples/helical-MHDturb")
    sim = pc.sim.get(rundir, quiet=True)

    if sim is False:
        raise RuntimeError(f"Could not get simulation in {rundir}")

    sim.compile(
        bashrc=False,
        cleanall=False,
        autoclean=True,
        raise_errors=True,
        )
    sim.run(bashrc=False, cleardata=True, raise_errors=True)

    return sim.datadir

@pytest.fixture
def temp_datadir():
    tmpdir = tempfile.mkdtemp("pc_test_write_h5_snap")
    yield pathlib.Path(tmpdir)
    shutil.rmtree(tmpdir)
