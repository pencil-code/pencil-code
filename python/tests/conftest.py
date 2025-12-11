import os
import pathlib
import pytest
import shutil
import tempfile

import pencil as pc
from test_utils import (
    compile_and_run_sample,
    _require_sample_markers,
    )

pencil_home = pathlib.Path(__file__).parent.parent.parent
os.environ["PENCIL_HOME"] = str(pencil_home)

#The following seems dirty, but is needed for subprocesses to be initialized with the correct Python path
os.environ["PYTHONPATH"] = str(pencil_home/"python")

def pytest_addoption(parser):
    parser.addoption(
        "--script-test-coverage",
        action="store_true",
        help="rerun the script autotests to capture code coverage data",
        default=False,
        )

#To populate data directories that are needed for integration tests.

@pytest.fixture(scope='session')
def datadir_conv_slab_noequi():
    sim = compile_and_run_sample("samples/conv-slab-noequi")
    return sim.datadir

@pytest.fixture(scope='session')
def datadir_conv_slab():
    sim = compile_and_run_sample("samples/conv-slab")

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
    sim = compile_and_run_sample("samples/conv-slab_cp_2")
    return sim.datadir

@pytest.fixture(scope='session')
def datadir_helical_MHDTurb_HDF5():
    sim = compile_and_run_sample("samples/helical-MHDturb_HDF5")
    return sim.datadir

@pytest.fixture(scope='session')
def datadir_helical_MHDTurb():
    sim = compile_and_run_sample("samples/helical-MHDturb")
    return sim.datadir

@pytest.fixture
def temp_datadir():
    tmpdir = tempfile.mkdtemp("pc_test_write_h5_snap")
    yield pathlib.Path(tmpdir)
    shutil.rmtree(tmpdir)
