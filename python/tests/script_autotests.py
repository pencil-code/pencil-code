import pytest
import subprocess
import pathlib
import os

import pencil as pc

from test_utils import get_rundir

samples_with_scripttests = [
    #each entry is a string, the path relative to PENCIL_HOME
    #to mark a failing test, use
    ##pytest.param(path, marks=pytest.mark.xfail)
    "samples/continuous-forcing-from-file",
    "samples/conv-slab-noequi",
    "samples/helical-MHDturb",
    "samples/helical-MHDturb_HDF5",
    "samples/power_xy/complex_iodist",
    "samples/power_xy/complex_nprocz_2",
    "samples/power_xy/complex_lpowerxyhdf5/nprocx_2",
    "samples/power_xy/complex_lpowerxyhdf5/nprocy_2",
    "samples/power_xy/integrate_shell",
    "samples/power_xy/complex",
    "samples/power_xy/integrate_shell_z",
    ]

#Assumes sourceme.sh has been run in the current context.
env = {}
for var in [
    "PATH",
    "PYTHONPATH",
    "PENCIL_HOME",
    "HOME",
    ]:
    env[var] = os.environ[var]

@pytest.mark.integration
@pytest.mark.pcautotest
@pytest.mark.parametrize("path", samples_with_scripttests)
def test_script_pcautotest(path):
    rundir = get_rundir(path)
    res = subprocess.run(
        f"pc_auto-test --auto-clean --script-tests=python '{rundir}'",
        env=env,
        shell=True,
        universal_newlines=True,
        stderr=subprocess.PIPE,
        )
    assert res.returncode == 0
    #Even if pc_auto-test fails, that is not reflected in its exit code
    assert res.stderr == ""
