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
    pytest.param("samples/helical-MHDturb", marks=pytest.mark.xfail(reason="broken after merging gputestv6")),
    "samples/continuous-forcing-from-file",
    pytest.param("samples/conv-slab-noequi", marks=pytest.mark.xfail(reason="broken after merging gputestv6")),
    "samples/power_xy/integrate_shell_z",
    "samples/power_xy/complex",
    "samples/power_xy/integrate_shell",
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
