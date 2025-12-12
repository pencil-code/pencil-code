import pytest
import subprocess
import pathlib
import os

import pencil as pc

from test_utils import (
    get_rundir,
    _require_sample_markers,
    )

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

def add_marks_to_sample(sample):
    """
    Adds the datadir dependency information for pytest-xdist.
    """
    if isinstance(sample, (pathlib.Path, str)):
        path = sample
        marks = []
    else:
        [path] = sample.values
        marks = sample.marks

    marks.extend(_require_sample_markers(path))
    return pytest.param(path, marks=marks)

_samples_with_scripttests = [add_marks_to_sample(sample) for sample in samples_with_scripttests]

#Assumes sourceme.sh has been run in the current context.
env = {}
for var in [
    "PATH",
    "PYTHONPATH",
    "PENCIL_HOME",
    "HOME",
    ]:
    env[var] = os.environ[var]

@pytest.mark.pcautotest
@pytest.mark.parametrize("path", map(add_marks_to_sample, samples_with_scripttests))
def test_script_pcautotest(path, pytestconfig):
    rundir = pathlib.Path(get_rundir(path))
    res = subprocess.run(
        f"pc_auto-test --no-lock --auto-clean --script-tests=python '{rundir}'",
        env=env,
        shell=True,
        universal_newlines=True,
        stderr=subprocess.PIPE,
        )

    if pytestconfig.option.script_test_coverage:
        for script in rundir.glob("tests/*.py"):
            subprocess.run(
                ["python", f"{script}"],
                capture_output=True,
                cwd=rundir/"tests",
                )

    assert res.returncode == 0
    #Even if pc_auto-test fails, that is not reflected in its exit code
    assert res.stderr == ""
