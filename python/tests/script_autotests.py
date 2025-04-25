import pytest
import subprocess
import pathlib
import os

import pencil as pc

pencil_home = pathlib.Path(__file__).parent/"../.."
samples_path = pencil_home/"samples"
samples_with_scripttests = set()
for p in samples_path.glob("**/tests/*.py"):
    samples_with_scripttests.add(p.parent.parent)

#Assumes sourceme.sh has been run in the current context.
env = {}
for var in ['PATH', 'PYTHONPATH', 'PENCIL_HOME']:
    env[var] = os.environ[var]

@pytest.mark.integration
@pytest.mark.parametrize("path", samples_with_scripttests)
def test_script_pcautotest(path):
    res = subprocess.run(
        ["pc_auto-test", "--auto-clean", "--script-tests=python", str(path)],
        env=env,
        shell=True,
        universal_newlines=True,
        stderr=subprocess.PIPE,
        )
    assert res.returncode == 0
    #Even if pc_auto-test fails, that is not reflected in its exit code
    assert res.stderr == ""
