$PENCIL_HOME/python/tests
=========================

# Running minimal tests

Run
```sh
$PENCIL_HOME/python/tests/test-python-modules.py
```
This option requires the [_Pytest_](https://pytest.org/) Python package to be
installed. Moreover, the Python packages that Pencil's Python module depends
on should be installed.

# Running full tests

While the above script only runs a minimal set of tests, the full set of tests
(including the script tests which are run by `pc_auto-test`) can be run in either
of the ways described below.

## Recommended

```sh
$PENCIL_HOME/python/tests/test-python-modules.py --full
```

### Requirements

- A Python interpreter should be in `$PATH` (e.g. `/usr/bin/python3`).
- [_Tox_](https://tox.wiki/) (Python package) should be installed.
- For the code coverage report, [_Coverage_](https://pypi.org/project/coverage/) (Python package) should also be installed.

### Notes

- You do not need to install any other Python packages (they will be automatically
installed in an isolated environment).
- If multiple Python versions are present on your system, the above command will
run tests using all of them.
- A HTML report, `report.html`, will be generated, with a link to the code coverage report.
- Use the `--outputdir` option to change the location of the reports.


## Alternative: directly calling pytest

Change into this directory and run
```sh
pytest
```
This requires `sourceme.sh` to have been sourced in the current shell.
Moreover, the Python packages that Pencil's Python module depends
on should be installed.

To run tests in parallel, install the `pytest-xdist` plugin, change into this
directory, and run
```sh
pytest -n 4 --dist loadgroup
```
where 4 is the number of tests to run at a time.

# Where can I add tests?

There are two kinds of tests, described below.

## Standalone tests of the Python module

These are executed by `test-python-modules`.

### Where are they defined?

Pytest will search for tests in the files given in the `python_files` key of
`pytest.ini`. Any function in these files whose name starts with `test` is
treated as a test.

### Using the output of a sample for your tests

If your test depends on the output of a sample, you can define a fixture that
will run Pencil before running the Python test. An example is the function
`datadir_helical_MHDTurb` in `conftest.py`. Such tests should be marked with the
`requires_sample` decorator, e.g.
```
@require_sample("samples/helical-MHDturb")
def test_read_var_2_trim(datadir_helical_MHDTurb):
	#test some stuff
	assert 1 == 1
```

## Tests that use pc_auto-test

There is also a set of tests that can be run by `pc_auto-test` by using its
`--script-tests=python` option. Alternatively, passing `-full` to
`test-python-modules` causes it to execute these tests in in addition to the
standalone tests described earlier.

### Where are they defined?

These tests are defined by Python scripts that
live in a subdirectory `tests` of the corresponding sample. Such tests are
expected to write output into a text file, which is then compared with a
reference file. For an example, see `samples/helical-MHDturb/tests/read_data.py`

### When should I use this option?

Use this option for tests that use the Python module to test the
output of the Fortran code (e.g. power spectra, averages).

# Historical notes

## Why not Proboscis?

As of 2024, [_Proboscis_](https://pythonhosted.org/proboscis/) has been
unmaintained for over a decade, and is being dropped from even the most
conservative Linux distros[^1].

- [^1] <https://tracker.debian.org/news/1575216/removed-1260-8-from-unstable/>
