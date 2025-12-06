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
- `sourceme.sh` should have been sourced in the current shell.

### Notes

- You do not need to install any other Python packages (they will be automatically
installed in an isolated environment).
- If multiple Python versions are present on your system, the above command will
run tests using all of them.
- Tox has also been configured to generate a code coverage report (`./htmlcov/index.html`).


## Alternative: directly calling pytest

Change into this directory and run
```sh
pytest
```
This also requires `sourceme.sh` to have been sourced in the current shell.
Moreover, the Python packages that Pencil's Python module depends
on should be installed.

To run tests in parallel, install the `pytest-xdist` plugin, change into this
directory, and run
```sh
pytest -n 4 --dist loadgroup
```
where 4 is the number of tests to run at a time.

# Adding tests

Pytest will search for tests in the files given in the `python_files` key of
`pytest.ini`. Any function in these files whose name starts with `test` is
treated as a test.

# Historical notes

## Why not Proboscis?

As of 2024, [_Proboscis_](https://pythonhosted.org/proboscis/) has been
unmaintained for over a decade, and is being dropped from even the most
conservative Linux distros[^1].

- [^1] <https://tracker.debian.org/news/1575216/removed-1260-8-from-unstable/>
