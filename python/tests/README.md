$PENCIL_HOME/python/tests
=========================

# Requirements

These tests depend on the [_Pytest_](https://pytest.org/) Python package.

# Running

## Minimal tests

```sh
$PENCIL_HOME/python/tests/test-python-modules.py
```

## Full tests

While the above script only runs a minimal set of tests, the full set of tests
(including the script tests which are run by `pc_auto-test`) can be run after
changing into this directory by simply calling
```sh
pytest
```
This requires `sourceme.sh` to have been sourced in the current shell.

## Parallelization

To run tests in parallel, install the `pytest-xdist` plugin, change into this
directory, and run
```sh
pytest -n 4 --dist loadgroup
```
where 4 is the number of tests to run at a time.

## Testing with multiple Python versions

A configuration file for [_Tox_](https://tox.wiki/) is provided. Change into
this directory and run
```sh
tox
```
While you will have to install the Python executables by yourself, the required
Python packages (`scipy`, `numpy`,...) will be automatically installed by `tox`.

# Adding tests

Pytest will search for tests in the files given in the `python_files` key of
`pytest.ini`. Any function in these files whose name starts with `test` is
treated as a test.

# Code coverage

After installing the `pytest-cov` Python package, simply run
```sh
pytest --cov=pencil --cov-report=html
```
to generate a HTML code coverage report. To view the report, open
`./htmlcov/index.html` in your browser.

# Historical notes

## Why not Proboscis?

As of 2024, [_Proboscis_](https://pythonhosted.org/proboscis/) has been
unmaintained for over a decade, and is being dropped from even the most
conservative Linux distros[^1].

- [^1] <https://tracker.debian.org/news/1575216/removed-1260-8-from-unstable/>
