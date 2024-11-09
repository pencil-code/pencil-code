$PENCIL_HOME/python/tests
=========================

# Introduction

Test Pencil Code Python modules, so we can feel better when refactoring the Python code.

```sh
$PENCIL_HOME/python/tests/test-python-modules.py
```

These tests depend on the [_Pytest_](https://pytest.org/) Python package.

# Adding tests

Pytest will search for tests in the files given in the `python_files` key of `pytest.ini`.
Any function in these files whose name starts with `test` is treated as a test.

# Testing with multiple Python versions

A configuration file for [_Tox_] is provided. Change into this directory and run

```sh
tox
```

# Historical notes

## Why not Proboscis?

As of 2024, [_Proboscis_](https://pythonhosted.org/proboscis/) has been unmaintained for over a decade, and is being dropped from even the most conservative Linux distros[^1].

- [^1] <https://tracker.debian.org/news/1575216/removed-1260-8-from-unstable/>
