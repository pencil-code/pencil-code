#!/usr/bin/python3
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

"""Run unit tests for Pencil Code Python modules."""


import argparse
import pathlib
import sys

def call_pytest():
    #Keep this import here so that call_tox works without pytest installed.
    import pytest

    sys.exit(pytest.main([
        '-c',
        str(pathlib.Path(__file__).parent/"pytest.ini"),
        "-m",
        "not integration",
        ]))

def call_tox():
    raise NotImplementedError

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--full",
        help = "If passed, run the full set of tests, generating HTML output and a code coverage report.",
        default = False,
        action = 'store_true',
        )

    args = parser.parse_args()

    if args.full:
        call_tox()
    else:
        call_pytest()
