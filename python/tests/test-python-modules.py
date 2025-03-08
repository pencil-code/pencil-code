#!/usr/bin/python3
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

"""Run unit tests for Pencil Code Python modules."""


import argparse
import sys
import pytest
import pathlib

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()

    sys.exit(pytest.main([
        '-c',
        str(pathlib.Path(__file__).parent/"pytest.ini"),
        "-m",
        "not integration",
        ]))
