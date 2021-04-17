#!/usr/bin/python3
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

"""Run unit tests for Pencil Code Python modules."""


import argparse

try:
    from proboscis import TestProgram
except ImportError:
    print("These tests work best with Proboscis installed:")
    print("  pip3 install proboscis")
    print("Continuing with dummy implementation")
    from proboscis_dummy import TestProgram

# Most likely the only place in the whole Pencil Code where importing '*'
# makes sense: Import all test cases from the following modules.
from read import *
from imports import *
from basics import *


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()

    TestProgram().run_and_exit()


if __name__ == "__main__":
    main()
