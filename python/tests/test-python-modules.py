#!/usr/bin/python3
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

"""Run unit tests for Pencil Code Python modules."""


import argparse

from test_utils import TestProgram

# Most likely the only place in the whole Pencil Code where importing '*'
# makes sense: Import all test cases from the following modules.
from read import *
from imports import *
from basics import *
from manual import *

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()

    TestProgram().run_and_exit()


if __name__ == "__main__":
    main()
