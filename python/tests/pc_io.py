#!/usr/bin/python3
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

"""
Test functions in pc.io
"""

import os
from pencil.io import get_value_from_file
from test_utils import make_test, assert_equal, get_rundir


def test_get_value_from_file() -> None:
    """Reading values from {start,run}.in"""
    rundir = get_rundir("samples/conv-slab")

    def get_s(quantity):
        return get_value_from_file("start.in", quantity, filepath=rundir)

    def get_r(quantity):
        return get_value_from_file("run.in", quantity, filepath=rundir)

    assert_equal(get_s("gravz"), -1)
    assert_equal(get_s("gravz_profile"), "const")
    assert_equal(get_s('bcz'), ['s','s','a','a2','a2:cT'])
    assert_equal(get_r("wcool"), 0.2)
    assert_equal(get_r("cool"), 15)
    assert_equal(get_r("gravz"), -1)
    assert_equal(get_r("gravz_profile"), "const")
    assert_equal(get_r("lupw_lnrho"), True)
