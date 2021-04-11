#!/usr/bin/python
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

"""
A minimal replacement for the proboscis test framework (which is not
available for Ubuntu Lucid).

"""

from __future__ import print_function

import sys

from typing import Any, Callable


functions = []
errors = []


def register(*args, **kwargs) -> None:
    functions.append(args[0])


def test(*args: Callable[..., None], **kwargs) -> Callable[..., None]:
    return register


def run_all_tests():
    for f in functions:
        print("\n********** Test: ", f, "**********")
        sys.stdout.flush()
        try:
            f()
        except Exception as e:
            errors.append((f, e))
    if errors:
        print("There were errors:", errors)
    else:
        print("Success")


def assert_equal(expected, actual, message=None):
    if message is None:
        message = "{} ≠ {}".format(expected, actual)
    assert expected == actual, message


def assert_true(actual, message=None):
    if message is None:
        assert bool(actual)
    else:
        assert bool(actual), message


class TestProgram(object):
    def run_and_exit(*args):
        run_all_tests()
