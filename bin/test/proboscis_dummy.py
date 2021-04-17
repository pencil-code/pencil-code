#!/usr/bin/python
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

"""
A minimal replacement for the proboscis test framework (which is not
available for Ubuntu Lucid).

"""

from __future__ import print_function

import sys

from typing import Any, Callable, Optional, TypeVar

F = TypeVar('F', bound=Callable[..., Any])


functions = []
errors = []




def register(function: F, **kwargs: Any) -> F:
    functions.append(function)
    return function


def test(function: F, **kwargs: Any) -> F:
    return register(function, **kwargs)


def run_all_tests() -> None:
    for f in functions:
        print("Testing {} ... ".format(_identify(f)), end="")
        sys.stdout.flush()
        try:
            f()
            print("ok")
        except Exception as e:
            errors.append((f, e))
    if errors:
        print("There were errors:", errors)
    else:
        print("Success")


def _identify(f: F) -> str:
    if f.__doc__ and len(f.__doc__) > 0:
        return f.__doc__.splitlines()[0]
    else:
        return f.__name__

def assert_equal(expected:Any , actual: Any, message: Optional[str]=None) -> None:
    if message is None:
        message = "{} ≠ {}".format(expected, actual)
    assert expected == actual, message


def assert_true(actual: bool, message: Optional[str]=None) -> None:
    if message is None:
        assert bool(actual)
    else:
        assert bool(actual), message


class TestProgram(object):
    """The class that runs the tests."""

    def run_and_exit(*args: Any) -> None:
        run_all_tests()
