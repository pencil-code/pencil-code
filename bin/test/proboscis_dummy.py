#!/usr/bin/python
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

"""
A minimal replacement for the proboscis test framework (which is not
available for Ubuntu Lucid).

"""

import sys


functions = []
errors = []


def register(*args, **kwargs):
    functions.append(args[0])


def test(*args, **kwargs):
    return register


def run_all_tests():
    for f in functions:
        print '\n********** Test: ', f, '**********'
        sys.stdout.flush()
        try:
            f()
        except Exception, e:
            errors.append((f, e))
    if errors:
        print 'There were errors:', errors
    else:
        print 'Success'


class TestProgram(object):

    def run_and_exit(*args):
        run_all_tests()
