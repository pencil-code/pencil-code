#!/usr/bin/python
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

"""
A minimal replacement for the proboscis test framework (which is not
available for Ubuntu Lucid).

"""

import sys


registry = []


def register (function, *args, **kwargs):
    registry.append(function)


def test(home=None, **kwargs):
    return register


class TestProgram(object):

    def run_and_exit(self):
        for test in registry:
            print '\n********** Test: ', test, '**********'
            sys.stdout.flush()
            test()
