#!/usr/bin/python
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

# Name:   proboscis.py
# Author: wd (wdobler [at] gmail [dot] com)
# Date:   13-Nov-2015
# Description:
#   A dummy Python module for running test+git-pc under Ubuntu Hardy.

functions = []
errors = []


def register(*args, **kwargs):
    functions.append(args[0])


def test(*args, **kwargs):
    return register


def run_all_tests():
    for f in functions:
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
