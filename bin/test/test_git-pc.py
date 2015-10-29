#!/usr/bin/python
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

'''Unit tests for the 'git pc' command.

These tests mostly check whether a 'git pc' subcommand runs without
throwing an error.
It would be nice (but much more work) to check for correctness of
operation.

For best results, install the packages
 - python-proboscis,
 - python-nose
and run these tests with

  ${PENCIL_HOME}/bin/test/test_git-pc.py

'''

import subprocess

from proboscis import test, TestProgram
#from proboscis.asserts import assert_equal, assert_not_equal,assert_true, assert_false


def main():

    TestProgram().run_and_exit()


@test()
def git_pc_in_path():
    '''Make sure 'git pc' is found'''
    run_system_cmd(['git', 'pc', '-h'])


@test(groups=['calling'])
def call_checkout():
    '''Run 'git pc checkout' (not really)'''
    # We don't really want to do a full checkout.
    # Is there something more baseic we can test?
    # print run_system_cmd(['git', 'pc', 'checkout'])


def run_system_cmd(cmd_list):
    '''Run a system command.

    We need to explicitly read and then print the output, in order to make
    nose's output capturing work.

    '''
    output = subprocess.check_output(cmd_list)
    print output


if __name__ == '__main__':
    main()
