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

To do:
 - Remove temporary directories for successful tests
 - Between consecutive git calls, increase GIT_AUTHOR_DATE /
   GIT_COMMITTER_DATE by ≥ 1 second in order to get a more readable
   history.

'''

import os
import shutil
import subprocess
import sys
import tempfile


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
    # Is there something more basic we can test?
    # print run_system_cmd(['git', 'pc', 'checkout'])


@test(groups=['tag-wip'])
def test_tag_wip():
    '''Tag unrecorded changes with 'git pc tag-wip\''''
    git = GitSandbox('tag_wip', initial_commit=True)
    file1 = 'committed-file'

    git.write_line_to(file1, 'Committed line.')
    git('add', file1)
    git.commit_all('Committing one line.')

    git.write_line_to(file1, 'Staged line.')
    git('add', file1)

    git.write_line_to(file1, 'Uncommitted line.')

    git.write_line_to(
        'uncommitted-file',
        'Uncommitted line in uncommitted file.'
        )

    git('pc', 'tag-wip')


def run_system_cmd(cmd_line):
    '''Run a system command, writing output to the terminal'''
    print ' '.join(cmd_line)
    print '\n'.join(run_system_cmd_get_output(cmd_line))


def run_system_cmd_get_output(cmd_line):
    '''Run a system command and return output as array of lines'''
    print ' '.join(cmd_line)
    try:
        output = subprocess.check_output(cmd_line)
        return output.splitlines()
    except subprocess.CalledProcessError, e:
        print e
        sys.exit(1)


class TmpDir(object):
    '''A temporary directory.

    After successful operation, that directory normally gets removed, so
    don't leave important files there.

    '''

    def __init__(self, parent_dir=None, name='test'):
        self.path = tempfile.mkdtemp(
            suffix='', prefix=name + '_', dir=parent_dir
            )

    def purge(self):
        '''Remove everything in this temporary directory.'''
        shutil.rmtree(self.path)


class GitSandbox(object):
    '''A directory associated with a git checkout

    Usage:
      git = GitSandbox('omni-fix')
      git('commit', '-m', 'Fix all problems')
      for l in git('status'):
          print 'Git: ', s
      files = git.system_cmd('ls', '-a')

    '''

    def __init__(
            self, name,
            bare=None, initial_commit=False, root_dir=None
            ):
        dir_basename = 'git-pc-test_' + name
        if root_dir:
            self.directory = TmpDir(root_dir, dir_basename)
        else:
            self.directory = TmpDir(None, dir_basename)
        os.chdir(self.directory.path)
        if bare:
            self.__call__('init', '--bare')
        else:
            self.__call__('init')
        if initial_commit:
            self.__call__(
                'commit', '--allow-empty',
                '-m', 'Initial commit in %s' % (dir_basename, )
                )

    def purge(self):
        if self.directory:
            self.directory.purge()

    def __call__(self, *args):
        cmd_list = ['git']
        cmd_list.extend(args)
        run_system_cmd(cmd_list)

    def commit_all(self, message):
        self.__call__('commit', '-a', '-m', message)

    def system_cmd(self, *args):
        return run_system_cmd_get_output(args)

    def write_line_to(self, filename, line):
        '''Create file FILENAME if necessary, and append the given line'''
        f = open(filename, 'a')
        f.write(line + '\n')
        f.close()


if __name__ == '__main__':
    main()
