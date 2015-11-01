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
 - Reduce output

'''

import datetime
import os
import shutil
import subprocess
import sys
import tempfile


from proboscis import test, TestProgram
#from proboscis.asserts import assert_equal, assert_not_equal,assert_true, assert_false


current_git_time = 0            # Our test commits all occurred in 1970...


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


@test(groups=['panic'])
def test_panic():
    '''Test 'git pc panic\''''
    git = GitSandbox('panic', initial_commit=True)

    for f in 'file1', 'file2', 'file3':
        # Commit file
        git.write_line_to(f, 'Committed line.')
        git('add', f)
        git('commit', f, '-m', 'Committing file %s.' % (f, ))

        # Stash another change
        git.write_line_to(f, 'Stashed line.')
        git('stash')

        # Forget about file
        git('reset', '--hard', 'HEAD~')

    git('pc', 'panic', '-l')
    git('pc', 'panic', '-g')
    git('pc', 'panic', '--full', '-g')


def run_system_cmd(cmd_line):
    '''Run a system command, writing output to the terminal'''
    print ' '.join(cmd_line)
    print '\n'.join(run_system_cmd_get_output(cmd_line))


def run_system_cmd_get_output(cmd_line):
    '''Run a system command and return output as array of lines'''
    print ' '.join(cmd_line)

    # Set the commit time.
    # We do this in order to have at least one second between different
    # git commits, because otherwise 'git log' and friends often show the
    # wrong time order for commits on different branches.
    global current_git_time
    current_git_time += 1
    dtime = datetime.datetime.fromtimestamp(current_git_time)
    time_string = dtime.ctime()
    os.environ['GIT_AUTHOR_DATE'] = time_string
    os.environ['GIT_COMMITTER_DATE'] = time_string
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
