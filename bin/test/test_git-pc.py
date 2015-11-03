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
 - Capture stderr ('git pc panic --full' writes to stderr, although maybe
   it shouldn't), but separately from stdout.

'''

import datetime
import os
import shutil
import subprocess
import sys
import tempfile


from proboscis import test, TestProgram
#from proboscis.asserts import assert_equal, \
#                              assert_not_equal, \
#                              assert_true, \
#                              assert_false


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


@test(groups=['ff-update'])
def test_ff_update():
    '''Test 'git pc ff-update\''''
    (server, git1, git2) = setup_git_with_server('ff-update')

    # git1: work on feature branch
    git1('checkout', '-b', 'feature-branch')
    git1.write_line_to('file1', 'Line added locally on feature branch')
    git1('add', 'file1')
    git1.commit_all('Commit file on feature branch')

    # git2: work on master
    git2.write_line_to('file2', 'Line added remotely on master')
    git2('add', 'file2')
    git2.commit_all('Commit file on master')
    git2('push')

    # git1: update master without checking it out
    git1('fetch')
    git1('pc', 'ff-update', 'master')


@test(groups=['reverse-merge'])
def test_reverse_merge():
    '''Test 'git pc reverse_merge\''''
    (server, git1, git2) = setup_git_with_server('reverse_merge')

    # git1: commit a change
    git1.write_line_to('file1', 'Line added locally')
    git1('add', 'file1')
    git1.commit_all('git1: Commit file locally')

    # git2: commit a change and push
    git2.write_line_to('file2', 'Line added on server')
    git2('add', 'file2')
    git2.commit_all('git2: Commit file and push to server')
    git2('push')

    # git1: fetch and reverse-merge
    git1('fetch')
    git1('pc', 'reverse-merge', 'master@{u}')


def run_system_cmd(cmd_line, dir=None):
    '''Run a system command, writing output to the terminal'''
    print ' '.join(cmd_line)
    print '\n'.join(run_system_cmd_get_output(cmd_line, dir))


def run_system_cmd_get_output(cmd_line, dir=None):
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
        pwd = os.getcwd()
        if dir:
            os.chdir(dir)
        output = subprocess.check_output(cmd_line)
        os.chdir(pwd)
        return output.splitlines()
    except subprocess.CalledProcessError, e:
        print e
        sys.exit(1)


class TmpDir(object):
    '''A temporary directory.

    After successful operation, that directory normally gets removed, so
    don't leave important files there.

    '''

    def __init__(self, parent_dir=None, name='test', suffix=''):
        self.path = tempfile.mkdtemp(
            suffix=suffix, prefix=name + '_', dir=parent_dir
            )

    def purge(self):
        '''Remove everything in this temporary directory.'''
        shutil.rmtree(self.path)


def setup_git_with_server(name, root_dir=None):
    '''Set up a server repo with two sandbox repos'''
    dir_basename = 'git-pc-test_' + name
    if root_dir:
        top_dir = TmpDir(root_dir, dir_basename)
    else:
        top_dir = TmpDir(None, dir_basename)
    root_dir = top_dir.path
    server = GitSandbox(
        'server', bare=True, root_dir=root_dir, create_tmp_dir=False
        )
    git1 = GitSandbox(
        'git1', root_dir=root_dir, create_tmp_dir=False,
        initial_commit=True
        )
    git2 = GitSandbox(
        'git2', root_dir=root_dir, create_tmp_dir=False
        )
    git1('remote', 'add', 'origin', server.directory)
    git1('remote')
    git1('push', '--set-upstream', 'origin', 'master')
    git2('remote', 'add', 'origin', server.directory)
    git2('fetch')
    git2('checkout', 'master')
    return (server, git1, git2)


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
            bare=None, initial_commit=False, root_dir=None,
            create_tmp_dir=True
            ):
        '''Arguments:
        name           -- the name of the repository
        bare           -- if true, create a bare repository
        initial_commit -- if true, add an initial (empty) commit
        root_dir       -- set the directory under which to create the
                          repository
        create_tmp_dir -- if true (default), create a temporary directory
                          based on a prefix + NAME.
                          Otherwise, create a directory called NAME.
        '''
        if create_tmp_dir:
            dir_basename = 'git-pc-test_' + name
        else:
            dir_basename = name
        if bare:
            suffix = '.git'
        else:
            suffix = ''
        if create_tmp_dir:
            if root_dir:
                self.tmp_dir = TmpDir(
                    root_dir, dir_basename, suffix=suffix
                    )
            else:
                self.tmp_dir = TmpDir(None, dir_basename, suffix=suffix)
            self.directory = self.tmp_dir.path
        else:
            if root_dir:
                self.directory = os.path.join(
                    root_dir, dir_basename + suffix
                    )
            else:
                self.directory = os.path.join(
                    tempfile.gettempdir(), dir_basename + suffix
                    )
            os.mkdir(self.directory)
        os.chdir(self.directory)
        if bare:
            self.__call__('init', '--bare')
        else:
            self.__call__('init')
        if initial_commit:
            self.__call__(
                'commit', '--allow-empty',
                '-m', 'Initial commit.'
                )

    def purge(self):
        if self.tmp_dir:
            self.tmp_dir.purge()

    def __call__(self, *args):
        cmd_list = ['git']
        cmd_list.extend(args)
        run_system_cmd(cmd_list, dir=self.directory)

    def commit_all(self, message):
        self.__call__('commit', '-a', '-m', message)

    def system_cmd(self, *args):
        return run_system_cmd_get_output(args, dir=self.directory)

    def write_line_to(self, filename, line):
        '''Create file FILENAME if necessary, and append the given line'''
        path = os.path.join(self.directory, filename)
        f = open(path, 'a')
        f.write(line + '\n')
        f.close()


if __name__ == '__main__':
    main()
