def is_sim_dir(path='.'):
    """Checks if a path is pointing at a pencil code simulation directory. Therefor, it checks the existance of start.in, run.in, src/cparam.local and src/Makefile.local"""

    from os.path import isdir as __isdir__
    from os.path import join as __join__
    from os.path import exists as __exists__

    if __exists__(__join__(path, 'run.in')) and __exists__(__join__(path, 'start.in')) and __exists__(__join__(path, 'src/cparam.local')) and __exists__(__join__(path, 'src/Makefile.local')):
        return True
    return False
