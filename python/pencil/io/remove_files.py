def remove_files(path, do_it=False, do_it_really=False):
    """ This method clears path COMPLETELY.
    Meaning, if you start this method my using path on your root dir
    your whole computer would be gone! So we implemented some safety mechanism.

    Args:
        path:                   path do clear
        do_it, do_it_really:    to activate pass True, True
    """

    import os, shutil
    import os
    from os.path import expanduser, isdir, realpath, relpath

    # safety checks
    # first get all sensible paths

    sens_paths = []
    for key in os.environ.keys():
        paths = os.environ[key].split(':')
        for p in paths:
            if isdir(p):
                sens_paths.append(realpath(p))

    sens_paths.append('/')
    sens_paths.append(expanduser('~'))

    for p in sens_paths:
        if not realpath(p) in sens_paths: sens_paths.append(realpath(p))
        while p != '/':
            p = os.path.dirname(p)
            if not realpath(p) in sens_paths: sens_paths.append(p)


    if realpath(path) in sens_paths:
        print('! ERROR: You better should not delete '+path)
        return False


    # Now I really remove files, or show what I would delete
    # This is super slow :(( Anyone an idea how to speed up?
    if os.path.exists(path):
        if do_it and do_it_really:
            os.system('rm -rf '+path)
        else:
            print('?? WARNING: Would remove: '+path)

    return True
