def pstalk(*args, **kwargs):
    """
    Read PSTALK files from Pencil Code using IDL.
    Uses IDL<->Python Bridge, this must be activated manually!

    Args:
        - data_dir      specify data_dir, default False
        - sim           specify simulation from which you want to read
        - swap_endian   change if needed to True, default False
        - quiet         verbosity, default False
    """

    var_tmp = ParticleStalkData(*args, **kwargs)
    return var_tmp


class ParticleStalkData(object):
    """
    ParticleStalkData -- holds Pencil Code PSTALK file data.
    """

    def __init__(self, data_dir=False, sim=False,
                 swap_endian=False, quiet=False):
        """
        Read PSTALK files from Pencil Code using IDL.
        Uses IDL<->Python Bridge, this must be activated manually!

        Args:
            - data_dir      specify data_dir, default False
            - sim           specify simulation from which you want to read
            - swap_endian   change if needed to True, default False
            - quiet         verbosity, default False


        """

        import numpy as np
        import os
        import pencilnew as pcn
        try:
            cwd = os.getcwd()
            from idlpy import IDL
            os.chdir(cwd)

        except:
            print('! ERROR: no idl<->python bridge found. Try whats written in pstalk-comment to fix that issue.')
            print('! ')
            print('! Use something like: (enshure you have IDL 8.5.1 or larger)')
            print('! export PYTHONPATH=$PYTHONPATH:$IDL_HOME/lib/bridges:$IDL_HOME/bin/bin.linux.x86_64')
            print('! export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64:$IDL_HOME/bin/bin.linux.x86_64')
            print('! in your .bashrc')
            print('! ')

        print('~ reading pstalk in IDL..')

        if data_dir == False:
            if sim == False:
                sim = pcn.get_sim()
            data_dir = sim.data_dir

        if quiet == False:
            quiet = '0'
        else:
            quiet = '1'

        if swap_endian == False:
            swap_endian = '0'
        else:
            swap_endian = '1'

        idl_call = ', '.join(['pc_read_pstalk', 'obj=pstalk', 'datadir="'+data_dir+'"', 'quiet='+quiet, 'swap_endian='+swap_endian])

        IDL.run(idl_call)

        print('~ parsing pstalk from IDL to python..')
        ps = IDL.pstalk

        for key in ps.keys():
            setattr(self, key, ps[key].T)
