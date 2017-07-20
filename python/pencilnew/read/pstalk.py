def pstalk(*args, **kwargs):
    """
    Read PSTALK files from Pencil Code using IDL.
    Uses IDL<->Python Bridge, this must be activated manually!

    Args:
        - datadir      specify datadir, default False
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

    def __init__(self, datadir=False, sim=False, noutmax='-1',
                 swap_endian=False, quiet=False):
        """
        Read PSTALK files from Pencil Code using IDL.
        Uses IDL<->Python Bridge, this must be activated manually!

        Args:
            - datadir      specify datadir, default False
            - sim           specify simulation from which you want to read
            - swap_endian   change if needed to True, default False
            - quiet         verbosity, default False


        """

        import numpy as np
        import os
        from os.path import join
        import pencilnew as pcn

        if datadir == False:
            if sim == False:
                sim = pcn.get_sim()
            datadir = sim.datadir

        if quiet == False:
            quiet = '0'
        else:
            quiet = '1'

        if swap_endian == False:
            swap_endian = '0'
        else:
            swap_endian = '1'

        try:
            cwd = os.getcwd()
            from idlpy import IDL
            os.chdir(cwd)

            print('~ reading pstalk in IDL..')

            idl_call = ', '.join(['pc_read_pstalk', 'obj=pstalk', 'datadir="'+datadir+'"', 'quiet='+quiet, 'swap_endian='+swap_endian, 'noutmax='+str(noutmax)])

            IDL.run(idl_call)

            print('~ parsing pstalk from IDL to python..')
            ps = IDL.pstalk

            for key in ps.keys():
                setattr(self, key, ps[key].T)

        except:
            print('! ERROR: no idl<->python bridge found. Try whats written in pstalk-comment to fix that issue.')
            print('! ')
            print('! Use something like: (enshure you have IDL 8.5.1 or larger)')
            print('! export PYTHONPATH=$PYTHONPATH:$IDL_HOME/lib/bridges:$IDL_HOME/bin/bin.linux.x86_64')
            print('! export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64:$IDL_HOME/bin/bin.linux.x86_64')
            print('! in your .bashrc')
            print('! ')
            print('! If you have it already installed, try: from idlpy import IDL and check for errors')
            print('! ')
            print('~ BACKUP SOLUTION: reading pstalk via pidly, starting IDL..')

            from pencilnew.backpack import pidly
            IDL = pidly.IDL(long_delay=0.05)	# start IDL engine
            from scipy.io.idl import readsav
            from pencilnew.io import mkdir

            ## read tstalk file
            print('## reading particle stalker file..')
            IDL('pc_read_pstalk, object=pstalk, datadir="'+sim.datadir+'"'+', quiet='+quiet+', swap_endian='+swap_endian)

            print('## transfering pstalk file from IDL to python..')
            mkdir(join(sim.pc_datadir,'tmp'))
            IDL('save, pstalk, filename="'+join(sim.pc_datadir,'tmp','pstalk.sav')+'"')
            ps = readsav(join(sim.pc_datadir,'tmp','pstalk.sav'))('pstalk')

            #from pencilnew.io import debug_breakpoint; debug_breakpoint()

            for key in set(ps.dtype.fields.keys()):
                key = key.lower()
                setattr(self, key, ps[key][0].T)
