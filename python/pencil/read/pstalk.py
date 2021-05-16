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

    def __init__(self, datadir=False, sim=False,
                 tmin=0, tmax=-1, noutmax='-1',
                 swap_endian=False, quiet=False, use_existing_pstalk_sav=False):
        """
        Read PSTALK files from Pencil Code using IDL.
        Uses IDL<->Python Bridge, this must be activated manually!

        Args:
            - datadir      specify datadir, default False
            - sim           specify simulation from which you want to read
            - swap_endian   change if needed to True, default False
            - quiet         verbosity, default False
            - use_existing_pstalk_sav
                            use existing <sim.datadir>/data/pc/tmp/pstalk.sav for speed up


        """

        import numpy as np
        import os
        from os.path import join
        from .. import get_sim

        if datadir == False:
            if sim == False:
                sim = get_sim()
            datadir = sim.datadir

        if quiet == False:
            quiet = '0'
        else:
            quiet = '1'

        if swap_endian == False:
            swap_endian = '0'
        else:
            swap_endian = '1'

        if use_existing_pstalk_sav == True:
            from scipy.io.idl import readsav

            print('~ reading existing pstalk..')

            ps = readsav(join(sim.pc_datadir,'tmp','pstalk.sav'))('pstalk')

            for key in set(ps.dtype.fields.keys()):
                if hasattr(self, key.lower()): continue
                setattr(self, key.lower(), ps[key][0].T)

        else:
            try:
                cwd = os.getcwd()
                from idlpy import IDL
                os.chdir(cwd)

                print('~ reading pstalk in IDL..')

                idl_call = ', '.join(['pc_read_pstalk', 'obj=pstalk', 'datadir="'+datadir+'"', 'it0='+str(tmin), 'it1='+str(tmax), 'quiet='+quiet, 'swap_endian='+swap_endian, 'noutmax='+str(noutmax)])

                IDL.run(idl_call)

                print('~ parsing pstalk from IDL to python..')
                ps = IDL.pstalk

                for key in set(ps.keys()):
                    if hasattr(self, key.lower()): continue
                    setattr(self, key.lower(), ps[key].T)

            except:
                print('! ERROR: no idl<->python bridge found. Try whats written in pstalk-comment to fix that issue.')
                print('! ')
                print('! Use something like: (ensure you have IDL 8.5.1 or larger)')
                print('! export PYTHONPATH=$PYTHONPATH:$IDL_HOME/lib/bridges:$IDL_HOME/bin/bin.linux.x86_64')
                print('! export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64:$IDL_HOME/bin/bin.linux.x86_64')
                print('! in your .bashrc')
                print('! ')
                print('! If you have it already installed, try: from idlpy import IDL and check for errors')
                print('! ')
                print('~ BACKUP SOLUTION: reading pstalk via pidly, starting IDL..')

                from ..backpack import pidly
                IDL = pidly.IDL(long_delay=0.05)	# start IDL engine
                from scipy.io.idl import readsav
                from ..pio import mkdir

                ## read tstalk file
                print('## reading particle stalker file..')
                IDL('pc_read_pstalk, object=pstalk, datadir="'+sim.datadir+'"'+', quiet='+quiet+', it0='+str(tmin)+', it1='+str(tmax)+', swap_endian='+swap_endian)

                print('## transfering pstalk file from IDL to python..')
                mkdir(join(sim.pc_datadir,'tmp'))
                IDL('save, pstalk, filename="'+join(sim.pc_datadir,'tmp','pstalk_'+str(tmin)+'_'+str(tmax)+'.sav')+'"')
                ps = readsav(join(sim.pc_datadir,'tmp','pstalk.sav'))('pstalk')

                #from pc.io import debug_breakpoint; debug_breakpoint()

                for key in set(ps.dtype.fields.keys()):
                    if hasattr(self, key.lower()): continue
                    setattr(self, key.lower(), ps[key][0].T)
