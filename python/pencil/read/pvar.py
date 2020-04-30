def pvar(*args, **kwargs):
    """
    Read PVAR files from Pencil Code using IDL. Does also work with block decomposition.
    Uses IDL<->Python Bridge, this must be activated manually!

    !! WARNING: SHAPE IS AS IN IDL: (X, Y, Z) !!

    Args:
        - varfile       put 'PVARXYZ' or just number here, 'VAR' will be replaced by 'PVAR' autom.
        - npar_max      maximal number of particles to be read in

        - datadir      specify datadir, default False
        - sim           specify simulation from which you want to read
        - proc          read from single proc, set number here
        - swap_endian   change if needed to True, default False
        - quiet         verbosity, default False

    If needed add manually to this script:
        - rmv, irmv, trmv, oldrmv are used for ???
        - solid_object is used for ???
        - theta_arr is used for ???
        - savefile is used for ???

    """

    var_tmp = ParticleData(*args, **kwargs)
    return var_tmp


class ParticleData(object):
    """
    Read PVAR files from Pencil Code using IDL.
    Uses IDL<->Python Bridge, this must be activated manually!

    !! WARNING: SHAPE IS AS IN IDL: (X, Y, Z) !!

    Args:
        - datadir      specify datadir, default False
        - sim           specify simulation from which you want to read
        - varfile       put 'PVARXYZ' or just number here, 'VAR' will be replaced by 'PVAR' autom.
        - npar_max      maximal number of particles to be read in

        - proc          read from single proc, set number here
        - swap_endian   change if needed to True, default False
        - quiet         verbosity, default False

    If needed add manually to this script:
        - rmv, irmv, trmv, oldrmv are used for ???
        - solid_object is used for ???
        - theta_arr is used for ???
        - savefile is used for ???

    """

    def __init__(self, varfile='pvar.dat', npar_max=-1,
                 datadir=False, sim=False, proc=-1, swap_endian=False, quiet=False, DEBUG=False):
        """
        Read PVAR files from Pencil Code using IDL.
        Uses IDL<->Python Bridge, this must be activated manually!

        Args:
            - datadir      specify datadir, default False
            - sim           specify simulation from which you want to read
            - varfile       put 'PVARXYZ' or just number here, 'VAR' will be replaced by 'PVAR' autom.
            - npar_max      maximal number of particles to be read in

            - proc          read from single proc, set number here
            - swap_endian   change if needed to True, default False
            - quiet         verbosity, default False

        """

        import numpy as np
        import os
        from .. import get_sim
        from ..math import is_number
        from sys import byteorder

        ####### interpret parameters
        if datadir == False:
            if sim == False:
                sim = get_sim()
        datadir = sim.datadir

        if os.path.exists(os.path.join(datadir,'grid.h5')):
            l_h5 = True
            import h5py 
        if not l_h5:
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
                return None

        if not l_h5:
            if quiet == False:
                quiet = '0'
            else:
                quiet = '1'

            if swap_endian == False:
                if byteorder == 'little': swap_endian = '0'
                elif byteorder == 'big': swap_endian = '1'
            else: print('? WARNING: Couldnt determine endianness!')

        ####### preparing IDL call
        # cleanup of varfile string
        if is_number(varfile): varfile = 'PVAR'+str(varfile)
        varfile = str(varfile)
        if varfile=='var.dat': varfile='pvar.dat'
        if varfile[:3]=='VAR': varfile='P'+varfile
        if l_h5:
            varfile = str.strip(varfile,'.dat')+'.h5'
            with h5py.File(os.path.join(datadir,'allprocs',varfile),'r') as hf:
                for key in hf['part'].keys():
                    setattr(self, key.lower(), hf['part'][key][()])
        #
        else:
            idl_call = ', '.join(['pc_read_pvar',
                                  'obj=pvar',
                                  'varfile="'+varfile+'"',
                                  'datadir="'+datadir+'"',
                                  'quiet='+quiet,
                                  'swap_endian='+swap_endian,
                                  'proc='+str(proc)
                                  ])

            # reduce number of particles to be read in
            if npar_max > 0: idl_call= idl_call+', npar_max='+str(npar_max)

            ####### show idl_call string if DEBUG
            if DEBUG == True: print('~ DEBUG: idl_call: '+idl_call)

            ###### read in var file in IDL
            print('~ reading '+varfile+' in IDL..')
            IDL.run(idl_call)

            ####### parse to python
            print('~ parsing PVAR from IDL to python..')
            pvar = IDL.pvar

            for key in pvar.keys():
                setattr(self, key.lower(), pvar[key])
            setattr(self, 'xp', pvar['XX'][0])
            setattr(self, 'yp', pvar['XX'][1])
            setattr(self, 'zp', pvar['XX'][2])
            setattr(self, 'vpx', pvar['VV'][0])
            setattr(self, 'vpy', pvar['VV'][1])
            setattr(self, 'vpz', pvar['VV'][2])
