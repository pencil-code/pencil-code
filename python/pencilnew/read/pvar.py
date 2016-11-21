def pvar(*args, **kwargs):
    """
    Read PVAR files from Pencil Code using IDL. Does also work with block decomposition.
    Uses IDL<->Python Bridge, this must be activated manually!

    Args:
        - varfile       put 'PVARXYZ' or just number here, 'VAR' will be replaced by 'PVAR' autom.
        - npar_max      maximal number of particles to be read in

        - data_dir      specify data_dir, default False
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

    Args:
        - data_dir      specify data_dir, default False
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
                 data_dir=False, sim=False, proc=-1, swap_endian=False, quiet=False, DEBUG=False):
        """
        Read PVAR files from Pencil Code using IDL.
        Uses IDL<->Python Bridge, this must be activated manually!

        Args:
            - data_dir      specify data_dir, default False
            - sim           specify simulation from which you want to read
            - varfile       put 'PVARXYZ' or just number here, 'VAR' will be replaced by 'PVAR' autom.
            - npar_max      maximal number of particles to be read in

            - proc          read from single proc, set number here
            - swap_endian   change if needed to True, default False
            - quiet         verbosity, default False

        """

        import numpy as np
        import os
        import pencilnew as pcn
        from pencilnew.math import is_number

        try:
            cwd = os.getcwd()
            from idlpy import IDL
            os.chdir(cwd)

        except:
            print('! ERROR: no idl<->python bridge found. Try whats written in pstalk-comment to fix that issue.')
            print('! ')
            print('! Use something like: ')
            print('! export PYTHONPATH=$PYTHONPATH:$IDL_HOME/lib/bridges:$IDL_HOME/bin/bin.linux.x86_64')
            print('! export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64:$IDL_HOME/bin/bin.linux.x86_64')
            print('! in your .bashrc')
            print('! ')
            break

        ####### interprate parameters
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

        ####### preparing IDL call
        # cleanup of varfile string
        if is_number(varfile): varfile = 'PVAR'+str(varfile)
        varfile = str(varfile)
        if varfile=='var.dat': varfile='pvar.dat'
        if varfile[:3]=='VAR': varfile='P'+varfile

        #
        idl_call = ', '.join(['pc_read_pvar',
                              'obj=pvar',
                              'varfile="'+varfile+'"',
                              'datadir="'+data_dir+'"',
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
            setattr(self, key, pvar[key])
