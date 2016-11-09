# simulation.py
#
# Create simulation object to operate on.
#
# Authors:
# A. Schreiber (aschreiber@mpia.de)
"""
Contains the simulation class which can be used to directly create, access and manipulate simulations.
"""


def simulation(*args, **kwargs):
    """
    Generate simulation object from parameters.
    Simulation objects are containers for simulations. pencil can work with several of them at once if stored in a list or dictionary.

    Args for Constructor:
        path:		path to simulation, default = '.'
        hidden:     set True to set hidden flag, default is False
        quiet:      suppress irrelevant output, default False

    Properties:

        self.name:          name of
        self.path:          path to simulation
        self.data_dir:      path to simulation data-dir (./data/)
        self.pc_dir:        path to simulation pc-dir (./.pc/)
        self.pc_data_dir:   path to simulation pendir in data_dir (data/.pc/)
        self.components:    list of files which are nessecarry components of the simulation
        self.optionals:     list of files which are optional components of the simulation
        self.hidden:        Default is False, if True this simulation will be ignored by pencil (but still is part of simdict)
        self.param:         list of param file
        self.grid:          grid object
    """

    return __Simulation__(*args, **kwargs)

class __Simulation__(object):
    """
    Simulation object.
    """

    def __init__(self, path='.', hidden=False, quiet=False):
        import os
        from os.path import join
        from os.path import exists
        from os.path import split
        #from pen.intern.hash_sim import hash_sim

        self.name = split(path)[-1]     # find out name and store it
        if self.name == '.' or self.name == '': self.name = split(os.getcwd())[-1]

        self.path = os.path.abspath(path)   # store paths
        if (not quiet): print('# Creating Simulation object for '+self.path)
        self.data_dir = join(self.path,'data')
        self.pc_dir = join(self.path,'.pc')
        self.pc_data_dir = join(self.path,'data','.pc')

        self.components = ['src/cparam.local', 'src/Makefile.local', 'start.in', 'run.in', 'print.in'] # core files of a simulation run
        self.optionals = ['*.in', '*.py', 'submit*']    # optinal files that should stick with the simulation when copied

        self.hidden = hidden                # hidden is default False
        self.update(self)                   # auto-update, i.e. read param.nml
        # Done


    def update(self, quiet=True):
        """Update simulation object:
                - read param.nml
                - read grid and ghost grid
        """
        from os.path import exists
        from os.path import join
        from pencilnew.read import param

        self.param = {}                     # read params into Simulation object
        try:
            if exists(join(self.data_dir,'param.nml')):
                param = param(quiet=quiet, data_dir=self.data_dir)
                # from pencil as import read_param
                # param = read_param(quiet=True, datadir=self.data_dir)
                for key in dir(param):
                    if key.startswith('__'): continue
                    self.param[key] = getattr(param, key)
            else:
                print('? WARNING: Couldnt find param.nml for '+self.path)
        except:
            print('! ERROR: while reading param.nml for '+self.path)

        self.grid = None
        self.ghost_grid = None
        try:                                # read grid
            self.grid = pencilnew.read.grid(data_dir=self.data_dir, trim=True, quiet=True)
            self.ghost_grid = pencilnew.read.grid(data_dir=self.data_dir, trim=False, quiet=True)
        except:
            if (not quiet): print('? WARNING: Couldnt load grid for '+self.path)

        self.export()


    def hide(self):
        """Set hide flag True for this simulation. """
        self.hidden = True; self.export()


    def unhide(self):
        """Set hide flag False for this simulation. """
        self.hidden = False; self.export()


    def export(self):
        """Export simulation object to its root/.pc-dir"""
        from pencilnew.io import save
        if self == False: print('! ERROR: Simulation object is bool object and False!')
        save(self, name='sim', folder=self.pc_dir)


    def started(self):
        """Returns whether simulation has already started.
        This is indicated by existing time_series.dat in data directory."""
        from os.path import exists
        from os.path import join
        return exists(join(self.path, 'data', 'time_series.dat'))


    def get_varlist(self, pos=False, particle=False):
        """Get a list of all existing VAR# file names.

        pos = False:                 give full list
        pos = 'last'/'first':        give latest/first var file
        pos = 'lastXXX' / 'firstXXX' give last/first XXX varfiles
        pos = list of numbers:       give varfiles at this positions
        particle = True:             return PVAR- instead of VAR-list"""

        import glob
        from os.path import join as join
        from os.path import basename
        from pencilnew.math import natural_sort

        key = 'VAR'
        if particle == True: key = 'PVAR'

        varlist = natural_sort([basename(i) for i in glob.glob(join(self.data_dir, 'proc0')+'/'+key+'*')])
        #if particle: varlist = ['P'+i for i in varlist]

        if pos == False: return varlist
        if pos == 'first': return [varlist[0]]
        if pos == 'last': return [varlist[-1]]
        if pos.startswith('last'): return varlist[-int(pos[4:]):]
        if pos.startswith('first'): return varlist[:int(pos[4:])]
        if type(pos) == type([]):
            if pos[0].startswith('VAR'): pos = [i[3:] for i in pos]
            if pos[0].startswith('PVAR'): pos = [i[3:] for i in pos]
            return [varlist[int(i)] for i in pos]
        return varlist


    def get_pvarlist(self, pos=False):
        """Same as get_varfiles(pos, particles=True). """
        return self.get_varfiles(pos=pos, particle=True)


    def get_lastvarfilename(self, particle=False):
        """Returns las varfile name as string."""
        return self.get_varfiles(pos='last', particle=particle)


    def get_value_from_file(self, filename, quantity, DEBUG=False):
        """ Use to read in a quantity from *.in or *.local files.

        Args:
            file:       can be "run.in", "start.in", "cparam.local"
            quantity:   variable to read in from file
        """
        from os.path import join
        from pencilnew.math import is_number
        from pencilnew.math import is_float
        from pencilnew.math import is_int

        if filename in ['run.in', 'start.in']:
            filepath = join(self.path, filename)
        elif filename in ['cparam.local']:
            filepath = join(self.path, 'src', filename)
        else:
            print('! Filename '+filename+' could not be interprated for '+filepath); return False

        with open(filepath, 'r') as f: data = f.readlines()     # open file and read content

        for id, line in enumerate(data):        # check lines for quantity
            if line.find(quantity) >= 0:
                # cleanup
                line = line.split('!')[0].strip()       # remove comments and empty spaces
                if not quantity in line:
                    print('? WARNING: Quantity "'+quantity+'" was found in a comment in'+filepath+' -> no action perfomed!')

                # special case double string tuple:
                if quantity in ['idiff', 'ivisc']:
                    line = line.replace(' ', '')
                    line = line.split(quantity+'=')[-1]
                    return [i.replace("'", '').replace('"', '') for i in line.split(',')]

                # special case triple float tuple:
                if quantity in ['Lxyz', 'xyz0', 'beta_glnrho_global']:
                    line = line.replace(' ', '')
                    line = line.split(quantity+'=')[-1]
                    return [float(i) for i in line.split(',')]

                # prepare default cases
                value_str = line.split(quantity+'=')[-1].split('=')[0]

                # case value is string
                if value_str[0] in ['"', "'"]: return value_str.replace("'", '').replace('"', '')

                # case value is bool
                if value_str[0] == 'T': return 'T'
                if value_str[0] == 'F': return 'F'

                # case value is number
                if is_number(value_str[0]):
                    print("NUMBER")
                    for i in value_str:      # iterate through string for extract number
                        if is_number(i) == False and not i in ['.','e', '-', '+']: break
                    value_str = value_str.split(i)[0]     # extract quantity value

                if is_int(value_str): return int(float(value_str))
                elif is_float(value_str): return float(value_str)
                else:
                    print('? WARNING: value from file is neither float nor int! value is: '+value_str)
                    print('? Check manually! ')
                    return

        else:
            print('! ERROR: Quantity "'+quantity+'" was not found in '+filepath)
            return False
