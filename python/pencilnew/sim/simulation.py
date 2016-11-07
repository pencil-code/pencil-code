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

    Properties:

        self.name:          name of
        self.path:          path to simulation
        self.data_dir:      path to simulation data-dir (./data/)
        self.pc_dir:        path to simulation pc-dir (./.pc/)
        self.pc_data_dir:   path to simulation pendir in data_dir (data/.pc/)
        self.components:    list of files which are components of the specific simulation
        self.hidden:        Default is False, if True this simulation will be ignored by pencil (but still is part of simdict)
        #self.status_hash:   Hash key representing the state of a simulation, generated from all components
        self.param:         list of param file
        self.grid:          grid object
    """

    return __Simulation__(*args, **kwargs)

class __Simulation__(object):
    """
    Simulation object
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

        #self.status_hash = hash_sim(path   # generate status hash identification
        self.hidden = hidden                # hidden is default False
        self.update(self)                   # auto-update, i.e. read param.nml


        # Done

    def update(self, quiet=True):
        from os.path import exists
        from os.path import join
        from pencilnew.read import param

        self.param = {}                     # read params into Simulation object
        # try:
        if exists(join(self.data_dir,'param.nml')):
            param = param(quiet=quiet, data_dir=self.data_dir)
            # from pencil as import read_param
            # param = read_param(quiet=True, datadir=self.data_dir)
            for key in dir(param):
                if key.startswith('__'): continue
                self.param[key] = getattr(param, key)
        else:
            print('? WARNING: Couldnt find param.nml for '+self.path)
        # except:
        #     print('! ERROR: while reading param.nml for '+self.path)

        self.grid = None
        self.ghost_grid = None
        try:                                # read grid
            self.grid = pencilnew.read.grid(data_dir=self.data_dir, trim=True, quiet=True)
            self.ghost_grid = pencilnew.read.grid(data_dir=self.data_dir, trim=False, quiet=True)
        except:
            if (not quiet): print('? WARNING: Couldnt load grid for '+self.path)

        self.export()

    def hide(self):
        self.hidden = True; self.export()

    def unhide(self):
        self.hidden = False; self.export()

    def export(self):
        """Export simulation object to its root/.pc-dir"""
        from pencilnew.io import save
        if self == False: print('! ERROR: Simulation object is bool object and False!')
        save(self, name='sim', folder=self.pc_dir)

    #   def unchanged(self):
    #     """Check if simulation hash has changed."""
    #     from pen.intern.hash_sim import hash_sim
      #
    #     old_hash = self.status_hash
    #     new_hash = hash_sim(self.path)
      #
    #     if old_hash == new_hash:
    #       return True
    #     else:
    #       return False

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

    def get_value_from_file(self, filename, quantity):
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
            filepath = join(self.dir, filename)
        elif filename in ['cparam.local']:
            filepath = join(self.dir, 'src', filename)
        else:
            print('!! Quantity '+quantity+' not found in '+filepath); return False


        with open(filepath, 'r') as f: data = f.readlines()     # open file and read content

        for id,line in enumerate(data):     # check lines for quantity
            if line.find(quantity) >= 0:
                # special cases:
                if quantity == 'Lxyz':
                    line = line.replace(' ', '')
                    line = line.split(quantity+'=')[-1]
                    return [float(i) for i in line.split(',')]

                # default case:
                for i in line.split(quantity+'=')[-1]:
                    if is_number(i) == False and not i in ['.','e', '-', '+']: break
                oldvalue = line.split(quantity+'=')[-1].split(i)[0]     # extract quantity value (not needed here)

                if is_int(oldvalue): return int(float(oldvalue))
                elif is_float(oldvalue): return float(oldvalue)
                else:
                    print('?? WARNING: value from file is neither float nor int! value is: '+oldvalue)
                    return oldvalue


        else:
            print('!! quantity '+quantity+' not found in '+filepath); return False
