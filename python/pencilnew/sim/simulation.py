
class Simulation:
    """Simulation objects are containers for simulations. pencil can work with several of them at once if stored in a list.

    Args for Constructor:
    path:		path to simulation, default = '.'

    Properties:
        self.name:		  name of simulation
        self.hidden:	  Default is False, if True this simulation will be ignored by pencil (but still is part of
        simdict)
        self.status_hash: Hash key representing the state of a simulation
        self.path:		  path to simulation
        self.dir          same as self.path
        self.data_dir:	  path to simulation data-dir (./data/)
        self.pc_dir:	      path to simulation pc-dir (./.pc/)
        self.pc_data_dir:  path to simulation pendir in data_dir (data/.pc/)
        self.param:		  list of param file
        self.grid:        grid object
    """

    def __init__(self, path, hidden=False, quiet=False):
        import os
        from os.path import join as __join__
        from os.path import exists as __exists__
        from os.path import split as __split__
        #from pen.intern.hash_sim import hash_sim
        from pencilnew.read import param as read_param

        # find out name and store it
        self.name = __split__(path)[-1]
        if self.name == '.' or self.name == '':
            self.name = __split__(os.getcwd())[-1]

        # store paths
        self.path = os.path.abspath(path)
        self.dir = self.path
        if (not quiet): print '# Creating Simulation object for '+self.path
        self.data_dir = __join__(self.path,'data')
        self.pc_dir = __join__(self.path,'.pc')
        self.pc_data_dir = __join__(self.path,'data','.pc')

        # generate status hash identification
        #self.status_hash = hash_sim(path)

        # hidden is default False
        self.hidden = hidden

        # read params into SIM object
        self.param = {}
        if __exists__(__join__(self.data_dir,'param.nml')):
              param = read_param(quiet=True, data_dir=self.data_dir)
              for key in dir(param):
                    if key.startswith('__'): continue
                    self.param[key] = getattr(param, key)
        else:
              print '?? WARNING: Couldnt find param.nml in simulation '+self.name+'! Simulation is now hidden from calculations!'
              self.param['UNSTARTED'] = True
              self.hidden=True

        try:
            self.grid = pencilnew.read.grid(data_dir=self.data_dir, trim=True, quiet=True)
            self.ghost_grid = pencilnew.read.grid(data_dir=self.data_dir, trim=False, quiet=True)
        except:
            self.grid = None
            self.ghost_grid = None


    def hide(self):
        self.hidden = True; self.export(); pen.refresh()

    def unhide(self):
        self.hidden = False; self.export(); pen.refresh()

    def export(self):
        """Export simulation object to its root-pendir"""
        from pencilnew.io import save
        if self == False:
            print '!! ERROR: Simulation object is bool object and False!'
        save(self, 'sim', folder=self.pc_dir)


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
        """Returns whether simulation has already started. This is indicated by existing time_Series.dat in data"""
        return ___exists___(___join___(self.path, 'data', 'time_series.dat'))

    def get_varfiles(self, pos=False, particle=False):
        """Get a list of all existing VAR# file names.

        pos = False:            give full list
        pos = 'last'/'first':       give newest/first var file
        post = list of numbers: give varfiles at this positions
        particle = True:        return PVAR isntead of VAR list"""
        import glob
        from os.path import join as __join__
        from os.path import basename
        from pencilnew.math import natural_sort
        varlist = natural_sort([basename(i) for i in glob.glob(__join__(self.data_dir, 'proc0')+'/VAR*')])
        if particle: varlist = ['P'+i for i in varlist]
        if pos == False:
            return varlist
        elif pos == 'first':
            return [varlist[0]]
        elif pos == 'last':
            return [varlist[-1]]
        elif pos.startswith('last'):
            return varlist[-int(pos[4:]):]
        elif pos.startswith('first'):
            return varlist[:int(pos[4:])]
        elif type(pos) == type([]):
            if pos[0].startswith('VAR'): pos = [i[3:] for i in pos]
            if pos[0].startswith('PVAR'): pos = [i[3:] for i in pos]
            return [varlist[i] for i in pos]
        else:
            return varlist

    def get_varlist(self, pos=False, particle=False):
        """Same as get_varfiles. """
        return self.get_varfiles()

    def get_lastvarfilename(self):
        """Returns las varfile name as string."""
        return self.get_varfiles()[-1].split('VAR')[-1]
