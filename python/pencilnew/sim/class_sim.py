
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
        self.datadir:	  path to simulation datadir (data/)
        self.pendir:	  path to simulation pendir (.pen/)
        self.pendatadir:  path to simulation pendir in datadir (data/.pen/)
        self.param:		  list of all simulation parameters
    """

    def __init__(self, path, hidden=False, quiet=False):
        ## modules
        import os
        from os.path import join
        #from pen.intern.hash_sim import hash_sim
        from pencil import read_param

        ## find out name and store it
        self.name = os.path.split(path)[-1]
        if self.name == '.' or self.name == '':
            self.name = os.path.split(os.getcwd())[-1]

        ## store paths
        self.path = os.path.abspath(path)
        self.dir = self.path
        if (not quiet):
            print '## Creating Simulation object for '+self.path
        self.datadir = join(self.path,'data')		# SIM.datadir
        self.pendir = join(self.path,'.pen')		# SIM.pendir
        self.pendatadir = join(self.path,'data','.pen')	# SIM.pendatadir

        ## generate status hash identification
        #self.status_hash = hash_sim(path)

        ## hidden is default False
        self.hidden = hidden

        ## read params into SIM object
        self.param = {}
        if os.path.exists(join(self.datadir,'param.nml')):
              param = read_param(quiet=True, datadir=self.datadir)
              for key in dir(param):
                    if key.startswith('__'):
                        continue
                    self.param[key] = getattr(param, key)
        else:
              print '?? WARNING: Couldnt find param.nml in simulation '+self.name+'! Simulation is now hidden from calculations!'
              self.param['UNSTARTED'] = True
              self.hidden=True


      ### METHODS ###
      def hide(self):
        self.hidden = True; self.export(); pen.refresh()

      def unhide(self):
        self.hidden = False; self.export(); pen.refresh()

      def export(self):
        """Export simulation object to its root-pendir"""
        from pen.intern.pkl_save import pkl_save
        if self == False:
          print '!! ERROR: Simulation object is bool object and False!'
        pkl_save(self, 'SIM', folder=self.pendir)


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

      def get_datadir(self):
        """Return path to data dir."""
        return self.datadir

      def started(self):
        """Returns whether simulation has already started. This is indicated by existing time_Series.dat in
        data"""
        from os.path import exists
        from os.path import join
        return exists(join(self.path, 'data', 'time_series.dat'))

      def get_varfiles(self, pos=False, particle=False):
        """Get a list of all existing VAR### files.

        pos = False:            give full lust
        pos = last/first:       give newest/first var file
        post = list of numbers: give varfiles at this positions
        particle = True:        return PVAR isntead of VAR list"""
        import glob
        from os.path import join
        from os.path import basename
        from pen.math import natural_sort
        varlist = natural_sort([basename(i) for i in glob.glob(join(self.get_datadir(), 'proc0')+'/VAR*')])
        if particle: varlist = ['P'+i for i in varlist]
        if pos == 'first':
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
          return self.get_varfiles()

      def get_lastvarfilename(self):
        return self.get_varfiles()[-1].split('VAR')[-1]
