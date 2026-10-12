"""
Common functions for all the test scripts in this folder
"""

import pencil as pc

class Sim(pc.sim.Simulation):
    def __init__(self, *args, quiet=True, **kwargs):
        super().__init__(*args, quiet=quiet, **kwargs)
        
        self.ts = self.get_ts()
        self.av = pc.read.aver(datadir=self.datadir, simdir=self.path, plane_list=['xy', 'z'])
        self.sl = self._read_slices()
        self.snaps = []
        self.p = pc.read.power(datadir=self.datadir, file_name='power_sp.dat')
        
        for fname in self.get_varlist():
            self.snaps.append(pc.read.var(var_file=fname, datadir=self.datadir, trimall=True, quiet=True))
        
        self.snaps.sort(key=lambda var: var.t)
    
    def _read_slices(self):
        io = self.param['io_strategy']
        if io == "dist":
            self.bash(
                "pc_build -t read_all_videofiles",
                bashrc=False,
                verbose=False,
                raise_errors=True,
                )
            self.bash(
                "src/read_all_videofiles.x",
                bashrc=False,
                verbose=False,
                raise_errors=True,
                )
        elif io == "HDF5":
            pass
        else:
            raise NotImplementedError(f"unsure whether I can directly read slices for io_strategy={io}")
        
        return pc.read.slices(datadir=self.datadir)
