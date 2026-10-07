"""
Utility to read the simulation and print the times at which various quantities are output
"""

import pencil as pc
import numpy as np

class Sim(pc.sim.Simulation):
	def __init__(self, *args, quiet=True, **kwargs):
		super().__init__(*args, quiet=quiet, **kwargs)
		
		self.ts = self.get_ts()
		self.av = pc.read.aver(datadir=self.datadir, simdir=self.path, plane_list=['xy', 'z'])
		self.sl = self._read_slices()
		self.snaps = []
		self.p = pc.read.power(file_name='power_sp.dat')
		
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

if __name__ == "__main__":
	sim = Sim(path=".")
	
	output = f"""
Time series:
	it:		{sim.ts.it}
	t:		{sim.ts.t}
	specialm:	{sim.ts.specialm}

Averages:
	t:			{sim.av.t}
	specialmz[:,0]:		{sim.av.xy.specialmz[:,0]}
	specialmxy[:,0,0]:	{sim.av.z.specialmxy[:,0,0]}

Snapshots:
	t:		{[float(var.t) for var in sim.snaps]}
	special[0,0,0]:	{[float(var.special[0,0,0]) for var in sim.snaps]}

Slices:
	t:			{sim.sl.t}
	xy.special[:,0,0]:	{sim.sl.xy.special[:,0,0]}

Power spectra:
	t:		{sim.p.t}
	sqrt(sp[:,0]):	{np.sqrt(sim.p.sp[:,0])}
"""
	print(output)
