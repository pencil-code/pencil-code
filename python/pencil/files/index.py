#$Id$
import os
from pencil import read_param, read_dim


def read_index(*args, **kwargs):
	"""Read index.pro and return an Index object.

	Params:
	------
	datadir='data/' (optionnal)
	"""
	return Index(*args, **kwargs)


class Index(dict):
	"""
	read index.pro and returns a read_index class which is composed
	of a python dictionnary

	Author: T. Gastine (tgastine@ast.obs-mip.fr)
	"""
	def __init__(self, datadir='data/', param=None, dim=None):
		"""Constructor:
		-----------

		Params:
		------
			datadir='data/' (optionnal)

		Returns:
        -------
		a read_index class
		"""
		if param is None:
			param = read_param(datadir=datadir, quiet=True)
		if dim is None:
			dim = read_dim(datadir=datadir)

		if param.lwrite_aux:
			totalvars = dim.mvar+dim.maux
		else:
			totalvars = dim.mvar

		f = open(os.path.join(datadir,'index.pro'))
		for line in f.readlines():
			clean = line.strip()
			name=clean.split('=')[0].strip().replace('[','').replace(']','')
			if (clean.split('=')[1].strip().startswith('intarr(')): continue
			if (clean.split('=')[1].strip().startswith('indgen(')):
				stri=clean.split('=')[1].strip()
				stri=stri.replace('(','( ')
				stri=stri.replace(')',' )')
				values = [int(s) for s in stri.split() if s.isdigit()]
				numbers= range(values[0],values[0]+values[1])
				for i in numbers:
					name=clean.split('=')[0].strip().replace('[','').replace(']','')
					name=name.lstrip('i')+str(i-values[1])
					name=name.lstrip('i')
					self[name]=i
					#~ print name, i
				continue
			else:
				val=int(clean.split('=')[1].strip())
           # val=int(float(clean.split('=')[1].strip())) # sean150513

            #            print name,val
            # need to compare val to totalvars as global indices
            # may be present in index.pro
            #            if (val != 0 and val <= totalvars and \
			if (val != 0  and val <= totalvars \
				and not name.startswith('i_') and name.startswith('i')):
				name=name.lstrip('i')
				if (name == 'lnTT' and param.ltemperature_nolog):
					name = 'tt'
			self[name] = val
