import numpy as np

def len_tsfile(datafile):
	#Calculates the length of the time series file
	tsinf = open(datafile, "r")

	line = tsinf.readline() # Skip the title line

	nstep = 0
	for line in tsinf:
		nstep = nstep + 1 

	tsinf.close()

	nstep = nstep - 1 #remove the extra step

	print nstep

	return nstep

def read_tsfile(datafile):
	#READ nstep
	#READ step t dt urms umax umin uxmin uxmax uymin uymax uzmin uzmax Marms Mamax
	nstep = len_tsfile(datafile)

	tsac = open(datafile, "r")

	line = tsac.readline()
	line = line.split()

	headers = line

	step = np.zeros((nstep+1), dtype=int)
	var_list = np.zeros((nstep+1, len(headers)-1))

	kk = 0
	for line in tsac:
		line = line.split()
		step[kk] = int(line[0])
		#print steps[kk]
		for ii in np.arange(len(headers)-1):
			var_list[kk, ii] = float(line[ii+1])
		kk = kk+1 

	tsac.close()

	#Discard the 'step' in headers 
	headers = headers[1:]

	return step, headers, var_list

class tsdata:
	'Includes all the data from the time series. Just needs the name of the datafile.'
	def __init__(self, tsfile):
		self.step, self.headers, self.var_list  = read_tsfile(tsfile)
		#The order of variable in the table: t dt urms uxrms uyrms uzrms rhorms umax rhomax
		for idx, name in enumerate(self.headers, start=0):
			if name == "t":
				self.t = self.var_list[:,idx]
			elif name == "dt":
				self.dt = self.var_list[:,idx]
			elif name == "urms":
				self.urms = self.var_list[:,idx]
			elif name == "uxrms":
				self.uxrms = self.var_list[:,idx]
			elif name == "uyrms":
				self.uyrms = self.var_list[:,idx]
			elif name == "uzrms":
				self.uzrms = self.var_list[:,idx]
			elif name == "uxmax":
				self.uxmax = self.var_list[:,idx]
			elif name == "uymax":
				self.uymax = self.var_list[:,idx]
			elif name == "uzmax":
				self.uzmax = self.var_list[:,idx]
			elif name == "rhorms":
				self.rhorms = self.var_list[:,idx]
			elif name == "umax":
				self.umax = self.var_list[:,idx]
                        elif name == "umin":
                                self.umin = self.var_list[:,idx]
			elif name == "rhomax":
				self.rhomax = self.var_list[:,idx]
                        elif name == "rhomin":
                                self.rhomin = self.var_list[:,idx]
                        elif name == "uxmin":
                                self.uxmin = self.var_list[:,idx]
                        elif name == "uymin":
                                self.uymin = self.var_list[:,idx]
                        elif name == "uzmin":
                                self.uzmin = self.var_list[:,idx]
			else:
				print "Error: " + name + " not recognized is a time series header!"

