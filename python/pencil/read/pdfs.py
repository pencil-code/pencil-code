"""
Reads the output of the pdf subroutine (power_spectrum.f90)

03-Nov-2025/Kishore: coded
"""

import numpy as np
import pathlib
import warnings

from pencil import read

class pdf:
	"""
	Reads the output of the power_xy subroutine
	"""
	def __init__(self, datadir="data", filename=None, debug=False):
		"""
		
		Arguments
			datadir: path to the data directory. Default: "data"
			filename: name of the file to read. Default: all files named pdf_*.dat
		
		Examples
			>>> p = pc.read.pdfs()
			>>> p.t #1D array, length nt
			>>> p.lnrho #PDFContainer instance
			>>> p.lnrho.bins: #1D array, length nb
			>>> p.lnrho.counts #2D array, size (nt,nb): counts in each bin at a particular timestep
			>>> plt.semilogy(p.lnrho.bin_centers, p.lnrho.calc_counts(t_min=10)) #plots the PDF of lnrho for t>10
		
		Attributes
			t: 1D array, times
			Remaining attributes are PDFContainer instances, named after the variables whose PDFs are being calculated.
		"""
		self._debug = debug
		
		datadir = pathlib.Path(datadir)
		
		if filename is None:
			files = self._get_pdf_files(datadir)
		else:
			files = [datadir/filename]
		
		param = read.param(datadir=datadir)
		
		for f in files:
			self._read_pdf_file(f, param)
	
	def keys(self):
		ks = []
		for k in self.__dict__.keys():
			if k[0] != '_':
				ks.append(k)
		return ks
	
	def _varname_from_file(self, f):
		filename = f.name
		filename = filename.removeprefix("pdf_")
		filename = filename.removesuffix(".dat")
		return filename
	
	def _get_pdf_files(self, datadir):
		files = list(datadir.glob("pdf_*.dat"))
		
		n_files = len(files)
		if n_files == 0:
			warnings.warn(f"no files found in {datadir.path}")
		else:
			if self._debug:
				print(f"Found {n_files} file(s)")
		
		return files
	
	def _is_logscale(self, var):
		"""
		For some reason, the output does not contain any information on whether logscale is T or F. We thus have this hack.
		"""
		if var in ["lncc", "lngcc"]:
			return True
		else:
			return False
	
	def _read_pdf_file(self, filename, param):
		if self._debug:
			print(f"Reading {filename}")
		
		var = self._varname_from_file(filename)
		logscale = self._is_logscale(var)
		
		ts = []
		counts = []
		
		with open(filename, 'r') as f:
			first_line = f.readline()
			metadata = [float(i) for i in first_line.split()]
			
			if logscale:
				_, n_pdf, pdf_dx, pdf_max_logscale, pdf_min_logscale, pdf_mean, pdf_rms = metadata
			else:
				_, n_pdf, pdf_dx, pdf_max, pdf_min, pdf_mean, pdf_rms = metadata
		
		n_till_now = n_pdf
		with open(filename, 'r') as f:
			for line in f:
				if n_till_now == n_pdf:
					t, *_ = [float(i) for i in line.split()]
					
					n_till_now = 0
					ts.append(t)
					counts.append([])
				elif n_till_now > n_pdf:
					raise RuntimeError(f"too many values obtained (got {n_till_now}, expected {n_pdf})")
				else:
					data = [float(i) for i in line.split()]
					n_till_now += len(data)
					counts[-1].extend(data)
		
		counts = np.array(counts)
		
		#Sanity checks
		assert counts.ndim == 2
		assert len(ts) == counts.shape[0]
		
		#calculate the bin edges
		#Note that the int intrinsic in Fortran rounds down (at least in gfortran)
		#TODO: the left and right bin edges should really be -np.inf and np.inf, but doing that would make the PDF zero in the first and last bins.
		if logscale:
			x = pdf_mean + pdf_rms * (pdf_min_logscale + 10**(pdf_dx*np.arange(n_pdf+1)))
		else:
			x = pdf_mean + pdf_rms*(pdf_min + pdf_dx*np.arange(n_pdf+1))
		
		self.t = np.array(ts)
		setattr(self, var, PDFContainer(bins=x, ts=self.t, counts=counts, logscale=logscale))

class PDFContainer:
	"""
	Attributes
		bins: 1D array, bin edges
		counts: 2D array, counts in bins (second axis) at different times (first axis)
	"""
	def __init__(self, bins, ts, counts, logscale):
		if not ((counts.ndim == 2) and (counts.shape == (len(ts), len(bins)-1))):
			raise RuntimeError(f"wrong shape for counts array: expected {(len(ts), len(bins)-1)}, got {counts.shape}")
		
		self.bins = bins
		self.counts = counts
		self._t = ts
		self._logscale = logscale
	
	@property
	def bin_centers(self):
		if self._logscale:
			return np.sqrt(self.bins[:-1] * self.bins[1:])
		else:
			return (self.bins[:-1] + self.bins[1:])/2
	
	def calc_counts(self, t_min=None, t_max=None):
		if t_min is None:
			it_min = 0
		else:
			it_min = np.argmin(np.abs(self._t - t_min))
		
		if t_max is None:
			it_max = len(self._t)
		else:
			it_max = np.argmin(np.abs(self._t - t_max))
		
		return np.sum(self.counts[it_min:it_max+1,:], axis=0)
	
	def calc_pdf(self, *args, **kwargs):
		"""
		Calculate the PDF of the variable 'variable' by summing over the counts saved between t_min and t_max.
		"""
		counts_allt = self.calc_counts(*args, **kwargs)
		
		dx = self.bins[1:] - self.bins[:-1]
		pdf = counts_allt/(dx*np.sum(counts_allt))
		return pdf
