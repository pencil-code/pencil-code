# --------------------------------------------------------------------------
# Author: Daniel Carrera (danielc@astro.lu.se)
# --------------------------------------------------------------------------
# Read index.pro and return it as a dictionary.
# 
# Usage:   read_index() -- Read data/index.pro
# --------------------------------------------------------------------------

function read_index(;datadir="data")
	
	index = Dict()
	lines = readlines( open(datadir * "/index.pro","r") )
	
	for line = lines
	    key, val = split(line, r" *= *")
	    index[ strip(key) ] = int(val)
	end
	
	k = split(" ixp  iyp  izp  ivpx  ivpy  ivpz  iap  iaps  inpswarm  irhopswarm")
	
	for i = 1:length(k)
	    if !haskey(index, k[i])  index[ k[i] ] = 0  end
	end
	
	return index
end
