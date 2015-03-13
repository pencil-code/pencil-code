# --------------------------------------------------------------------------
#     read_dim()        -- Read data/dim.dat
#     read_dim(proc=-1) -- Read data/dim.dat
# --------------------------------------------------------------------------
#     read_pdim()        -- Read data/pdim.dat
#     read_pdim(proc=-1) -- Read data/pdim.dat
#     read_pdim(proc=10) -- Read data/proc10/pdim.dat
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# Author: Daniel Carrera (danielc@astro.lu.se)
#         Based on Python dim.py and pdim.py
# 
# Date:   December 2013.
# --------------------------------------------------------------------------
#>doc
# ## Read dimensions
# 
# Source: `julia/src/dimensions.jl`
# 
# ### Provides:
# 
# Functions to parse `dim.dat` and `pdim.dat` and return the result as a
# `Dict()`. The first file contains general dimensions like `nx`, `mx`,
# `nghostx` (where `mx == nx + 2*nghostx`), `mvar`, `nprocx`, etc. The
# second file contains particle dimensions like `npar`.
# 
# | Function      | Summary                                  |
# |---------------|------------------------------------------|
# | `read_dim()`  | Reads `dim.dat` and returns a `Dict()`.  |
# | `read_pdim()` | Reads `pdim.dat` and returns a `Dict()`. |
# 
# ### Tutorial:
# 
# ```python
#   julia> using Pencil
#   
#   #
#   # Read data/dim.dat
#   #
#   julia> dim = read_dim();
#   
#   julia> typeof(dim)
#   Dict{Any,Any} (constructor with 2 methods)
#   
#   julia> dim["keys"]
#   33-element Array{Any,1}:
#    "mvar"     
#    "ny"       
#    "m2"       
#    "mx"       
#    "my"       
#    "nghosty"  
#    "nygrid"   
#    "nghostx"  
#    "precision"
#    "nprocx"   
#    "ipz"      
#    "nxgrid"   
#    "ipy"      
#    ...
#    "nprocz"   
#    "mxgrid"   
#    "mz"       
#    "nz"       
#    "maux"     
#    "ipx"      
#    "mzgrid"   
#   
#   #
#   # Compare data/dim.dat  vs data/proc10/dim.dat
#   #
#   
#   julia> read_dim()        # Read data/dim.dat
#   
#   julia> dim["nprocz"]
#   16
#   
#   julia> dim["ipz"]
#   -1
#   
#   julia> read_dim(proc=10) # Read data/proc10/dim.dat
#   
#   julia> dim["nprocz"]
#   -1
#   
#   julia> dim["ipz"]
#   10
#   
#   #
#   # read_pdim() has a similar API.
#   #
#   julia> pdim = read_pdim();
#   
#   julia> typeof(pdim)
#   Dict{Any,Any} (constructor with 2 methods)
#   
#   julia> pdim["keys"]
#   3-element Array{Any,1}:
#    "mpvar"     
#    "npar"      
#    "npar_stalk"
#   
#   julia> pdim["npar"]
#   16384
# ```
# 
# ### Optional parameters:
# 
# | Option          | Description                                                     |
# | --------------- |-----------------------------------------------------------------|
# | `datadir="xxx"` | Path to the data directory (default: "data").                   |
# | `proc=n`        | Read the `dim` or `pdim` file from the `data/proc$n` directory. |
# 
#>end

function read_pdim(;datadir="data",proc=-1)
	
	if proc < 0
		filename = datadir * "/pdim.dat"  # Global.
	else
		filename = datadir * "/proc$proc/pdim.dat" # Local.
	end
	
	pdim = Dict()
	line = open(readline, filename,"r")
	
	pdim["npar"],pdim["mpvar"],pdim["npar_stalk"] = int(split(line))
	
	pdim["keys"] = collect(keys(pdim))
	
	return pdim
end

function read_dim(;datadir="data",proc=-1)
	
	if proc < 0
		filename = datadir * "/dim.dat"  # Global.
	else
		filename = datadir * "/proc$proc/dim.dat" # Local.
	end
	
	dim   = Dict()
	lines = open(readlines, filename,"r")
	
	#
	# Read lines.
	#
	dim["mx"],dim["my"],dim["mz"],dim["mvar"],dim["maux"] = int(split(lines[1]))[1:5]
	dim["nghostx"], dim["nghosty"], dim["nghostz"]        = int(split(lines[3]))[1:3]
	dim["mglobal"]   = length(split(lines[1])) == 6
	dim["precision"] = chomp(lines[2])
	
	if proc < 0
    	dim["nprocx"],dim["nprocy"],dim["nprocz"],dim["iprocz_slowest"] = int(split(lines[4]))
	    dim["ipx"] = dim["ipy"] = dim["ipz"] = -1
	else
    	dim["nprocx"] = dim["nprocy"] = dim["nprocz"] = dim["iprocz_slowest"] = -1
	    dim["ipx"], dim["ipy"], dim["ipz"] = int(split(lines[4]))
	end
	
	#
	# Derived quantities.
	#
	dim["nx"] = dim["mx"] - 2*dim["nghostx"]
	dim["ny"] = dim["my"] - 2*dim["nghosty"]
	dim["nz"] = dim["mz"] - 2*dim["nghostz"]
	dim["mw"] = dim["mx"] * dim["my"] * dim["mz"]
	dim["l1"] = dim["nghostx"]
	dim["m1"] = dim["nghosty"]
	dim["n1"] = dim["nghostz"]
	dim["l2"] = dim["mx"] - dim["nghostx"] - 1
	dim["m2"] = dim["my"] - dim["nghosty"] - 1
	dim["n2"] = dim["mz"] - dim["nghostz"] - 1
	
	if proc < 0
		# Global
		dim["nxgrid"] = dim["nx"]
		dim["nygrid"] = dim["ny"]
		dim["nzgrid"] = dim["nz"]
		dim["mxgrid"] = dim["nxgrid"] + 2*dim["nghostx"]
		dim["mygrid"] = dim["nygrid"] + 2*dim["nghosty"]
		dim["mzgrid"] = dim["nzgrid"] + 2*dim["nghostz"]
	else
		# Local
		dim["nxgrid"] = dim["nygrid"] = dim["nzgrid"] = 0
		dim["mxgrid"] = dim["mygrid"] = dim["mzgrid"] = 0
	end
	
	dim["keys"] = collect(keys(dim))
	
	return dim
end

export read_dim
export read_pdim
