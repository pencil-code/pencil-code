# --------------------------------------------------------------------------
# Author: Daniel Carrera (danielc@astro.lu.se)
#         Based on Perl pc_plot_ts
# 
# Date:   December 2013.
# --------------------------------------------------------------------------
#>doc
# ## Read time series
# 
# Source: `julia/src/timeseries.jl`
# 
# ### Provides:
# 
# A function to read the time series from `time_series.dat` as a `Dict()`,
# with one key for each variable in `time_series.dat`.
# 
#
# | Function     | Summary                                           |
# |--------------|---------------------------------------------------|
# | `read_ts()`  | Returns a `Dict()` with keys like `rhopmax`, etc. |
# 
# ### Tutorial:
# 
# ```python
#   julia> using Pencil
#   
#   julia> ts = read_ts();
#   
#   julia> typeof(ts)
#   Dict{Any,Any} (constructor with 2 methods)
#   
#   #
#   # List of available keys / variables found in timeseries.dat
#   #
#   julia> ts["keys"]
#   35-element Array{String,1}:
#    "it"       
#    "t"        
#    "dt"       
#    "walltime" 
#    "nblockmin"
#    "nblockmax"
#    "nmigmax"  
#    "nparmin"  
#    "nparmax"  
#    "nparbmax" 
#    "rhom"     
#    "rhomax"   
#    "rhomin"   
#    "dtdragp"  
#    "urms"     
#    ...          
#    "uxuym"    
#    "TTm"      
#    "TTmax"    
#    "vpxm"     
#    "vpx2m"    
#    "vpz2m"    
#    "xpm"      
#    "xp2m"     
#    "zpm"      
#    "zp2m"     
#    "npmax"    
#    "rhopm"    
#    "rhopmax"  
#    "dedragp"  
#    "decollp"  
#   
#   julia> typeof(ts["rhopmax"])
#   Array{Float64,1}
#   
#   julia> length(ts["rhopmax"])
#   9528
#   
#   #
#   # Plotting with Python's Matplotlib.
#   #
#   julia> using PyPlot
#   
#   julia> plot( ts["t"] , ts["rhopmax"] )
# ```
# 
# ### Optional parameters:
# 
# | Option          | Description                                    |
# | --------------- |------------------------------------------------|
# | `datadir="xxx"` | Path to the data directory (default: "data").  |
# 
#>end

function read_ts(;datadir="data")
    # 
    # I could have used a DataFrame but I wanted to minimize dependencies.
    #
    ts = Dict()
    
    lines = open(readlines, datadir * "/time_series.dat", "r")
    
    #
    # The first row is a header.
    #
    header = split(lines[1],r"-+")[2:end-1]
    ncols  = length(header)
    
    #
    # There may be more headers, so I define data conservatively.
    #
    data = zeros( length(lines) , ncols );
    
    nrows = 0
    for line = lines
        if !beginswith(line,"#")
            nrows += 1
            data[nrows,:] = float(split(line))
        end
    end
    for j = 1:ncols
        ts[ header[j] ] = data[1:nrows,j]
    end
    
    #
    # Add the dust-to-gas ratio and solid-mass-fraction at each step.
    #
    param = read_param(datadir=datadir)
    dtog0 = parsefloat(param["eps_dtog"])
    ts["dtog"] = dtog0 * ts["rhom"][1] ./ ts["rhom"]
    ts["Z"] = 1 ./ (1 + 1./ts["dtog"])
    
    #
    # Finish up.
    #
    ts["keys"] = [ header , "dtog" ]
    
    return ts
end

export read_ts
