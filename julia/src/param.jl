# --------------------------------------------------------------------------
# Author: Daniel Carrera (danielc@astro.lu.se)
# 
# Date:   April 2014
# --------------------------------------------------------------------------
#>doc
# ## Read simulation parameters
# 
# Source: `julia/src/param.jl`
# 
# ### Provides:
# 
# A function to read the simulation parameters from `param.nml` and `param2.nml`
# as a `Dict()`, with one key for each parameter.
# 
#>end

function read_param(;datadir="data")
    
    param = Dict()
    
    
    #
    # Only care about param2.nml; that file contains the runtime values.
    # If there is no param2.nml, I prefer to get an error, as there is no run.
    #
    file  = open(datadir * "/param2.nml", "r")
    lines = readlines(file)
    
    #
    # Process each line.
    #
    for l in lines
        if ismatch(r"=", l)
            key, val = split(l, "=")
            key = lowercase(strip(key))
            val = strip(val)
            
            if endswith(val, ",") val = chop(val) end
            
            param[key] = val
        end
    end
    
    param["keys"] = collect(keys(param))
    
    return param
end

export read_param
