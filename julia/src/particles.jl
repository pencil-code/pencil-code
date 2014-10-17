# --------------------------------------------------------------------------
# Author: Daniel Carrera (danielc@astro.lu.se)
#         With help from Anders Johansen.
# 
# Date:   December 2013.
# --------------------------------------------------------------------------
#>doc
# ## Read particles
# 
# Source: `julia/src/particles.jl`
# 
# ### Provides:
# 
# Two functions to read particle data from `pvar.dat` or the `PVAR*` files.
# The difference is only in the output: `read_pvar` is more consistent with
# most functions in this module in that in returns a `Dict` with keys like
# `x` and `vx`, while `read_particles` returns an array of `Particle` objects
# which is more cache-friendly for numerical computation.
# 
# | Function           | Summary                                         |
# |--------------------|-------------------------------------------------|
# | `read_pvar()`      | Returns a `Dict()` with keys like `x` and `vx`. |
# | `read_particles()` | Returns an array of `Particle` objects.         |
# 
# ### Tutorial:
# 
# ```python
#   julia> using Pencil
#   
#   julia> pars = read_particles();
#   INFO: mpvar = 6
#   INFO: Read 16384 particles
#   
#   #
#   # "pars" is an array of 16384 particles.
#   #
#   julia> typeof(pars)
#   Array{Particle,1}
#   
#   julia> size(pars)
#   (16384,)
#   
#   #
#   # Each particle has a position (x,y,z) and velocity (u,v,w)
#   #
#   julia> pars[1].x
#   -0.058057133821926316
#   
#   julia> pars[1].u
#   0.0003632957506802519
#   
#   #############
#   
#   julia> pvar = read_pvar();
#   INFO: mpvar = 6
#   INFO: Read 16384 particles
#   
#   julia> typeof(pvar)
#   Dict{Any,Any} (constructor with 2 methods)
#   
#   #
#   # The "keys" key has a list of dictionary keys / variables read (from index.pro)
#   #
#   julia> pvar["keys"]
#   6-element Array{Any,1}:
#   "y" 
#   "vx"
#   "vy"
#   "x" 
#   "z" 
#   "vz"
#   
#   julia> typeof(pvar["x"])
#   Array{Float64,1}
#   
#   julia> size(pvar["x"])
#   (16384,)
# ```
# 
# ### Optional parameters:
# 
# | Option          | Description                                    |
# | --------------- |------------------------------------------------|
# | `datadir="xxx"` | Path to the data directory (default: "data").  |
# | `snapshot=n`    | If n > 0, return data from file `PVAR$n`.      |
# 
# ### TODO:
# 
# Need to implement the `proc` parameter.
# 
#>end

#################
#
# TRADITIONAL API
#
#################

function read_pvar(;datadir="data",proc=-1,snapshot=-1)
    #
    #  Read dim.dat
    #
    dim = read_dim(datadir=datadir)
    
    dtype = dim["precision"] == "D" ? Float64 : Float32
    nproc = dim["nprocx"] * dim["nprocy"] * dim["nprocz"]
    
    #
    # Number of variables stored per particle -- from pdim.dat
    #
    pdim  = read_pdim(datadir=datadir)
    mpvar = pdim["mpvar"]
    npar  = pdim["npar"]
    
    info("mpvar = $mpvar")
    
    #
    # File to read
    #
    filename = snapshot < 0 ? "pvar.dat" : "PVAR$snapshot";
    
    #
    # Read the PVAR files.
    #
    fp = zeros(npar,mpvar)
    i = 0
    for proc = 0:(nproc-1)
        stream = open( datadir * "/proc$proc" * "/" * filename , "r" )
        
        npar_loc = Fortran.read(stream, Int32)[1] # Number of particles in this proc.
        
        if npar_loc > 0
            ipar     = Fortran.read(stream, Int32, npar_loc) # Particle indices.
            fp[i+1:i+npar_loc,:] = Fortran.read(stream, dtype, (npar_loc, mpvar)) # Particle data.
            i += npar_loc
        end
        
        close(stream)
    end
    
    if i == npar
        info("Read $i particles")
    else
        warn("Expected $npar particles but got $i")
    end
    
    
    #
    # Index of each variable -- from index.pro
    #
    index = read_index(datadir=datadir)
    pvar = Dict()
    ivar = split("ixp  iyp  izp  ivpx  ivpy  ivpz  iap  iaps")
     var = split( "x    y    z    vx    vy    vz    a    as" )
    
    for i = 1:length(ivar)
        if index[ ivar[i] ] > 0
            pvar[  var[i] ] = fp[ : , index[ ivar[i] ] ]
        end
    end
    pvar["keys"] = collect(keys(pvar))
    
    return pvar
end


#################
#
# ALTERNATE API
#
#################

type Particle
    x ::FloatingPoint # x position.
    y ::FloatingPoint # y position.
    z ::FloatingPoint # z position.
    vx::FloatingPoint # vx velocity.
    vy::FloatingPoint # vy velocity.
    vz::FloatingPoint # vz velocity.
    a ::FloatingPoint # particle radius.
    as::FloatingPoint # particle sink radius.
end

function read_particles(;datadir="data",proc=-1,snapshot=-1)
    #
    # Re-arrange data from read_pvar().
    #
    pvar = read_pvar(datadir=datadir,proc=proc,snapshot=snapshot)
    npar = read_pdim(datadir=datadir)["npar"]
    
    particles = Array(Particle,npar)
    
    for i = 1:npar
        particles[i] = Particle( haskey(pvar, "x" ) ? pvar["x" ][i] : 0.0,
                                 haskey(pvar, "y" ) ? pvar["y" ][i] : 0.0,
                                 haskey(pvar, "z" ) ? pvar["z" ][i] : 0.0,
                                 haskey(pvar, "vx") ? pvar["vx"][i] : 0.0,
                                 haskey(pvar, "vy") ? pvar["vy"][i] : 0.0,
                                 haskey(pvar, "vz") ? pvar["vz"][i] : 0.0,
                                 haskey(pvar, "a" ) ? pvar["a" ][i] : 0.0,
                                 haskey(pvar, "as") ? pvar["as"][i] : 0.0)
    end
    return particles
end


export Particle
export read_pvar
export read_particles

