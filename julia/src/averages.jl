# --------------------------------------------------------------------------
# Author: Daniel Carrera (danielc@astro.lu.se)
#         Based on Python xyaver.py and yaver.py
# 
# Date:   December 2013.
# --------------------------------------------------------------------------
#>doc
# ## Read averages
# 
# Source: `julia/src/averages.jl`
# 
# ### Provides:
# 
# Functions to read 2D and 1D averages. These functions return all variable averages
# as a `Dict()` unless the `var` option is specified. They also return all time steps,
# unless the option `it` (for "*iteration*") is specified.
# 
# | Function        | Summary                                                      |
# |-----------------|--------------------------------------------------------------|
# | `read_xyaver()` | Read `xyaverages.dat` and return the variables as a `Dict()` |
# | `read_xzaver()` | Read `xzaverages.dat` and return the variables as a `Dict()` |
# | `read_yzaver()` | Read `yzaverages.dat` and return the variables as a `Dict()` |
# | `read_xaver()`  | Read `xaverages.dat`  and return the variables as a `Dict()` |
# | `read_yaver()`  | Read `yaverages.dat`  and return the variables as a `Dict()` |
# | `read_zaver()`  | Read `zaverages.dat`  and return the variables as a `Dict()` |
# 
# ### Tutorial:
# 
# ```python
#   julia> using Pencil
#   
#   julia> xyaver = read_xyaver();
#   
#   #
#   # Self-inspection. List of keys in xyaver.
#   #
#   julia> collect( keys(xyaver) )
#    "uxmz"  
#    "ux2mz" 
#    "oumz"  
#    "uz2mz" 
#    "uymz"  
#    "keys"  
#    "uy2mz" 
#    "uxuymz"
#    "rhomz" 
#    "uzmz"  
#    "t"     
#    "rhopmz"
#   
#   #
#   # The key "keys" gives you a list of the available keys / variables read.
#   #
#   julia> xyaver["keys"]
#    "t"     
#    "uxmz"  
#    "uymz"  
#    "uzmz"  
#    "ux2mz" 
#    "uy2mz" 
#    "uz2mz" 
#    "uxuymz"
#    "oumz"  
#    "rhomz" 
#    "rhopmz"
#   
#   #
#   # size() returns (rows,columns). (3416,) is a column vector.
#   #
#   julia> size( xyaver["t"] ) # (rows,columns)
#   (3416,)
#   
#   #
#   # nit    == number of iterations
#   # ngridz == number of cells in z
#   # 
#   julia> nit, ngridz = size( xyaver["rhomz"] )
#   (3416,128)
# ```
# 
# ### Optional parameters:
# 
# | Option          | Description                                          |
# | --------------- |------------------------------------------------------|
# | `var="rhopmax"` | Read only one variable instead of the entire Dict(). |
# | `datadir="xxx"` | Path to the data directory (default: "data").        |
# | `varfile="xxx"` | Name of the data file (default: "xyaverages.dat").   |
# | `infile="xxx"`  | Variables names file (default: "xyaver.in").         |
# | `it=12`         | If n > 0, return only that iteration (starts at 1).  |
# 
# ### TODO:
# 
# Need to implement `read_yaver(var="rhopmz")` and `read_xyaver(it=10)`
# 
#>end


#################
#
# TWO-DIMENSIONAL AVERAGES
#
#################

function read_xyaver(;datadir="data",varfile="xyaverages.dat",infile="xyaver.in",var="")
    read_2d_aver("z",datadir,varfile,infile,var)
end

function read_xzaver(;datadir="data",varfile="xzaverages.dat",infile="xzaver.in",var="")
    read_2d_aver("y",datadir,varfile,infile,var)
end

function read_yzaver(;datadir="data",varfile="yzaverages.dat",infile="yzaver.in",var="")
    read_2d_aver("x",datadir,varfile,infile,var)
end

function read_2d_aver(axis,datadir,varfile,infile,var)
    #
    # Best to avid readlines( ... varfile ... ) to conserve memory.
    #
	aver = Dict()
	vars = map( strip , readlines( open(datadir * "/../" * infile,"r") ) )
	dim  = read_dim( datadir=datadir )
	
	ngrid   = dim["n"*axis]                       # Size of the grid.
	nvars   = length(vars)                        # Number of variables in infile.
	nlines  = countlines(datadir * "/" * varfile) # Number of lines in varfile.
	allvars = var == ""                           # Return all varibles as a Dict.
	
	
	# Each record has 1 line + data in 8 columns.
	#
	# WARNING: This formula is different from the one in xyaver.py
	#          Why does xyaver.py use "nz % 8" instead of "(n_var*nz) % 8" ?
	reclength::Int = ceil( nvars*ngrid/8 )
	nrecords::Int  = nlines / (1 + reclength)
	
	########
	#
	# Read the records
	#
	########
	aver["keys"] = [ "t", vars ]
	aver["t"] = zeros(nrecords)
	
	#
	# Allocate only the memory we need.
	#
	if allvars
    	for var = vars
	        aver[var] = zeros(nrecords,ngrid)
    	end
	else
	    aver[var] = zeros(nrecords,ngrid)
	end
	
	#
	# Store only the data we need.
	#
	file = open(datadir * "/" * varfile  ,"r")
	for i = 1:nrecords
	    aver["t"][i] = float( readline( file ) )
	    
	    str = ""
	    for k = 1:reclength
	        str = str * chomp(readline( file ))
	    end
	    
	    dat = float(split( str ))
	    
	    for k = 1:nvars
	        if allvars | (var == vars[k])
    	        j = 1 + (k-1)*ngrid
	            aver[ vars[k] ][i,:] = dat[j:k*ngrid]
	        end
	    end
	end
	
    #
    # Note: I choose to not transpose the data. I realize that this means an
    #       extra character when plotting:
    #
    #           plot( x, rhopmx[t,:]' )
    #
    #       But transposing would mean putting the time on the second index:
    #
    #           plot( x, rhopmx[:,t] )
    #
    #       I think the two options are equally good and equally bad.
    # 
	return allvars ? aver : aver[var]
end

export read_xyaver
export read_xzaver
export read_yzaver



#################
#
# ONE-DIMENSIONAL AVERAGES
#
#################


function read_xaver(;datadir="data",varfile="xaverages.dat",infile="xaver.in",it::Int=-1)
    read_1d_aver("y","z","x",datadir,varfile,infile,it)
end
function read_yaver(;datadir="data",varfile="yaverages.dat",infile="yaver.in",it::Int=-1)
    read_1d_aver("x","z","y",datadir,varfile,infile,it)
end
function read_zaver(;datadir="data",varfile="zaverages.dat",infile="zaver.in",it::Int=-1)
    read_1d_aver("x","y","z",datadir,varfile,infile,it)
end

#
#  read_xaver -->  (y,z ; x)  == (a,b ; c)
#  read_yaver -->  (x,z ; y)  == (a,b ; c)
#  read_zaver -->  (x,y ; z)  == (a,b ; c)
#
function read_1d_aver(a,b,c,datadir,varfile,infile,it::Int)
    #
    # Global dims.
    #
    dim = read_dim(datadir=datadir)
    data_type = dim["precision"] == "D" ? Float64 : Float32
    head_type = Int32  # Header type.
    
    #
    # Variables.
    #
    vars  = map( chomp , readlines( open(datadir * "/../" * infile,"r") ) )
    nvars = length(vars)
    
    #
    # Derived parameters.
    #
    na = dim["n$a"]
    nb = dim["n$b"]
    nproc = dim["nproc$a"] * dim["nproc$b"] * dim["nproc$c"]
    size_a::Int = na / dim["nproc$a"] # Size of a a-slice (in grid cells).
    size_b::Int = nb / dim["nproc$b"] # Size of a b-slice (in grid cells).
    
    #
    # Size of time and data segments, and an entire step/record (in bytes).
    #
    head_len = sizeof(head_type)
    time_len = sizeof(data_type) * 1
    data_len = sizeof(data_type) * size_a * size_b * nvars
    step_len = head_len * 4  +  time_len  +  data_len
    
    nsteps::Int = filesize(datadir * "/proc0/" * varfile) / step_len
    
    #
    # Read the data.
    #
    aver = Dict()
    for k = 1:nvars
        #   ------>   IF (it < 1) THEN (read all iterations) ELSE (read only iteration "it")
        aver[vars[k]] =   it < 1    ?  zeros(na, nb, nsteps)   :    zeros(na, nb)
    end
    
    
    for proc = 0:(nproc-1)
        
        dim = read_dim(datadir=datadir, proc=proc)
        
        a1 = dim["ip$a"] * size_a + 1 # Start this a-slice.
        b1 = dim["ip$b"] * size_b + 1 # Start this b-slice.
        a2 = a1 + size_a - 1          # End this a-slice.
        b2 = b1 + size_b - 1          # End this b-slice.
        
        stream = open(datadir * "/proc$proc/" * varfile,"r")
        
        if it < 1
            #
            # Read all slices.
            #
            for i = 1:nsteps
                # Skip time.
                skip(stream, time_len + 2*head_len)
                
                # Does the header give the data length that I expect?
                rec_len = read(stream, head_type)
                if data_len != rec_len
                    error("Shape does not match record length. Expected $data_len b, received $rec_len b")
                end
                
                # Read the data.
                for k = 1:nvars
                     aver[vars[k]][a1:a2,b1:b2,i] = read(stream,data_type,size_a,size_b)
                end
                
                # Skip repeated header
                skip(stream, sizeof(head_type) )
            end
        else
            #
            # Read only one time step.
            #
            
            # Example: it == 1 reads the very first record; skips only the time record.
            skip(stream, step_len*(it-1) + time_len + 2*head_len)
            
            # Does the header give the data length that I expect?
            rec_len = read(stream, head_type)
            if data_len != rec_len
                error("Shape does not match record length. Expected $data_len b, received $rec_len b")
            end
            
            # Read the data.
            for k = 1:nvars
                aver[ vars[k]][a1:a2,b1:b2] = read(stream,data_type,size_a, size_b)
            end
        end
        
        close(stream)
    end
    
    #
    # Transpose the data.
    #
#    for k = 1:nvars
#        aver[vars[k]] = aver[vars[k]]'
#    end
    
    #
    # List of keys.
    #
    aver["keys"] = vars
    
    return aver
end


export read_xaver
export read_yaver
export read_zaver

