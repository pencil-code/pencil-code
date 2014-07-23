# --------------------------------------------------------------------------
# Module to read and write Fortran binary files.
# 
# Usage:
# 
#     import FortranIO
#     
#     file = open( "data/yaverages.dat", "r" )
#     
#     a = Fortran.read( file )
#     b = Fortran.read( file, Float64 )
#     c = Fortran.read( file, Float64, 1 )
#     d = Fortran.read( file, Float64, (nvars,nz,nx) )
#
# --------------------------------------------------------------------------
# Author: Daniel Carrera (danielc@astro.lu.se)
#         Inspired by npfile.py
# --------------------------------------------------------------------------

module Fortran

function read(stream, datatype=Float64, shape=1 ; headtype=Int32, endian::String="host")
    
    reclength = Main.read(stream, headtype)
    
    if prod(shape)*sizeof(datatype) != reclength
        error("Shape does not match record length. Expected $(prod(shape)), got $reclength")
    end
    
    # Strangely, Main.read() will not accept a shape of (Int32,Int32)
    shape = map(int64 , shape)
    
    slice = Main.read(stream, datatype, shape)
    if endian == "le"  slice = map(ltoh, Main.read(stream, datatype, shape))  end
    if endian == "be"  slice = map(ntoh, Main.read(stream, datatype, shape))  end
    
    # Fortran records end with the header repeated. Skip this.
    skip(stream, sizeof(headtype))
    
    return slice
end

end

