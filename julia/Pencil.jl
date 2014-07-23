module Pencil
	
	import Fortran
	
	#
	# Files that provide the external API.
	#
	include("src/param.jl")
	include("src/averages.jl")
	include("src/particles.jl")
	include("src/dimensions.jl")
	include("src/timeseries.jl")
	
	#
	# Files that provide only internal functions.
	#
	include("src/index.jl")
	
end
