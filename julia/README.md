# Julia module for the Pencil Code

**Author**:   Daniel Carrera (danielc@astro.lu.se)

**Date**:   Last modified on July 2014.

## Introduction

Julia is a high-level, high-performance dynamic programming language
for technical computing, with syntax that is familiar to users of other
technical computing environments. You can obtain Julia from [http://julialang.org],
or if you are using Ubuntu/Debian, you can install it with `apt-get install julia`.

This is the documentation for the Julia module for the Pencil Code. This module
contains convenience functions for post-processing data files from the Pencil
Code. To use this module, add this to your `~/juliarc.jl` so that Julia can find
the module.

    push!(LOAD_PATH,  ENV["PENCIL_HOME"] * "/julia")


At this point you can load the `Pencil` module:

	ubuntu$ julia        # Start the Julia interpreter.
	...
	julia> using Pencil  # Load the module.

**NOTE:** At the present time, you also need to add the `push!` line at the top of stand-alone programs.


## Plotting and graphics

Julia has several plotting packages. The one I like is
[PyPlot](https://github.com/stevengj/PyPlot.jl) which uses Python's
[Matplotlib](http://matplotlib.org/index.html) library. To install
`PyPlot` run the following:

```
	ubuntu$ sudo apt-get install python-matplotlib
	
	ubuntu$ julia
	...
	julia> Pkg.add("PyPlot")
```

Brief `PyPlot` tutorial:

````python
	julia> using PyPlot
	
	julia> x = linspace(0,2*pi,1000);
	julia> y = sin(3*x + 4*cos(2*x));
	
	julia> plot(x, y, color="red", linewidth=2.0, linestyle="--")

	julia> title("A sinusoidally modulated sinusoid")
	
	#
	# Save figure as a PNG:
	#
	julia> savefig("myfigure.png")
	
	#
	# LaTeX in labels and titles.
	#
	julia> title(L"Plot of $\Gamma_3(x)$")  # L for LaTeX.
	
	#
	# Colour mesh: plot a 2D grid.
	#
	julia> y = [1:128] * ones(128)';  # Col vector x Row vector
	
	julia> r2 = (y - 64).^2 + (y' - 64).^2;
	
	julia> pcolormesh(r2)
	
	julia> axis([0,128,0,128]) # [xmin, xmax, ymin, ymax]
	
	julia> savefig("colmesh.png")
	
	#
	# 3D plotting
	#
	julia> surf(r2)
	
	julia> mesh(r2)
````

The `PencilPlot` module provides a color map called "density" that may
be useful when plotting a density variable like `rhopmxz`. `PencilPlot`
loads `PyPlot`. Example usage:

```
	ubuntu$ cd /my/simulation/dir
	ubuntu$ julia
	...
	julia> using Pencil
	
	julia> using PencilPlot
	
	julia> rhopmxz = read_yaver(it=10)["rhopmxz"];
	
	julia> axis([0, 128, 0, 128])
	
	#
	# PencilPlot defines the "density" color map.
	#
	julia> pcolormesh(rhopmxz, cmap=ColorMap("density") )
	
	julia> savefig("rhopmxz.png")
```




