module PencilPlot
	
	using PyPlot

	#
	# Define a ColorMap that can be used for density plots.
	#
	
	colors = [
		[ 0   0   0  ],
		[ 0   0   100],
		[ 0   0   255],
		[ 155 0   155],
		[ 255 0   0  ],
		[ 195 195 50 ],
		[ 255 255 255]
	]
	colors = [
		[ 0   0   0  ],
		[ 0   0   255],
		[ 155 0   155],
		[ 255 0   0  ],
		[ 195 195 50 ],
		[ 255 255 255]
	]
	colors = [
		[ 0   0   0  ],
		[ 0   0   0  ],
		[ 0   0   255],
		[ 65  0   190],
		[ 130 0   125],
		[ 195 0   60 ],
		[ 255 0   0  ],
		[ 200 100 0  ],
		[ 150 100 0  ],
		[ 150 150 50 ],
		[ 250 250 250]
	]
	
	colors = colors ./ 255
	
	nc = size(colors, 1)
	
    r = Array((Float64,Float64,Float64), nc)
    g = Array((Float64,Float64,Float64), nc)
    b = Array((Float64,Float64,Float64), nc)
	
	for i = 1:nc
		j = 1.0 * (i-1) / (nc-1)
		
		r[i] = (j, colors[i,1],  colors[i,1])
		g[i] = (j, colors[i,2],  colors[i,2])
		b[i] = (j, colors[i,3],  colors[i,3])
	end
	
	cdensity = ColorMap("density", r, g, b)
	
	register_cmap("density", cdensity)
	
	export cdensity
end
