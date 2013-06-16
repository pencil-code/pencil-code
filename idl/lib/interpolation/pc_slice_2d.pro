;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_slice_2d.pro     ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  $Id$
;
;  Description:
;   Extracts any 2D slice from a given 3D data array by interpolation.
;   This method implements the idea of Hermann Amandus Schwarz et al. (1865),
;   and extends it from a 2D minimal bilinear surface solution to an
;   interpolation in 3D between source data points given on a cartesian grid.
;   Even though we cheat here a bit for simplicity, a real minimal bilinear
;   surface would be more complicated than this edge-interpolation method.
;
;  Parameters:
;   * in             1D or multi-dimensional data array
;   * source         original grid structure
;   * slice_gird     (optional) returns the grid coordinates of the slice
;   * anchor         anchor point grid coordinates [array of three components]
;   * theta          rotation angle around Z-axis (0°-360°)
;   * phi            rotation angle around the horizontal axis (0°-360°)
;
;   The horizontal axis for phi-rotation is defined by the cut of the XZ-plane
;   after theta-rotation with the XY-plane going though the anchor point.
;   For phi={0°|180°} the cut is {parallel|anti-parallel} to the Z-axis.
;   
;
;  Example:
;  ========
;
;   IDL> pc_read_var_raw, obj=var, tags=tags, grid=grid, dim=dim
;   IDL> Temp = pc_get_quantity ('Temp', var, tags)
;   IDL> anchor = [ 1.2, 3.4, 5.6 ]
;   IDL> theta = 75.6
;   IDL> phi = 12.3
;   IDL> Temp_slice = pc_slice_2d (Temp, grid, anchor, theta, phi, slice_grid=slice_grid, dim=dim)
;


; Extract a 2D slice from 3D data cube.
function pc_slice_2d, in, source, anchor, theta, phi, slice_grid=target, grid=in_grid, dim=dim, zoom=zoom

	if ((n_elements (in) eq 0) or (n_elements (source) eq 0)) then begin
		; Print usage
		print, "USAGE:"
		print, "======"
		print, "slice = pc_slice_2d (data, source_grid, anchor, theta, phi, slice_grid=slice_grid)"
		print, ""
		return, -1
	end

	if (n_elements (dim) eq 0) then pc_read_dim, obj=dim
	if (n_elements (zoom) eq 0) then zoom = 1

	NaN = !Values.D_NaN
	pi_180 = double (!Pi) / 180

	in_size = size (in)
	n_dim = in_size[0]
	if (n_dim ne 3) then message, "ERROR: the input data must be 3-dimensional."

	x_size = in_size[1]
	y_size = in_size[2]
	z_size = in_size[3]

	if (n_elements (in_grid) eq 0) then in_grid = { x:dindgen (x_size), y:dindgen (y_size), z:dindgen (z_size), lequidist:[1,1,1], lperi:[0,0,0], dx_1:1.0, dy_1:1.0, dz_1:1.0 }

	if (in_size[1] eq dim.nxgrid) then nghostx = 0 else nghostx = dim.nghostx
	if (in_size[2] eq dim.nygrid) then nghosty = 0 else nghosty = dim.nghosty
	if (in_size[3] eq dim.nzgrid) then nghostz = 0 else nghostz = dim.nghostz

	ax = anchor[0]
	ay = anchor[1]
	az = anchor[2]

	; Construct output slice
	max_size = max ([ x_size, y_size, z_size ])
	Lx = in_grid.x[x_size-1-nghostx] - in_grid.x[nghostx] + in_grid.lperi[0] * mean (in_grid.dx)
	Ly = in_grid.y[y_size-1-nghosty] - in_grid.y[nghosty] + in_grid.lperi[1] * mean (in_grid.dy)
	Lz = in_grid.z[z_size-1-nghostz] - in_grid.z[nghostz] + in_grid.lperi[2] * mean (in_grid.dz)
	L_diagonal = sqrt (Lx^2 + Ly^2 + Lz^2)
	d_min = max ([ mean (in_grid.dx), mean (in_grid.dy), mean (in_grid.dz) ]) / zoom
	num_points = ceil (L_diagonal / d_min)
	num_points += 1 - (num_points mod 2)
	out = make_array (num_points, num_points, /double, value=NaN)

	; Construct equidistant grid coordinates of input data
	Lx_i = in_grid.x[x_size-1-nghostx] - in_grid.x[nghostx]
	Ly_i = in_grid.y[y_size-1-nghosty] - in_grid.y[nghosty]
	Lz_i = in_grid.z[z_size-1-nghostz] - in_grid.z[nghostz]
	L_max = max ([ Lx_i, Ly_i, Lz_i ])
	mid = (num_points - 1) / 2
	grid = (dindgen (num_points) - mid) / (max_size-1) * L_max

	; Construct grid coordinates of output slice
	target = make_array ([ num_points, num_points, 3 ], /double, value=NaN)

	; Rotate XZ-plane around Z-axis (theta)
	cos_theta = cos (theta * pi_180)
	sin_theta = sin (theta * pi_180)
	x = grid * cos_theta
	y = grid * sin_theta
	for pv = 0, num_points-1 do begin
		; X and Y coordinates
		target[*,pv,0] = x
		target[*,pv,1] = y
	end

	; Rotate theta-rotated XZ-plane around theta-rotated x-axis (phi)
	cos_phi = cos (phi * pi_180)
	sin_phi = sin (phi * pi_180)
	add_x = grid * sin_phi * sin_theta
	add_y = grid * sin_phi * cos_theta
	z = grid * cos_phi
	for ph = 0, num_points-1 do begin
		; X and Y coordinates
		target[ph,*,0] -= add_x
		target[ph,*,1] += add_y
		; Z coordinates
		target[ph,*,2] = z
	end

	; Shift to the anchor position
	target[*,*,0] += ax
	target[*,*,1] += ay
	target[*,*,2] += az

	; soap-film interpolation to the target slice grid coordinates
	for ph = 0, num_points-1 do begin
		for pv = 0, num_points-1 do begin

			; Find grid coordinates closest to the target coordinate
			px = pc_find_index (target[ph,pv,0], in_grid.x, num=x_size, /round)
			py = pc_find_index (target[ph,pv,1], in_grid.y, num=y_size, /round)
			pz = pc_find_index (target[ph,pv,2], in_grid.z, num=z_size, /round)

			; Target position must still be inside the source data array
			if ((px lt 0) or (py lt 0) or (pz lt 0)) then continue
			if ((px ge x_size-1) or (py ge y_size-1) or (pz ge z_size-1)) then continue

			; Find lower neighbouring grid coordinates of the target coordinate
			if ((px ge 1) and (target[ph,pv,0] lt in_grid.x[px])) then px -= 1
			if ((py ge 1) and (target[ph,pv,1] lt in_grid.y[py])) then py -= 1
			if ((pz ge 1) and (target[ph,pv,2] lt in_grid.z[pz])) then pz -= 1

			; Upper plane can be interpolated by using an inner grid cell
			if (target[ph,pv,0] eq in_grid.z[x_size-1]) then px -= 1
			if (target[ph,pv,1] eq in_grid.z[y_size-1]) then py -= 1
			if (target[ph,pv,2] eq in_grid.z[z_size-1]) then pz -= 1

			; Target position must still be inside the source data array
			if ((px lt 0) or (py lt 0) or (pz lt 0)) then continue

			; Fraction of the target position to the side-length of the grid cell
			fx = (target[ph,pv,0] - in_grid.x[px]) / (in_grid.x[px+1] - in_grid.x[px])
			fy = (target[ph,pv,1] - in_grid.y[py]) / (in_grid.y[py+1] - in_grid.y[py])
			fz = (target[ph,pv,2] - in_grid.z[pz]) / (in_grid.z[pz+1] - in_grid.z[pz])

			; Find interpolated values of the target projection to the 12 edges
;			xl_yl = (1.0 - fz) * in[px  ,py  ,pz  ] + fz * in[px  ,py  ,pz+1]
;			xl_yu = (1.0 - fz) * in[px  ,py+1,pz  ] + fz * in[px  ,py+1,pz+1]
;			xu_yl = (1.0 - fz) * in[px+1,py  ,pz  ] + fz * in[px+1,py  ,pz+1]
;			xu_yu = (1.0 - fz) * in[px+1,py+1,pz  ] + fz * in[px+1,py+1,pz+1]
			xl_zl = (1.0 - fy) * in[px  ,py  ,pz  ] + fy * in[px  ,py+1,pz  ]
			xl_zu = (1.0 - fy) * in[px  ,py  ,pz+1] + fy * in[px  ,py+1,pz+1]
			xu_zl = (1.0 - fy) * in[px+1,py  ,pz  ] + fy * in[px+1,py+1,pz  ]
			xu_zu = (1.0 - fy) * in[px+1,py  ,pz+1] + fy * in[px+1,py+1,pz+1]
;			yl_zl = (1.0 - fx) * in[px  ,py  ,pz  ] + fx * in[px+1,py  ,pz  ]
;			yl_zu = (1.0 - fx) * in[px  ,py  ,pz+1] + fx * in[px+1,py  ,pz+1]
;			yu_zl = (1.0 - fx) * in[px  ,py+1,pz  ] + fx * in[px+1,py+1,pz  ]
;			yu_zu = (1.0 - fx) * in[px  ,py+1,pz+1] + fx * in[px+1,py+1,pz+1]

			; Find interpolated values of the target projection to the 6 sides
			; For the soap film, there are two identical solutions for interpolation (commented out duplicate ones)
			xy_l = (1.0 - fx) * xl_zl + fx * xu_zl
			xy_u = (1.0 - fx) * xl_zu + fx * xu_zu
;			xy_l = (1.0 - fy) * yl_zl + fy * yu_zl
;			xy_u = (1.0 - fy) * yl_zu + fy * yu_zu
;			xz_l = (1.0 - fx) * xl_yl + fx * xu_yl
;			xz_u = (1.0 - fx) * xl_yu + fx * xu_yu
;			xz_l = (1.0 - fz) * yl_zl + fz * yl_zu
;			xz_u = (1.0 - fz) * yu_zl + fz * yu_zu
;			yz_l = (1.0 - fy) * xl_yl + fy * xu_yl
;			yz_u = (1.0 - fy) * xl_yu + fy * xu_yu
;			yz_l = (1.0 - fz) * xl_zl + fz * xl_zu
;			yz_u = (1.0 - fz) * xl_zl + fz * xl_zu

			; Find interpolated value of the target point
			out[ph,pv] = (1.0 - fz) * xy_l + fz * xy_u
;			out[ph,pv] = (1.0 - fy) * xz_l + fy * xz_u
;			out[ph,pv] = (1.0 - fx) * yz_l + fx * yz_u
;			*** WORK HERE:
;			*** Need to check if these three values really are identical (theoretically, they should).
		end
	end

	return, out

end

