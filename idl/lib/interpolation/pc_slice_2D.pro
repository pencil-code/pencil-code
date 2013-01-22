;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_slice_2D.pro     ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  $Id$
;
;  Description:
;   Extracts any 2D slice from a given 3D data array by interpolation.
;   This method implements the idea of Hermann Amandus Schwarz et al. (1865),
;   and extends it from a 2D minimal bilinear surface solution to an
;   interpolation in 3D between source data points given on a cartesian grid.
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
;   IDL> Temp_slice = pc_slice_2D (Temp, grid, anchor, theta, phi, slice_grid=slice_grid, dim=dim)
;


; Extract a 2D slice from 3D data cube.
function pc_slice_2D, in, source, anchor, theta, phi, slice_grid=target, dim=dim

	if ((n_elements (in) eq 0) or (n_elements (source) eq 0)) then begin
		; Print usage
		print, "USAGE:"
		print, "======"
		print, "slice = pc_slice_2D (data, source_grid, anchor, theta, phi, slice_grid=slice_grid)"
		print, ""
		return, -1
	end

	if (n_elements (dim) eq 0) then pc_read_dim, obj=dim

	NaN = !Values.D_NaN
	pi_180 = double (!Pi) / 180

	in_size = size (in)
	n_dim = in_size[0]
	if (n_dim ne 3) then message, "ERROR: the input data must be 3-dimensional."

	x_size = in_size[1]
	y_size = in_size[2]
	z_size = in_size[3]

	; Construct equidistant grid coordinates of input data
	xl = source.x[dim.nghostx]
	xu = source.x[dim.nghostx+x_size-1]
	yl = source.y[dim.nghosty]
	yu = source.y[dim.nghosty+y_size-1]
	zl = source.z[dim.nghostz]
	zu = source.z[dim.nghostz+z_size-1]
	x = dindgen (x_size) / (x_size-1.0) * (xu - xl) + xl
	y = dindgen (y_size) / (y_size-1.0) * (yu - yl) + yl
	z = dindgen (z_size) / (z_size-1.0) * (zu - zl) + zl

	ax = anchor[0]
	ay = anchor[1]
	az = anchor[2]

	; Construct output slice
	num_horiz = max ([ x_size, y_size ])
	num_vert = z_size
	out = make_array (num_horiz, num_vert, type=in_size[n_dim+1])

	; Construct grid coordinates of output slice
	target = make_array ([ num_horiz, num_vert, 3 ], /double, value=NaN)

	; Rotate XZ-plane around Z-axis (theta)
	for ph = 0, num_horiz-1 do begin
		; X and Y coordinates
		target[ph,*,0] = (x[ph] - ax) * cos (theta * pi_180) + ax
		target[ph,*,1] = (x[ph] - ax) * sin (theta * pi_180) + ay
		; Z remains unchanged
		target[ph,*,2] = z
	end

	; Rotate theta-rotated XZ-plane around horizontal axis (phi)
	for ph = 0, num_horiz-1 do begin
		; X and Y coordinates
		target[ph,*,0] -= (target[ph,*,2] - az) * sin (phi * pi_180) * sin (theta * pi_180)
		target[ph,*,1] += (target[ph,*,2] - az) * sin (phi * pi_180) * cos (theta * pi_180)
	end
	for pv = 0, num_vert-1 do begin
		; Z coordinates
		target[*,pv,2] -= (z[pv] - az) * (1.0 - sin ((phi+90) * pi_180))
	end

	; soap-film interpolation to the target slice grid coordinates
	for ph = 0, num_horiz-1 do begin
		for pv = 0, num_vert-1 do begin

			; Find neighbouring grid coordinates of the target coordinate
			px = pc_find_index (target[ph,pv,0], x, num=x_size, /round)
			py = pc_find_index (target[ph,pv,1], y, num=y_size, /round)
			pz = pc_find_index (target[ph,pv,2], z, num=z_size, /round)
			if (target[ph,pv,0] lt x[px]) then px -= 1
			if (target[ph,pv,1] lt y[py]) then py -= 1
			if (target[ph,pv,2] lt z[pz]) then pz -= 1
			; Upper plane can be interpolated by using an inner grid cell
			if (target[ph,pv,0] eq x[x_size-1]) then px -= 1
			if (target[ph,pv,1] eq y[y_size-1]) then py -= 1
			if (target[ph,pv,2] eq z[z_size-1]) then pz -= 1

			; Target position must still be inside the source data array
			if ((px lt 0) or (py lt 0) or (pz lt 0)) then continue
			if ((px ge x_size-1) or (py ge y_size-1) or (pz ge z_size-1)) then continue

			; Fraction of the target position to the side-length of the grid cell
			fx = (target[ph,pv,0] - x[px]) / (x[px+1] - x[px])
			fy = (target[ph,pv,1] - y[py]) / (y[py+1] - y[py])
			fz = (target[ph,pv,2] - z[pz]) / (z[pz+1] - z[pz])

			; Find interpolated values of the target projection to the 12 edges
			xl_yl = (1.0 - fz) * in[px  ,py  ,pz  ] + fz * in[px  ,py  ,pz+1]
			xl_yu = (1.0 - fz) * in[px  ,py+1,pz  ] + fz * in[px  ,py+1,pz+1]
			xu_yl = (1.0 - fz) * in[px+1,py  ,pz  ] + fz * in[px+1,py  ,pz+1]
			xu_yu = (1.0 - fz) * in[px+1,py+1,pz  ] + fz * in[px+1,py+1,pz+1]
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
			xz_l = (1.0 - fx) * xl_yl + fx * xu_yl
			xz_u = (1.0 - fx) * xl_yu + fx * xu_yu
;			xz_l = (1.0 - fz) * yl_zl + fz * yl_zu
;			xz_u = (1.0 - fz) * yu_zl + fz * yu_zu
			yz_l = (1.0 - fy) * xl_yl + fy * xu_yl
			yz_u = (1.0 - fy) * xl_yu + fy * xu_yu
;			yz_l = (1.0 - fz) * xl_zl + fz * xl_zu
;			yz_u = (1.0 - fz) * xl_zl + fz * xl_zu

			; Find interpolated value of the target point
			out[ph,pv] = (1.0 - fz) * xy_l + fz * xy_u
;			out[ph,pv] = (1.0 - fy) * xz_l + fy * xz_u
;			out[ph,pv] = (1.0 - fx) * yz_l + fx * yz_u
;			*** WORK HERE:
;			*** Need to check if these three values really are identical.
		end
	end
	

	return, out

end

