;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_get_streamline.pro     ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  $Id$
;
;  Description:
;   Calculation of coordinates along a traced streamline.
;
;  Parameters:
;   * field          Data cube of a vector field (4-dimensional: [nx,ny,nz,3]).
;   * anchor         Anchor point in grid coordinates.
;   * grid           Grid structure (Default: equidistant grid spacing of unit length 1.0).
;   * direction      Direction (1: along, -1: against the field, Default: both).
;   * periodic       3D-array of periodicity flags (Default: no periodicity, if no grid is given).
;   * precision      Precision of streamline tracing between grid points, a value
;                    of 0.1 results in 10 interpolations per grid distance (Default: 0.1).
;   * max_length     Maximum length of a streamline.
;   * length         Returns the full length of a traced streamline.
;   * num            Returns the number of points of the traced streamline.
;   * origin         Returns the array index of the origin of the traced streamline.
;   * coords         Returns an array of grid coordinates of the traced streamline.
;   * distances      Returns an array of the distances from the anchor point along the streamline.
;   * cache          Cache vector field data cube for later use.
;
;  Returns:
;   * indices        Array of data indices of the traced streamline.
;
;  Example:
;  ========
;
;   Load varfile and extract Temperature along a magnetic filedline:
;   IDL> pc_read_var_raw, obj=var, tags=tags, grid=grid
;   IDL> B = pc_get_quantity ('B', var, tags)
;   IDL> Temp = pc_get_quantity ('Temp', var, tags)
;   IDL> indices = pc_get_streamline (B, anchor=[2.0, 3.5, 1.2], grid=grid, distances=distances, length=length)
;


; Calculation of streamline coordinates.
function pc_get_streamline, field, anchor=anchor, grid=grid, distances=distances, coords=coords, direction=dir, periodic=periodic, precision=precision, length=length, num=num, origin=origin, max_length=max_length, cache=cache

	common pc_get_streamline_common, data, nx, ny, nz, mx, my, mz, Box_xyz_lower, Box_xyz_upper

	default, precision, 0.1
	default, nghost, 3
	default, nbox, 3.0

	; Periodicity
	if (keyword_set (grid)) then periodic = grid.lperi
	default, periodic, [ 0, 0, 0 ]
	if (n_elements (periodic) eq 1) then periodic = [ periodic, periodic, periodic ]
	periodic = periodic eq 1

	if (n_elements (anchor) eq 0) then message, "ERROR: no anchor point for streamline given."
	if (n_elements (field) eq 0) then begin
		if (n_elements (data) le 1) then message, "ERROR: no vector field given."
	end else begin
		; Size of new given data array
		s = size (field)
		nx = s[1]
		ny = s[2]
		nz = s[3]
		mx = nx + 2*nghost
		my = ny + 2*nghost
		mz = nz + 2*nghost

		; Box indices for lower and upper corner
		Box_xyz_lower = [ 0, 0, 0 ]
		Box_xyz_upper = [ nx, ny, nz ] + (periodic - 1)

		if (not any (periodic)) then begin
			data = field
		end else begin
			; Extend data with one ghost layer in each periodic direction
			data = dblarr (nx+periodic[0], ny+periodic[1], nz+periodic[2], 3)
			data[0:nx-1,0:ny-1,0:nz-1,*] = field
			if (periodic[0]) then data[nx,0:ny-1,0:nz-1,*] = field[0,*,*,*]
			if (periodic[1]) then data[0:nx-1,ny,0:nz-1,*] = field[*,0,*,*]
			if (periodic[2]) then data[0:nx-1,0:ny-1,nz,*] = field[*,*,0,*]
			if (periodic[0] and periodic[1]) then data[nx,ny,0:nz-1,*] = field[0,0,*,*]
			if (periodic[1] and periodic[2]) then data[0:nx-1,ny,nz,*] = field[*,0,0,*]
			if (periodic[0] and periodic[2]) then data[nx,0:ny-1,nz,*] = field[0,*,0,*]
			if (all (periodic)) then data[nx,ny,nz,*] = field[0,0,0,*]
		end
	end

	default, dir, 0
	if (dir eq 0) then begin
		; Combine forward and backward streamlines from starting point
		along = pc_get_streamline (field, anchor=anchor, grid=grid, distances=distances, coords=coords, direction=1, periodic=periodic, precision=precision, length=length, num=num, origin=origin, max_length=max_length, /cache)
		against = pc_get_streamline (field, anchor=anchor, grid=grid, distances=d2, coords=against_coords, direction=-1, periodic=periodic, precision=precision, length=l2, num=n2, max_length=max_length, /cache)
		if (n2 le 1) then return, along
		length += l2
		num += n2 - 1
		origin = n2 - 1
		distances = [ -reverse (d2[1:*]), distances ]
		against = against[*,1:*]
		against_coords = against_coords[*,1:*]
		if (size (against, /n_dim) eq 2) then against = reverse (against, 2)
		if (size (against_coords, /n_dim) eq 2) then against_coords = reverse (against_coords, 2)
		coords = [ [against_coords], [coords] ]
		return, [ [against], [along] ]
	end
	if (dir lt 0) then dir = -1 else dir = 1

	; Grid coordinates
	if (keyword_set (grid)) then begin
		x = grid.x
		y = grid.y
		z = grid.z
		dx = 1.0 / grid.dx_1
		dy = 1.0 / grid.dy_1
		dz = 1.0 / grid.dz_1
		if (any (strcmp (tag_names (grid), "nghost", /fold_case))) then begin
			nghost = grid.nghost
			mx = nx + 2*nghost
			my = ny + 2*nghost
			mz = nz + 2*nghost
		end
		dr = sqrt (mean (dx)^2 + mean (dy)^2 + mean (dz)^2)
	end else begin
		if (periodic[0]) then x = dindgen (mx) - nghost + 0.5 else x = dindgen (mx) / (mx-1) * (mx) - nghost
		if (periodic[1]) then y = dindgen (my) - nghost + 0.5 else y = dindgen (my) / (my-1) * (my) - nghost
		if (periodic[2]) then z = dindgen (mz) - nghost + 0.5 else z = dindgen (mz) / (mz-1) * (mz) - nghost
		dx = replicate (1.0, mx)
		dy = replicate (1.0, my)
		dz = replicate (1.0, mz)
		dr = 1.0
	end

	; Maximum streamline length
	if (n_elements (max_length) eq 0) then max_length = total (Box_xyz_upper - Box_xyz_lower) * nbox * dr
	if (max_length le 0.0) then message, "ERROR: 'max_length' must be positive definite."

	; Starting position
	pos = anchor
	if (keyword_set (grid)) then begin
		if (nx+2*nghost ne n_elements (x)) then message, "ERROR: the field data doesn't fit to the X-grid coordinates."
		if (ny+2*nghost ne n_elements (y)) then message, "ERROR: the field data doesn't fit to the Y-grid coordinates."
		if (nz+2*nghost ne n_elements (z)) then message, "ERROR: the field data doesn't fit to the Z-grid coordinates."
		; Convert anchor point into equidistant unit grid coordinates
		pos[0] = pc_find_index (anchor[0], x, num=mx) - nghost
		pos[1] = pc_find_index (anchor[1], y, num=my) - nghost
		pos[2] = pc_find_index (anchor[2], z, num=mz) - nghost
	end

	; Iterate finding points on the streamline
	origin = 0
	num = 0
	length = 0.0
	indices = [ [pos] ]
	coords = [ [anchor] ]
	distances = [ length ]
	last = anchor
	done = 0
	while (all (((pos ge 0) and (pos le ([nx, ny, nz]-1))) or periodic) and (length lt max_length) and not done) do begin

		; Interpolate data
		i = floor (pos)
		residual = pos - i
		tmp_data = data[i[0]:i[0]+1,i[1]:i[1]+1,i[2]:i[2]+1,*]
		vector_x = interpolate (tmp_data[*,*,*,0], residual[0], residual[1], residual[2])
		vector_y = interpolate (tmp_data[*,*,*,1], residual[0], residual[1], residual[2])
		vector_z = interpolate (tmp_data[*,*,*,2], residual[0], residual[1], residual[2])
		vector_abs = sqrt (vector_x^2 + vector_y^2 + vector_z^2)

		; Find projected step size
		dx_local = interpolate (dx, pos[0] + nghost)
		dy_local = interpolate (dy, pos[1] + nghost)
		dz_local = interpolate (dz, pos[2] + nghost)
		step_x = (vector_x / vector_abs) * precision * dir * (dr / dx_local)
		step_y = (vector_y / vector_abs) * precision * dir * (dr / dy_local)
		step_z = (vector_z / vector_abs) * precision * dir * (dr / dz_local)

		; Reduce projected step size to non-periodic boundaries
		reduce = 0.0
		if (not periodic[0]) then reduce = max ([ reduce, -(pos[0] + step_x), pos[0] + step_x - (nx-1) ] > 0.0) / abs (step_x)
		if (not periodic[1]) then reduce = max ([ reduce, -(pos[1] + step_y), pos[1] + step_y - (ny-1) ] > 0.0) / abs (step_y)
		if (not periodic[2]) then reduce = max ([ reduce, -(pos[2] + step_z), pos[2] + step_z - (nz-1) ] > 0.0) / abs (step_z)
		if (reduce gt 0.0) then begin
			step_x *= (1.0 - reduce)
			step_y *= (1.0 - reduce)
			step_z *= (1.0 - reduce)
			done = 1
			if (reduce ge 1.0) then continue
		end

		; Calculate new position
		pos += [ step_x, step_y, step_z ]

		; Add step length in the given grid coordinates
		delta_x = last[0] - interpolate (x, pos[0] + nghost)
		delta_y = last[1] - interpolate (y, pos[1] + nghost)
		delta_z = last[2] - interpolate (z, pos[2] + nghost)
		delta = sqrt (delta_x^2 + delta_y^2 + delta_z^2)

		; Apply periodic boundaries
		pos += [ nx, ny, nz ] * periodic * (pos lt Box_xyz_lower)
		pos -= [ nx, ny, nz ] * periodic * (pos ge Box_xyz_upper)

		; Compute the path length along the traced streamline
		length += delta
		distances = [ distances, length ]

		; Add indices and grid coordinates to the list of traced streamline points
		indices = [ [indices], [pos] ]
		point = pos
		if (keyword_set (grid)) then begin
			point[0] = interpolate (x, pos[0] + nghost)
			point[1] = interpolate (y, pos[1] + nghost)
			point[2] = interpolate (z, pos[2] + nghost)
		end
		coords = [ [coords], [point] ]
		num++
		last = point
	end

	if (not keyword_set (cache)) then begin
		data = 0
		nx = 0
		ny = 0
		nz = 0
		mx = 0
		my = 0
		mz = 0
	end

	return, indices

end

