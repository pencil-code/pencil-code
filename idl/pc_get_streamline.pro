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
;   * periodic       3D-array of periodicity flags (Default: no periodicity).
;   * precision      Precision of streamline tracing between grid points, a value
;                    of 0.1 results in 10 interpolations per grid distance (Default: 0.1).
;   * max_length     Maximum length of a streamline.
;   * length         Returns the full length of a traced streamline.
;   * coords         Returns an array of grid coordinates of the traced streamline.
;   * distance       Returns an array of the distances from the anchor point along the streamline.
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
;   IDL> indices = pc_get_streamline (B, anchor=[2.0, 3.5, 1.2], grid=grid, distance=distance, length=length)
;


; Calculation of streamline coordinates.
function pc_get_streamline, field, anchor=anchor, grid=grid, distance=distance, coords=coords, direction=dir, periodic=periodic, precision=precision, length=length, max_length=max_length, cache=cache

	common pc_get_streamline_common, data, nx, ny, nz, mx, my, mz

	default, precision, 0.1
	default, nghost, 3
	default, nbox, 3.0

	if (n_elements (anchor) eq 0) then message, "ERROR: no anchor point for streamline given."
	if (n_elements (field) eq 0) then begin
		if (n_elements (data) le 1) then message, "ERROR: no vector field given."
	end else begin
		; Size of data array
		s = size (field)
		nx = s[1]
		ny = s[2]
		nz = s[3]
		mx = nx + 2*nghost
		my = ny + 2*nghost
		mz = nz + 2*nghost
		; Extent data with periodic boundary
		data = dblarr (nx+1, ny+1, nz+1, 3)
		data[0:nx-1,0:ny-1,0:nz-1,*] = field
		data[nx,0:ny-1,0:nz-1,*] = field[0,*,*,*]
		data[0:nx-1,ny,0:nz-1,*] = field[*,0,*,*]
		data[0:nx-1,0:ny-1,nz,*] = field[*,*,0,*]
		data[nx,ny,0:nz-1,*] = field[0,0,*,*]
		data[0:nx-1,ny,nz,*] = field[*,0,0,*]
		data[nx,0:ny-1,nz,*] = field[0,*,0,*]
		data[nx,ny,nz,*] = field[0,0,0,*]
	end

	default, dir, 0
	if (dir eq 0) then begin
		; Combine forward and backward streamlines from starting point
		along = pc_get_streamline (field, anchor=anchor, grid=grid, distance=d1, coords=along_coords, direction=1, periodic=periodic, precision=precision, length=l1, max_length=max_length, /cache)
		against = pc_get_streamline (field, anchor=anchor, grid=grid, distance=d2, coords=against_coords, direction=-1, periodic=periodic, precision=precision, length=l2, max_length=max_length, /cache)
		length = l1 + l2
		if (n_elements (d2) le 1) then begin
			distance = d1
			coords = along_coords
			return, along
		endif
		distance = [ -reverse (d2[1:*]), d1 ]
		against = against[*,1:*]
		against_coords = against_coords[*,1:*]
		if (size (against, /n_dim) eq 2) then against = reverse (against, 2)
		if (size (against_coords, /n_dim) eq 2) then against_coords = reverse (against_coords, 2)
		coords = [ [against_coords], [along_coords] ]
		return, [ [against], [along] ]
	end
	if (dir lt 0) then dir = -1 else dir = 1

	; Periodicity
	if (keyword_set (grid)) then periodic = grid.lperi
	default, periodic, [ 0, 0, 0 ]
	if (n_elements (periodic) eq 1) then periodic = [ periodic, periodic, periodic ]
	periodic = periodic eq 1

	; Box coordinates for lower and upper corner
	Box_xyz_lower = [ 0, 0, 0 ] - periodic * [ 0.5, 0.5, 0.5 ]
	Box_xyz_upper = [ nx, ny, nz ] - periodic * [ 0.5, 0.5, 0.5 ]

	; Grid coordinates
	if (keyword_set (grid)) then begin
		x = grid.x
		y = grid.y
		z = grid.z
		dx = 1.0 / grid.dx_1
		dy = 1.0 / grid.dy_1
		dz = 1.0 / grid.dz_1
		dr = sqrt (grid.dx^2 + grid.dy^2 + grid.dz^2)
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
	max_len = max_length / dr

	; Apply periodic boundaries to starting position
	anchor += [ nx, ny, nz ] * periodic * (anchor lt Box_xyz_lower)
	anchor -= [ nx, ny, nz ] * periodic * (anchor ge Box_xyz_upper)

	; Starting position
	pos = anchor
	if (keyword_set (grid)) then begin
		if (nx+6 ne n_elements (x)) then message, "ERROR: the field data doesn't fit to the X-grid coordinates."
		if (ny+6 ne n_elements (y)) then message, "ERROR: the field data doesn't fit to the Y-grid coordinates."
		if (nz+6 ne n_elements (z)) then message, "ERROR: the field data doesn't fit to the Z-grid coordinates."
		; Convert anchor point into equidistant unit grid coordinates
		pos[0] = pc_find_index (anchor[0], x, mx) - nghost
		pos[1] = pc_find_index (anchor[1], y, my) - nghost
		pos[2] = pc_find_index (anchor[2], z, mz) - nghost
	end

	; Iterate finding points on the streamline
	length = 0.0
	indices = [ [pos] ]
	coords = [ [anchor] ]
	distance = [ length ]
	last = anchor
	done = 0
	while (all (((pos ge 0) and (pos le ([nx, ny, nz]-1))) or periodic) and (length lt max_len) and not done) do begin

		; Interpolate data
		vector_x = interpolate (data[*,*,*,0], pos[0] + nghost, pos[1] + nghost, pos[2] + nghost)
		vector_y = interpolate (data[*,*,*,1], pos[0] + nghost, pos[1] + nghost, pos[2] + nghost)
		vector_z = interpolate (data[*,*,*,2], pos[0] + nghost, pos[1] + nghost, pos[2] + nghost)
		vector_abs = sqrt (vector_x^2 + vector_y^2 + vector_z^2)

		; Find projected step size
		dx_local = interpolate (dx, pos[0] + nghost)
		dy_local = interpolate (dy, pos[1] + nghost)
		dz_local = interpolate (dz, pos[2] + nghost)
		step_x = (vector_x / vector_abs) * precision * dir * (dx_local / dr)
		step_y = (vector_y / vector_abs) * precision * dir * (dy_local / dr)
		step_z = (vector_z / vector_abs) * precision * dir * (dz_local / dr)

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
		distance = [ distance, length ]

		; Add indices and grid coordinates to the list of traced streamline points
		indices = [ [indices], [pos] ]
		point = pos
		if (keyword_set (grid)) then begin
			point[0] = interpolate (x, pos[0] + nghost)
			point[1] = interpolate (y, pos[1] + nghost)
			point[2] = interpolate (z, pos[2] + nghost)
		end
		coords = [ [coords], [point] ]
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

