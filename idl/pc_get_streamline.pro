;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_get_streamline.pro     ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  $Id$
;
;  Description:
;   Calculation of coordinates along a traced streamline.
;
;  Parameters:
;   * data           Data cube of a vector field (4-dimensional: [nx,ny,nz,3]).
;   * anchor         Anchor point in grid coordinates (2-dimensional: [3,num_lines]).
;   * grid           Grid structure (Default: equidistant grid spacing of unit length 1.0).
;   * direction      Direction (1: along, -1: against the vector field, Default: both).
;   * periodic       3D-array of periodicity flags (Default: no periodicity, if no grid is given).
;   * precision      Precision of streamline tracing between grid points, a value
;                    of 0.1 results in 10 interpolations per grid distance (Default: 0.1).
;   * select         Save only every select-th point along the traced streamline.
;   * max_length     Maximum length of a streamline.
;   * length         Returns the full length of each traced streamline.
;   * num_lines      Returns the number of traced streamlines.
;   * num_points     Returns the number of points in each traced streamline.
;   * origin         Returns the array index of the origin of each traced streamline.
;   * coords         Returns an array of grid coordinates of each traced streamline.
;   * distances      Returns an array of the distances from the anchor point along each streamline.
;   * return_indices Return an array of indices of the streamline (Default: return streamlines structure).
;
;  Returns:
;   * streamlines    Streamlines structure containing all relevant data (Default).
;   * indices        Array of data indices of the traced streamline (if /return_indices is set).
;
;  Examples:
;  =========
;
;   Load varfile and extract Temperature along a magnetic filedline:
;   IDL> pc_read_var_raw, obj=var, tags=tags, grid=grid
;   IDL> B = pc_get_quantity ('B', var, tags)
;   IDL> Temp = pc_get_quantity ('Temp', var, tags)
;   IDL> indices = pc_get_streamline (B, anchor=[2.0, 3.5, 1.2], grid=grid, distances=distances, length=length, /return_indices)
;   IDL> Temp_streamline = pc_extract_streamline (Temp, indices)
;
;   Load varfile and extract Temperature along several magnetic filedlines:
;   IDL> pc_read_var_raw, obj=var, tags=tags, grid=grid
;   IDL> B = pc_get_quantity ('B', var, tags)
;   IDL> Temp = pc_get_quantity ('Temp', var, tags)
;   IDL> seeds = pc_seed_points (grid)
;   IDL> streamlines = pc_get_streamline (B, anchor=seeds, grid=grid)
;   IDL> Temp_streamlines = pc_extract_streamline (Temp, streamlines, name='Temperature')
;


; Calculation of streamline coordinates.
function pc_get_streamline, data, anchor=anchor, grid=grid, distances=distances, coords=coords, direction=dir, periodic=periodic, precision=precision, select=select, length=length, num_lines=num_lines, num_points=num_points, origin=origin, max_length=max_length, return_indices=return_indices

	default, dir, 0
	default, precision, 0.1
	default, select, 5
	default, nghost, 3
	default, nbox, 3.0
	default, max_packet_length, 1000000L

	; Periodicity
	if (keyword_set (grid)) then periodic = grid.lperi
	default, periodic, [ 0, 0, 0 ]
	if (n_elements (periodic) eq 1) then periodic = [ periodic, periodic, periodic ]
	periodic = periodic eq 1

	if (n_elements (anchor) eq 0) then message, "ERROR: no anchor point for streamline given."

	; Size of new given data array
	s = size (data)
	nx = s[1]
	ny = s[2]
	nz = s[3]
	mx = nx + 2*nghost
	my = ny + 2*nghost
	mz = nz + 2*nghost

	; Box indices for lower and upper corner
	Box_xyz_lower = [ 0, 0, 0 ]
	Box_xyz_upper = [ nx, ny, nz ] + (periodic - 1)

	if (size (anchor, /n_dimensions) gt 1) then begin
		; Iterate though list anchor points
		num_lines = long (n_elements (anchor[0,*]))
		line_pos = 0L
		first = lonarr (num_lines)
		last = lonarr (num_lines)
		num_points = lonarr (num_lines)
		packet = 0L
		packet_length = 0L
		indices_packet = dblarr (3, max_packet_length)
		coords_packet = dblarr (3, max_packet_length)
		distances_packet = dblarr (max_packet_length)
		origin = lonarr (num_lines)
		length = dblarr (num_lines)
		for pos = 0L, num_lines - 1L do begin
			stream = pc_get_streamline (data, anchor=anchor[*,pos], grid=grid, direction=dir, periodic=periodic, precision=precision, select=select, max_length=max_length)
			num_points[pos] = stream.num_points
			if ((packet_length + num_points[pos]) gt max_packet_length) then begin
				if (packet_length gt 0L) then begin
					if (packet eq 0L) then begin
						indices_stack = indices_packet[*,0L:packet_length-1L]
						coords_stack = coords_packet[*,0L:packet_length-1L]
						distances_stack = distances_packet[0L:packet_length-1L]
					end else begin
						indices_stack = [ [indices_stack], [indices_packet[*,0L:packet_length-1L]] ]
						coords_stack = [ [coords_stack], [coords_packet[*,0L:packet_length-1L]] ]
						distances_stack = [ distances_stack, distances_packet[0L:packet_length-1L] ]
					end
					packet++
					print, "Packet: ", packet, pos
				end
				if (num_points[pos] gt max_packet_length) then begin
					max_packet_length = num_points[pos] * 10L
					indices_packet = dblarr (3, max_packet_length)
					coords_packet = dblarr (3, max_packet_length)
					distances_packet = dblarr (max_packet_length)
				end
				packet_length = 0L
			end
			indices_packet[*,packet_length:packet_length+num_points[pos]-1L] = stream.indices
			coords_packet[*,packet_length:packet_length+num_points[pos]-1L] = stream.coords
			distances_packet[packet_length:packet_length+num_points[pos]-1L] = stream.distances
			origin[pos] = stream.origin
			length[pos] = stream.length
			first[pos] = line_pos
			line_pos += num_points[pos]
			last[pos] = line_pos - 1L
			packet_length += num_points[pos]
		end
		if (packet_length gt 0L) then begin
			if (packet eq 0L) then begin
				indices_stack = indices_packet[*,0L:packet_length-1L]
				coords_stack = coords_packet[*,0L:packet_length-1L]
				distances_stack = distances_packet[0L:packet_length-1L]
			end else begin
				indices_stack = [ [indices_stack], [indices_packet[*,0L:packet_length-1L]] ]
				coords_stack = [ [coords_stack], [coords_packet[*,0L:packet_length-1L]] ]
				distances_stack = [ distances_stack, distances_packet[0L:packet_length-1L] ]
			end
			packet++
			print, "Packet: ", packet, pos
		end
		if (keyword_set (return_indices)) then return, indices_stack
		return, { indices:indices_stack, coords:coords_stack, distances:distances_stack, num_points:num_points, length:length, origin:origin, num_lines:num_lines, first:first, last:last }
	end

	if (dir eq 0) then begin
		; Combine forward and backward streamlines from starting point
		along = pc_get_streamline (data, anchor=anchor, grid=grid, distances=distances, coords=coords, direction=1, periodic=periodic, precision=precision, select=select, length=length, num_points=num_points, origin=origin, max_length=max_length, /return_indices)
		against = pc_get_streamline (data, anchor=anchor, grid=grid, distances=d2, coords=against_coords, direction=-1, periodic=periodic, precision=precision, select=select, length=l2, num_points=n2, max_length=max_length, /return_indices)
		if (n2 le 1) then begin
			if (keyword_set (return_indices)) then return, along
			return, { indices:along, coords:coords, distances:distances, num_points:num_points, length:length, origin:origin, num_lines:1L, first:[ 0L ], last:[ num_points-1L ] }
		end
		length += l2
		num_points += n2 - 1
		origin = n2 - 1
		distances = [ -reverse (d2[1:*]), distances ]
		against = against[*,1:*]
		against_coords = against_coords[*,1:*]
		if (size (against, /n_dim) eq 2) then against = reverse (against, 2)
		if (size (against_coords, /n_dim) eq 2) then against_coords = reverse (against_coords, 2)
		coords = [ [against_coords], [coords] ]
		indices = [ [against], [along] ]
		if (keyword_set (return_indices)) then return, indices
		return, { indices:indices, coords:coords, distances:distances, num_points:num_points, length:length, origin:origin, num_lines:1L, first:[ 0L ], last:[ num_points-1L ] }
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
		if (nx+2*nghost ne n_elements (x)) then message, "ERROR: the data doesn't fit to the X-grid coordinates."
		if (ny+2*nghost ne n_elements (y)) then message, "ERROR: the data doesn't fit to the Y-grid coordinates."
		if (nz+2*nghost ne n_elements (z)) then message, "ERROR: the data doesn't fit to the Z-grid coordinates."
		; Convert anchor point into equidistant unit grid coordinates
		pos[0] = pc_find_index (anchor[0], x, num=mx) - nghost
		pos[1] = pc_find_index (anchor[1], y, num=my) - nghost
		pos[2] = pc_find_index (anchor[2], z, num=mz) - nghost
	end

	; Iterate finding points on the streamline
	origin = 0
	num_points = 1L
	length = 0.0d0
	indices = [ [pos] ]
	coords = [ [anchor] ]
	distances = [ length ]
	last = anchor
	num_total = 0L
	done = 0
	local_data = data[0:1,0:1,0:1,*]
	while (all (((pos ge 0) and (pos le ([nx, ny, nz]-1))) or periodic) and (length lt max_length) and not done) do begin

		; Interpolate data
		int_pos = (floor (pos) < (Box_xyz_upper-1)) > Box_xyz_lower
		residual = pos - int_pos
		lower = int_pos
		upper = int_pos + 1
		if (periodic[0] and (upper[0] ge nx)) then upper[0] = 0
		if (periodic[1] and (upper[1] ge ny)) then upper[1] = 0
		if (periodic[2] and (upper[2] ge nz)) then upper[2] = 0
		local_data[0,0,0,*] = data[lower[0],lower[1],lower[2],*]
		local_data[0,0,1,*] = data[lower[0],lower[1],upper[2],*]
		local_data[0,1,0,*] = data[lower[0],upper[1],lower[2],*]
		local_data[0,1,1,*] = data[lower[0],upper[1],upper[2],*]
		local_data[1,0,0,*] = data[upper[0],lower[1],lower[2],*]
		local_data[1,0,1,*] = data[upper[0],lower[1],upper[2],*]
		local_data[1,1,0,*] = data[upper[0],upper[1],lower[2],*]
		local_data[1,1,1,*] = data[upper[0],upper[1],upper[2],*]
		vector_x = interpolate (local_data[*,*,*,0], residual[0], residual[1], residual[2])
		vector_y = interpolate (local_data[*,*,*,1], residual[0], residual[1], residual[2])
		vector_z = interpolate (local_data[*,*,*,2], residual[0], residual[1], residual[2])
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

		; Find new position in grid coordinates
		point = pos
		if (keyword_set (grid)) then begin
			point[0] = interpolate (x, pos[0] + nghost)
			point[1] = interpolate (y, pos[1] + nghost)
			point[2] = interpolate (z, pos[2] + nghost)
		end
		last = point

		; Skip unselected points
		if (((num_total++ mod select) ne 0) and not done) then continue

		; Add indices and grid coordinates to the list of traced streamline points
		distances = [ distances, length ]
		indices = [ [indices], [pos] ]
		coords = [ [coords], [point] ]
		num_points++
	end

	if (keyword_set (return_indices)) then return, indices
	return, { indices:indices, coords:coords, distances:distances, num_points:num_points, length:length, origin:origin, num_lines:1L, first:[ 0L ], last:[ num_points-1L ] }

end

