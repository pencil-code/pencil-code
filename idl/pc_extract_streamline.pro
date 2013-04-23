;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_extract_streamline.pro     ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  $Id$
;
;  Description:
;   Extraction of any quantiy along previously traced streamlines.
;
;  Parameters:
;   * data           Data cube of a scalar field or vector field (3- or 4-dimensional).
;   * indices        Data indices of the traced streamline returned by 'pc_get_streamline'.
;   * name           Name of the quantity to be extraced.
;   * label          Tag-name label of the quantity inside the returned structure.
;   * packet_size    Size of streamlines packet to be used, if a streamlines structure is given.
;   * return_values  Return the extracted data values as array instead of inside a structure.
;
;  Returns:
;   * quantity       Extracted (interpolated) data along the given streamlines.
;
;  Examples:
;  =========
;
;   Load varfile and extract Temperature along a magnetic filedline:
;   IDL> pc_read_var_raw, obj=var, tags=tags, grid=grid
;   IDL> B = pc_get_quantity ('B', var, tags)
;   IDL> Temp = pc_get_quantity ('Temp', var, tags)
;   IDL> indices = pc_get_streamline (B, anchor=[2.0, 3.5, 1.2], grid=grid, distances=distances, length=length, /return_indices)
;   IDL> Temp_streamline = pc_extract_streamline (Temp, indices, /return_values)
;   IDL> B_streamline = pc_extract_streamline (B, indices, /return_values)
;   IDL> B_abs = B_streamline[0,*]^2 + B_streamline[1,*]^2 + B_streamline[2,*]^2
;   IDL> plot, distances, Temp_streamline, xtitle="coordinate along streamline", ytitle="temperature"
;   IDL> plot, distances, B_abs, xtitle="coordinate along streamline", ytitle="magnetic field", /ylog
;
;   Load varfile and extract Temperature along several magnetic filedlines, using caching:
;   IDL> pc_read_var_raw, obj=var, tags=tags, grid=grid
;   IDL> B = pc_get_quantity ('B', var, tags)
;   IDL> Temp = pc_get_quantity ('Temp', var, tags)
;   IDL> streamline_1 = pc_get_streamline (B, anchor=[2.0, 3.5, 1.2], grid=grid, /cache)
;   IDL> streamline_2 = pc_get_streamline (anchor=[2.2, 3.3, 1.1], grid=grid, /cache)
;   IDL> streamline_3 = pc_get_streamline (anchor=[2.1, 3.4, 1.0], grid=grid)
;   IDL> streamlines = { num:3, set_1:streamline_1, set_2:streamline_2, set_3:streamline_3 }
;   IDL> Temp_streamlines = pc_extract_streamline (Temp, streamlines, name='Temp')
;   IDL> B_streamlines = pc_extract_streamline (B, streamlines, name='B')
;   IDL> Temp_streamline = pc_select_streamline (Temp_streamline, 1)
;   IDL> B_streamline = pc_select_streamline (B_streamline, 1)
;   IDL> B_abs = B_streamline[0,*]^2 + B_streamline[1,*]^2 + B_streamline[2,*]^2
;   IDL> plot, distances, Temp_streamline, xtitle="coordinate along streamline", ytitle="temperature"
;   IDL> plot, distances, B_abs, xtitle="coordinate along streamline", ytitle="magnetic field", /ylog
;


; Calculation of streamline coordinates.
function pc_extract_streamline, data, streamlines, name=name, label=label, precision=precision, packet_size=packet_size, return_values=return_values

	; Default settings:
	default_name = 'Quantity'
	if (not keyword_set (packet_size)) then packet_size = 100L
	if (not keyword_set (precision)) then precision = 'D'
	if (precision ne 'F') then precision = 'D'

	if (n_elements (data) eq 0) then message, "ERROR: no data array given."
	if (n_elements (streamlines) eq 0) then message, "ERROR: no streamline(s) given."
	if (size (streamlines, /type) ne 8) then begin
		num_points = long (n_elements (streamlines[0,*]))
		streamlines = { num:1L, set_1:{ indices:streamlines, num_points:num_points, num_lines:1L, first:[ 0L ], last:[ num_points-1L ] } }
	end
	if (not any (strcmp (tag_names (streamlines.(1)), 'indices', /fold_case))) then message, "ERROR: no indices in given streamlines structure."
	if (not keyword_set (name)) then name = default_name
	if (size (name, /type) ne 7) then name = default_name
	name = strtrim (name, 2)
	if (not keyword_set (label)) then label = name
	if (size (label, /type) ne 7) then label = name
	label = idl_validname (label, /convert_all)

	n_dim = size (data, /n_dim)
	if ((n_dim lt 3) or (n_dim gt 4)) then message, "ERROR: data array dimension is invalid."

	s = size (data)
	nx = s[1]
	ny = s[2]
	nz = s[3]
	if (n_dim eq 3) then num = 1 else num = s[4]

	; Iterate over streamlines
	quantity = { name:name }
	for set = 1L, streamlines.num do begin

		indices = streamlines.(set).indices
		num_points = total (streamlines.(set).num_points)
		if (precision eq 'D') then extract = dblarr (num, num_points) else extract = fltarr (num, num_points)

		; Follow the streamline
		for pos = 0L, num_points - 1L do begin
			int_pos = (floor (indices[*,pos]) < ([nx, ny, nz] - 2)) > 0
			residual = indices[*,pos] - int_pos
			loc_data = data[int_pos[0]:int_pos[0]+1,int_pos[1]:int_pos[1]+1,int_pos[2]:int_pos[2]+1,*]
			; Iterate over the data components
			for comp = 0, num - 1 do begin
				extract[comp,pos] = interpolate (loc_data[*,*,*,comp], residual[0], residual[1], residual[2])
			end
		end
		quantity = create_struct (quantity, label+'_'+strtrim (set, 2), reform (extract))
	end

	if (keyword_set (return_values)) then begin
		; Return data values in an array
		num_points = 0L
		for set = 1L, streamlines.num do num_points += total (streamlines.(set).num_points)
		if (precision eq 'D') then quantity_stack = dblarr (num, num_points) else quantity_stack = fltarr (num, num_points)
		start = 0L
		for pos = 1L, streamlines.num do begin
			actual_set = streamlines.(pos)
			num_transfer = actual_set.num_points
			if (num eq 1) then begin
				quantity_stack[start:start+num_transfer-1L] = quantity.(pos)
			end else begin
				quantity_stack[*,start:start+num_transfer-1L] = quantity.(pos)
			end
			start += num_transfer
		end
		return, quantity_stack
	end

	; Return data values set structure
	return, quantity

end

