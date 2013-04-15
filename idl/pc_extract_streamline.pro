;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_extract_streamline.pro     ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  $Id$
;
;  Description:
;   Extraction of any quantiy along a previously traced streamline.
;
;  Parameters:
;   * data           Data cube of a scalar field or vector field (3- or 4-dimensional).
;   * indices        Data indices of the traced streamline returned by 'pc_get_streamline'.
;   * name           Tag-name of quantity to be extraced, if a streamlines structure is given.
;   * packet_size    Size of streamlines packet to be used, if a streamlines structure is given.
;
;  Returns:
;   * quantity       Extracted (interpolated) data values along the traced streamline.
;
;  Examples:
;  =========
;
;   Load varfile and extract Temperature along a magnetic filedline:
;   IDL> pc_read_var_raw, obj=var, tags=tags, grid=grid
;   IDL> B = pc_get_quantity ('B', var, tags)
;   IDL> Temp = pc_get_quantity ('Temp', var, tags)
;   IDL> indices = pc_get_streamline (B, anchor=[2.0, 3.5, 1.2], grid=grid, distances=distances, length=length)
;   IDL> Temp_streamline = pc_extract_streamline (Temp, indices)
;   IDL> B_streamline = pc_extract_streamline (B, indices)
;   IDL> B_abs = B_streamline[0,*]^2 + B_streamline[1,*]^2 + B_streamline[2,*]^2
;   IDL> plot, distances, Temp_streamline, xtitle="coordinate along streamline", ytitle="temperature"
;   IDL> plot, distances, B_abs, xtitle="coordinate along streamline", ytitle="magnetic field", /ylog
;
;   Load varfile and extract Temperature along several magnetic filedlines, using caching:
;   IDL> pc_read_var_raw, obj=var, tags=tags, grid=grid
;   IDL> B = pc_get_quantity ('B', var, tags)
;   IDL> Temp = pc_get_quantity ('Temp', var, tags)
;   IDL> indices = pc_get_streamline (B, anchor=[2.0, 3.5, 1.2], grid=grid, coords=coords, num=num, origin=origin, /cache)
;   IDL> streamline_1 = { indices:indices, coords:coords, num:num, origin:origin }
;   IDL> indices = pc_get_streamline (anchor=[2.2, 3.3, 1.1], grid=grid, coords=coords, num=num, origin=origin, /cache)
;   IDL> streamline_2 = { indices:indices, coords:coords, num:num, origin:origin }
;   IDL> indices = pc_get_streamline (anchor=[2.1, 3.4, 1.0], grid=grid, coords=coords, num=num, origin=origin)
;   IDL> streamline_3 = { indices:indices, coords:coords, num:num, origin:origin }
;   IDL> streamlines = { num:3, streamline_1:streamline_1, streamline_2:streamline_2, streamline_3:streamline_3 }
;   IDL> Temp_streamlines = pc_extract_streamline (Temp, streamlines, name='Temp')
;


; Calculation of streamline coordinates.
function pc_extract_streamline, data, streamlines, name=name, packet_size=packet_size

	; Default settings:
	if (not keyword_set (packet_size)) then packet_size = 100L

	if (n_elements (data) eq 0) then message, "ERROR: no data array given."
	if (n_elements (streamlines) eq 0) then message, "ERROR: no streamline(s) given."
	if (size (streamlines, /type) ne 8) then streamlines = { num:1L, streamline_1:{ indices:streamlines, num:n_elements (streamlines[0,*]) } }
	if (not any (strcmp (tag_names (streamlines), 'indices', /fold_case)) then message, "ERROR: no indices in given streamlines structure."
	if (not keyword_set (name)) then name = 'QUANTITY'
	if (size (name, /type) ne 7) then name = 'QUANTITY'
	name = idl_validname (name, /convert_all)

	n_dim = size (data, /n_dim)
	if ((n_dim lt 3) or (n_dim gt 4)) then message, "ERROR: data array dimension is invalid."

	s = size (data)
	nx = s[1]
	ny = s[2]
	nz = s[3]
	if (n_dim eq 3) then num = 1 else num = s[4]

	; Iterate over streamlines
	quantity = { name:name }
	for line = 1L, streamlines.num do begin

		indices = streamlines.(line).indices
		length = streamlines.(line).num
		extract = dblarr (num, length)

		; Follow the streamline
		for pos = 0, length - 1 do begin
			int_pos = (floor (indices[*,pos]) < ([nx, ny, nz] - 2)) > 0
			residual = indices[*,pos] - int_pos
			loc_data = data[int_pos[0]:int_pos[0]+1,int_pos[1]:int_pos[1]+1,int_pos[2]:int_pos[2]+1,*]
			; Iterate over the data components
			for comp = 0, num - 1 do begin
				extract[comp,pos] = interpolate (loc_data[*,*,*,comp], residual[0], residual[1], residual[2])
			end
		end

		if ((line mod packet_size) eq 1L) then begin
			; Create streamlines packet
			packet = create_struct (name+'_'+strtrim (line, 2), reform (extract))
		end elseif (((line mod packet_size) eq 0L) or (line eq streamlines.num)) then begin
			; Transfer packet to streamlines structure
			quantity = create_struct (quantity, packet)
			packet = 0
		end else begin
			; Add streamline to packet
			packet = create_struct (packet, name+'_'+strtrim (line, 2), reform (extract))
		end
	end

	if (streamlines.num eq 1) then return, quantity.(1)
	return, quantity

end

