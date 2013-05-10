;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_find_streamline.pro     ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  $Id$
;
;  Description:
;   Find nearest streamline point by a given coordinate.
;
;  Parameters:
;   * coord          Coordinate.
;   * streamlines    Streamlines structure containing the traced streamlines.
;   * nearest        Array index of the point nearest to the given coordinate.
;
;  Returns:
;   * num            Structure index number of the nearest streamline.
;
;  Example:
;  ========
;
;   Load varfile and extract Temperature along a magnetic filedline:
;   IDL> pc_read_var_raw, obj=var, tags=tags, grid=grid
;   IDL> B = pc_get_quantity ('B', var, tags)
;   IDL> Temp = pc_get_quantity ('Temp', var, tags)
;   IDL> indices = pc_get_streamline (B, anchor=[2.0, 3.5, 1.2], grid=grid, coords=coords, distances=distances, length=length)
;   IDL> streamlines = { num:1, streamline_1:{ indices:indices, coords:coords } }
;   IDL> nearest_num = pc_find_streamline ([1.9, 3.6, 1.2], streamlines, nearest=nearest_point)
;


; Search for the nearest streamline point.
function pc_find_streamline, coord, streamlines, nearest=nearest, distance=distance

	; Iterate over streamline sets
	nearest = lonarr (streamlines.num.lines)
	distances_stream = dblarr (streamlines.num.lines)
	num_comp = n_elements (coord)
	for set = 1L, streamlines.num.sets do begin
		; Find nearest point within streamline set
		distances = streamlines.(set).coords
		for comp = 0, num_comp - 1 do distances[comp,*] -= coord[comp]
		distances = total (distances^2, 1)
		distances_stream[set-1L] = min (distances)
		nearest[set-1L] = where (distances eq distances_stream[set-1L])
	end

	; Find nearest streamline
	distance = min (distances_stream)
	num = (where (distances_stream eq distance))[0]
	nearest = nearest[num]
	distance = sqrt (distance)

	return, (num + 1L)

end

