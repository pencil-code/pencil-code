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

	; Iterate over streamlines
	nearest = intarr (streamlines.num)
	distances_stream = dblarr (streamlines.num)
	for pos = 1, streamlines.num do begin
		; Find nearest point within streamline
		distances = streamlines.(pos).coords
		for i = 0, streamlines.(pos).num - 1 do distances[*,i] -= coord
		distances = total (distances^2, 1)
		distances_stream[pos-1] = min (distances)
		nearest[pos-1] = where (distances eq distances_stream[pos-1])
	end

	; Find nearest streamline
	distance = min (distances_stream)
	num = (where (distances_stream eq distance))[0]
	nearest = nearest[num]
	distance = sqrt (distance)

	return, (num + 1)

end

