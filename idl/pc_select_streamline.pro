;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_select_streamline.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  $Id$
;
;  Description:
;   Select a streamline by number from a given streamlines structure.
;
;  Parameters:
;   * streamlines    Streamlines structure.
;   * num            Number of streamline to return.
;
;  Returns:
;   * streamline     Returns the num-th streamline of the given streamlines structure.
;
;  Example:
;  ========
;
;   Load varfile and extract Temperature along the second magnetic filedline:
;   IDL> pc_read_var_raw, obj=var, tags=tags, grid=grid
;   IDL> B = pc_get_quantity ('B', var, tags)
;   IDL> Temp = pc_get_quantity ('Temp', var, tags)
;   IDL> seeds = pc_seed_points (grid)
;   IDL> streamlines = pc_get_streamline (B, anchor=seeds, grid=grid)
;   IDL> Temp_streamlines = pc_extract_streamline (Temp, streamlines, name='Temp')
;   IDL> stream_2 = pc_select_streamline (streamlines, 2)
;   IDL> Temp_2 = pc_select_streamline (Temp_streamlines, 2, streamlines=streamlines)
;   IDL> plot, stream_2.distances, Temp_2, xtitle="coordinate along streamline", ytitle="temperature"
;


; Find streamline in a structure of sets of streamlines.
function pc_select_streamline, data, num, streamlines=streamlines

	num_sets = data.num
	if (num_sets lt 1L) then message, "ERROR: The given streamlines structure doesn't contain any streamlines."

	is_quantity = 0
	if (size (data.(1), /type) ne 8) then begin
		if (not keyword_set (streamlines)) then message, "ERROR: For selection of a quantity data streamline, the corresponding streamlines structure must begiven as parameter."
		; This is a quantity structure, extract the needed data values
		is_quantity = 1
		num_components = size (actual_set, /n_dim)
	end

	; Iterate though sets of streamlines
	passed_num = 0L
	for set = 1L, num_sets do begin
		if (keyword_set (streamlines)) then actual_set = streamlines.(set) else actual_set = data.(set)
		num_lines = actual_set.num_lines
		if (num le (passed_num + num_lines)) then begin
			if (num_lines eq 1) then begin
				; Streamline is within actual set
				if (is_quantity) then return, data.(set) else return, actual_set
			end
			; Need to extract data from a set of multiple streamlines
			pos = num - passed_num - 1L
			if (is_quantity) then begin
				if (num_components eq 1) then return, (data.(pos+1L))[first[pos]:last[pos]] else return, (data.(pos+1L))[*,first[pos]:last[pos]]
			end
			indices = actual_set.indices[*,actual_set.first[pos]:actual_set.last[pos]]
			coords = actual_set.coords[*,actual_set.first[pos]:actual_set.last[pos]]
			distances = actual_set.distances[actual_set.first[pos]:actual_set.last[pos]]
			num_points = actual_set.num_points[pos]
			length = actual_set.length[pos]
			origin = actual_set.origin[pos]
			return, { indices:indices, coords:coords, distances:distances, num_points:num_points, length:length, origin:origin, num_lines:1L, first:[ 0L ], last:[ num_points-1L ] }
		end
		passed_num += num_lines
	end

	message, "ERROR: The given streamlines structure has only "+strtrim (passed_num, 2)+" elements."
	return, !Values.D_NaN

end

