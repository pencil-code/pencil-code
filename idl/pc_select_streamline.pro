;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_select_streamline.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  $Id$
;
;  Description:
;   Select a streamline by number from a given streamlines structure.
;
;  Parameters:
;   * data           Streamlines structure or quantity data array.
;   * num            Number of streamline to select and return.
;   * streamlines    Streamlines structure, only needed if "data" is a quantity data array.
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

	if (data.num.sets lt 1L) then message, "ERROR: The given streamlines structure doesn't contain any streamlines."

	; Iterate though sets of streamlines
	passed_num = 0L
	for set = 1L, data.num.sets do begin
		if (keyword_set (streamlines)) then num_lines = streamlines.(set).num_lines else num_lines = data.(set).num_lines
		if (num le (passed_num + num_lines)) then begin
			; Streamline lies within actual set
			if (num_lines eq 1) then return, data.(set)
			; Need to extract data from a set of multiple streamlines
			pos = num - passed_num - 1L
			if (keyword_set (streamlines)) then begin
				first = streamlines.(set).first
				last = streamlines.(set).last
				; Return selected part of the quantity streamline
				if (size (data.(1), /n_dimensions) eq 1) then begin
					return, (data.(pos+1L))[first[pos]:last[pos]]
				end else begin
					return, (data.(pos+1L))[*,first[pos]:last[pos]]
				end
			end
			; Return selected part of the streamline
			first = data.(set).first
			last = data.(set).last
			indices = data.(set).indices[*,first[pos]:last[pos]]
			coords = data.(set).coords[*,first[pos]:last[pos]]
			distances = data.(set).distances[first[pos]:last[pos]]
			num_points = data.(set).num_points[pos]
			length = data.(set).length[pos]
			origin = data.(set).origin[pos]
			return, { indices:indices, coords:coords, distances:distances, num_points:num_points, length:length, origin:origin, num_lines:1L, first:[ 0L ], last:[ num_points-1L ] }
		end
		passed_num += num_lines
	end

	message, "ERROR: The given streamlines structure has only "+strtrim (passed_num, 2)+" elements."
	return, !Values.D_NaN

end

