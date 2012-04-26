;;;;;;;;;;;;;;;;;;;;;;
;;;   cslice.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;
;;;  $Id$
;;;
;;;  Description:
;;;   Fast and powerful to use tool to view and compare slices of 3D data
;;;  To do:
;;;   Add more comments


; Compatibility mode => calls simple interface
pro cslice, cube, limits, units=units, coords=coords, scaling=scaling

	resolve_routine, "pc_gui_companion", /COMPILE_FULL_FILE, /NO_RECOMPILE

	; some error checking
	if (n_elements (cube) le 0) then begin
		print, "Sorry, the data is undefined."
		return
	end

	if (size (cube, /n_dimensions) eq 4) then begin
		n_dims = (size (cube))[4]
		dims_str = "cube_1:cube[*,*,*,0]"
		if (n_dims gt 1) then begin
			for dim=1, n_dims-1 do begin
				dims_str += ", cube_"+strtrim (dim+1,2)+":cube[*,*,*,"+strtrim (dim,2)+"]"
			end
		end
		res = execute ("set = { "+dims_str+" }")
		if (not res) then begin
			print, "Could not create 4-dimensional dataset!"
			stop
		end
	end else begin
		set = { cube:cube }
	end
	cmp_cslice, set, limits=limits, units=units, coords=coords, scaling=scaling

	return
end

