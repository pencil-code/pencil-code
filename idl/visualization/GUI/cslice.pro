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

	set = { cube:cube }
	cmp_cslice, set, limits=limits, units=units, coords=coords, scaling=scaling

	return
end

