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

	resolve_routine, "analyse_companion", /COMPILE_FULL_FILE, /NO_RECOMPILE

	; scaling factor for visualisation
	dim = size (cube)
	default, scaling, fix (128 / max (dim[1:3]))
	if (n_elements (scaling) eq 1) then if (scaling lt 1) then scaling = 1

	; setup limits, if necessary
	if (n_elements (limits) eq 0) then limits = lindgen ((size (cube))[1], (size (cube))[2], (size (cube))[3])

	set = { cube:cube }
	cmp_cslice, set, limits, units=units, scaling=scaling

	return
end

