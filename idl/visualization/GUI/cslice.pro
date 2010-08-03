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

	if (n_elements (limits) eq 0) then limits = lindgen ((size (cube))[1], (size (cube))[2], (size (cube))[3])

	set = { cube:cube }
	cmp_cslice, set, limits, units=units, scaling=scaling

	return
end

