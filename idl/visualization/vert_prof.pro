;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   vert_prof.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  $Id$
;;;
;;;  Description:
;;;    Plots a vertical profile of the given 3D quantity.
;;;

; Calculate and draw a vertical profile of the given 3D data
;
; data:   3D data cube (can be including ghost cells)
; coord:  coordinates for the vertical position of data in the cube,
;         asumed to be in the center of the data cube (eg. without ghost cells).
;         If omitted, the index numbers are used as coordinates.
; title:  title string for the plot
; log:    set this to use a logarithmic scale for data display
;
pro vert_prof, data, coord=coord, title=title, log=log

	num_z = (size (data))[3]
	if (n_elements (coord) eq 0) then coord = findgen (num_z)
	num_coord_z = (size (coord))[3]
	start_z = (num_z - num_coord_z) / 2

	prof_mean = dblarr (num_z)
	prof_min = dblarr (num_z)
	prof_max = dblarr (num_z)

	for z = 0, num_coord_z - 1 do begin
		prof_mean[z] = mean (data[*,*,start_z + z], /double)
		tmp = minmax (data[*,*,start_z + z])
		prof_min[z] = tmp[0]
		prof_max[z] = tmp[1]
	end
	data_min = min (data)
	data_max = max (data)

	plot, prof_mean, coord, xlog=log, xs=2, ys=1, xrange=[data_min, data_max], yrange=[2*coord[0]-coord[1], 2*coord[num_coord_z-1]-coord[num_coord_z-2]], title=title
	oplot, prof_mean, coord, psym=3, color=200
	oplot, prof_min, coord, linestyle=2
	oplot, prof_max, coord, linestyle=2
end

