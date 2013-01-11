;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   cmp_cslice.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  $Id$
;;;
;;;  Description:
;;;   Fast and powerful to use tool to view and compare slices of 3D data
;;;  To do:
;;;   Add more comments


; Simple interface to cmp_cslice_cache without caching mechanism
pro cmp_cslice, sets, limits=limits, units=units, coords=coords, scaling=scaling

	resolve_routine, "pc_gui_companion", /COMPILE_FULL_FILE, /NO_RECOMPILE

	; DEFAULT SETTINGS:
	min_size = 128

	; setup cube dimensions
	dims = size (sets.(0))
	dimensionality = dims[0]
	num_x = dims[1]
	if (dimensionality ge 2) then num_y = dims[2] else num_y = 1
	if (dimensionality ge 3) then num_z = dims[3] else num_z = 1

	set_names = tag_names (sets)
	num_names = n_elements (set_names)

	exe_1 = "varsets = { "
	exe_2 = "set = { "
	for i = 0, num_names-1 do begin
		exe_1 += set_names[i]+":reform(sets."+set_names[i]+","+strtrim(num_x,2)+","+strtrim(num_y,2)+","+strtrim(num_z,2)+"), "
		exe_2 += set_names[i]+":'"+set_names[i]+"', "
	end
	exe_1 = strmid (exe_1, 0, strlen (exe_1)-2)
	exe_2 = strmid (exe_2, 0, strlen (exe_2)-2)
	exe_1 += " }"
	exe_2 += " }"

	res = execute (exe_1)
	if (not res) then begin
		print, "ERROR: could not generate varsets!"
		print, "("+exe_1+")"
		return
	end

	res = execute (exe_2)
	if (not res) then begin
		print, "ERROR: could not generate set!"
		print, "("+exe_2+")"
		return
	end

	varfiles = { title:"N/A", time:0.0d0, loaded:1, number:0, precalc_done:1 }

	; setup coordinates, if necessary
	if (n_elements (coords) eq 0) then begin
		coords = { x:findgen(num_x), y:findgen(num_y), z:findgen(num_z), dx:1.0, dy:1.0, dz:1.0, dx_1:replicate(1.0, num_x), dy_1:replicate(1.0, num_y), dz_1:replicate(1.0, num_z), nx:num_x, ny:num_y, nz:num_z, orig_nx:num_x, orig_ny:num_y, orig_nz:num_z, x_off:0, y_off:0, z_off:0, l1:0, l2:num_x-1, m1:0, m2:num_y-1, n1:0, n2:num_z-1, lequidist:[1,1,1], lperi:[0,0,0], ldegenerated:[0,0,0], nghost:0 }
	end

	; setup a scaling factor to have a minimum size, if necessary
	default, scaling, 1
	dims = (size (varsets.(0)))[1:size (varsets.(0), /n_dimensions)] * scaling
	if (not any (dims ge min_size)) then begin
		scaling *= ceil (min_size / double (max (dims)))
	end
	if (n_elements (scaling) eq 1) then scaling = [ scaling, scaling, scaling ]

	cmp_cslice_cache, set, set_content=varsets, set_files=varfiles, limits=limits, units=units, coords=coords, scaling=scaling

	return
end

