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

	common varset_common, set, overplot, oversets, unit, coord, varsets, varfiles, datadir, sources, param, run_param

	resolve_routine, "pc_gui_companion", /COMPILE_FULL_FILE, /NO_RECOMPILE

	; DEFAULT SETTINGS:
	min_size = 128

	set_names = tag_names (sets)
	num_names = n_elements (set_names)

	exe_1 = "varsets = { "
	exe_2 = "set = { "
	for i=0, num_names-1 do begin
		if ((size (sets.(i)))[0] ne 3) then continue
		exe_1 += set_names[i]+":sets."+set_names[i]+", "
		exe_2 += set_names[i]+":'"+set_names[i]+"', "
	end
	exe_1 = strmid (exe_1, 0, strlen (exe_1)-2)
	exe_2 = strmid (exe_2, 0, strlen (exe_2)-2)
	exe_1 += " }"
	exe_2 += " }"

	res = execute (exe_1)
	if (not res) then begin
		print, "ERROR: could not generate varsets!"
		return
	end

	res = execute (exe_2)
	if (not res) then begin
		print, "ERROR: could not generate set!"
		return
	end

	varfiles = { title:"N/A", loaded:1, number:0, precalc_done:1 }

	; setup a scaling factor to have a minimum size, if necessary
	default, scaling, 1
	dims = (size (varsets.(0)))[1:size (varsets.(0), /n_dimensions)]
	if (not any (dims ge min_size)) then begin
		scaling = ceil (min_size / double (max (dims)))
	end
	if (n_elements (scaling) eq 1) then scaling = [ scaling, scaling, scaling ]

	; setup limits, if necessary
	if (n_elements (limits) eq 0) then begin
		dims = size (varsets.(0))
		num_x = dims[1]
		num_y = dims[2]
		num_z = dims[3]
		limits = reform (lindgen (dims[1], dims[2], dims[3]), num_x, num_y, num_z)
	end

	cmp_cslice_cache, set, limits, units=units, coords=coords, scaling=scaling

	return
end

