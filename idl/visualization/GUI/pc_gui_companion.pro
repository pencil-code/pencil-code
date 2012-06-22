;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_gui_companion.pro     ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  $Id$
;;;
;;;  Description:
;;;   Framework for precalculation and comparision of output in pencil units.
;;;   Companion procedures needed by 'pc_gui.pro'.
;;;
;;;  To do:
;;;   Add more comments


; Prepares the varset
pro pc_gui_prepare_varset, num, units, coords, varset, overset, dir, params, run_params, idlvar_list

	common varset_common, set, overplot, oversets, unit, coord, varsets, varfiles, datadir, sources, param, run_param, var_list

	datadir = dir

	unit = units
	coord = coords
	param = params
	run_param = run_params
	var_list = idlvar_list

	varfiles = { title:"-", time:0.0d0, loaded:0, number:-1, precalc_done:0 }
	varfiles = replicate (varfiles, num)

	varsets = replicate (varset, num)
	oversets = replicate (overset, num)
end


; Precalculates a data set and loads data, if necessary
pro pc_gui_precalc, i, number=number, varfile=varfile, datadir=dir, dim=dim, param=par, run_param=run_par, varcontent=varcontent, allprocs=allprocs, show_aver=show_aver, time=time, cut_x=cut_x, cut_y=cut_y, cut_z=cut_z

	common varset_common, set, overplot, oversets, unit, coord, varsets, varfiles, datadir, sources, param, run_param, var_list

	; Default settings
	default, show_aver, 0
	default, cut_x, -1
	default, cut_y, -1
	default, cut_z, -1
	default, number, i
	default, dir, pc_get_datadir()
	default, datadir, dir
	default, time, 0.0d0
	if (keyword_set (par)) then param = par
	if (keyword_set (run_par)) then run_param = run_par

	if (varfiles[i].number le 0) then varfiles[i].number = number

	if (varfiles[i].loaded eq 0) then begin
		default, varfile, "var.dat"
		if (n_elements (vars) eq 0) then begin
			print, 'Reading: ', varfile, ' ... please wait!'
			if (total([cut_x, cut_y, cut_z] < 0) ge -2) then begin
				pc_read_slice_raw, varfile=varfile, var_list=var_list, object=vars, tags=tags, datadir=datadir, param=param, par2=run_param, varcontent=varcontent, time=time, quiet=(i ne 0), cut_x=cut_x, cut_y=cut_y, cut_z=cut_z
			end else begin
				pc_read_var_raw, varfile=varfile, object=vars, tags=tags, datadir=datadir, dim=dim, param=param, par2=run_param, varcontent=varcontent, allprocs=allprocs, time=time, quiet=(i ne 0)
			end
			sources = varcontent.idlvar
			sources = sources[where (varcontent.idlvar ne 'dummy')]
			pc_gui_precalc_data, number, vars, tags, dim, grid
			print, 'Ready.'
		end
		varfiles[i].title = varfile
		varfiles[i].loaded = 1
		varfiles[i].precalc_done = 1
		varfiles[i].time = time * unit.time / unit.default_time
		vars = 0
	end

	if (show_aver) then draw_averages, number
	if (keyword_set (par)) then par = param
	if (keyword_set (run_par)) then run_par = run_param
end


; Precalculates a data set
pro pc_gui_precalc_data, i, vars, index, dim, gird

	common varset_common, set, overplot, oversets, unit, coord, varsets, varfiles, datadir, sources, param, run_param, var_list

	; First and last physical value, excluding ghost cells
	l1 = coord.l1
	l2 = coord.l2
	m1 = coord.m1
	m2 = coord.m2
	n1 = coord.n1
	n2 = coord.n2

	; Compute all desired quantities from available source data
	tags = tag_names (varsets[i])
	num = n_elements (tags)
	for pos = 0, num-1 do begin
		tag = tags[pos]
		last = (pos eq num-1)
		varsets[i].(pos) = pc_get_quantity (vars, index, tag, unit=unit, dim=dim, grid=grid, param=param, run_param=run_param, datadir=datadir, /cache, clean=last)

		; Divide by default units, where applicable.
		if (any (strcmp (tag, ['u_abs', 'u_x', 'u_y', 'u_z'], /fold_case))) then $
			varsets[i].(pos) /= unit.default_velocity
		if (any (strcmp (tag, ['Temp'], /fold_case))) then $
			varsets[i].(pos) /= unit.default_temperature
		if (any (strcmp (tag, ['rho'], /fold_case))) then $
			varsets[i].(pos) /= unit.default_density
		if (any (strcmp (tag, ['ln_rho'], /fold_case))) then $
			varsets[i].(pos) -= alog (unit.default_density)
		if (any (strcmp (tag, ['log_rho'], /fold_case))) then $
			varsets[i].(pos) -= alog10 (unit.default_density)
		if (any (strcmp (tag, ['rho_u_z'], /fold_case))) then $
			varsets[i].(pos) /= unit.default_density * unit.default_velocity
		if (any (strcmp (tag, ['b_x', 'b_y', 'b_z'], /fold_case))) then $
			varsets[i].(pos) /= unit.default_magnetic_field
		if (any (strcmp (tag, ['Spitzer_dt'], /fold_case))) then $
			varsets[i].(pos) /= unit.default_time
		if (any (strcmp (tag, ['j_abs'], /fold_case))) then $
			varsets[i].(pos) /= unit.default_current_density
	end

	; Compute all desired overplot quantities from available source data
	tags = tag_names (oversets[i])
	num = n_elements (tags)
	for pos = 0, num-1 do begin
		last = (pos eq num-1)
		oversets[i].(pos) = float (pc_get_quantity (vars, index, tags[pos], unit=unit, dim=dim, grid=grid, param=param, run_param=run_param, datadir=datadir, /cache, clean=last))

		; Divide by default units, where applicable.
		if (any (strcmp (tag, ['uu'], /fold_case))) then $
			oversets[i].(pos) /= float (unit.default_velocity)
		if (any (strcmp (tag, ['b'], /fold_case))) then $
			oversets[i].(pos) /= float (unit.default_magnetic_field)
	end
end


; Dummy routine
pro pc_gui_companion

	pc_gui_companion_loaded = 1
end

