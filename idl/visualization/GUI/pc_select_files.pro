;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_select_files.pro      ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  $Id$
;;;
;;;  Description:
;;;   GUI for selecting a subset of available snapshot files.
;;;   The returned array contains a list of selected filenames.
;;;   Optional parameters are:
;;;   * num_selected (returns the number of selected files)
;;;   * pattern (contains the file-search pattern, default: "VAR[0-9]*")
;;;   * varfile (contains the varfile loaded by default, default: "var.dat" or "var.h5")
;;;   * addfile (contains an additional filename, default: "crash.dat" or "crash.h5")
;;;   * datadir (contains the datadir, default: pc_get_datadir)
;;;   * allprocs (contains the IO-strategy parameter, default: automatic)
;;;   * procdir (contains procdir based on the chosen IO-strategy)
;;;   * varcontent (contains or returns the varcontent structure)
;;;   * quantities (contains or returns the selected physical quantities)
;;;   * overplots (contains or returns the selected overplot quantities)
;;;   * cut_x (contains the pixel value in x of the yz-slice)
;;;   * cut_y (contains the pixel value in y of the xz-slice)
;;;   * cut_z (contains the pixel value in z of the xy-slice)
;;;   * xs (contains the starting pixel value in x of the sub-volume)
;;;   * ys (contains the starting pixel value in y of the sub-volume)
;;;   * zs (contains the starting pixel value in z of the sub-volume)
;;;   * xe (contains the ending pixel value in x of the sub-volume)
;;;   * ye (contains the ending pixel value in y of the sub-volume)
;;;   * ze (contains the ending pixel value in z of the sub-volume)
;;;   If an optional parameter is given as undefined, its default is returned.
;;;
;;;   Examples:
;;;   IDL> pc_select_files, obj=files, num_selected=num_selected
;;;   IDL> pc_select_files, obj=files, pattern="PVAR1[5-9]*"
;;;   IDL> pc_select_files, obj=files, varfile="myvar.dat", datadir=datadir
;;;   IDL> pc_select_files, obj=files, datadir="mydata", procdir=procdir
;;;   IDL> pc_select_files, obj=files, addfile="crash.dat", /allprocs
;;;   IDL> pc_select_files, obj=files, pattern="VAR[1-9]*", allprocs=allprocs
;;;   IDL> pc_select_files, obj=files, varcontent=varcontent, quantities=quantities
;;;   IDL> pc_select_files, obj=files, cut_x=cut_x, cut_y=cut_y, cut_z=cut_z
;;;


; Update a list of given quantities.
pro pc_select_files_update_list, list, all, indices, default=default, avail=avail_list

	common select_files_common, num_files, selected, num_selected, var_selected, add_selected, sources, sources_selected, num_cont, cont_selected, quant, quant_selected, quant_list, all_quant, quant_avail, over, over_selected, over_list, all_over, over_avail, cut_pos, max_pos, slice, skipping, stepping, data_dir, units, run_par, start_par, gb_per_file, cont_corr, subvol_corr, subvol_xs, subvol_ys, subvol_zs, subvol_xe, subvol_ye, subvol_ze, scaling_x, scaling_y, scaling_z, nx, ny, nz, nghost, skip_last, n_slice

	if (not keyword_set (default)) then default = all

	avail_list = "[N/A] ("+list+")"
	indices = -1
	if (any (cont_selected lt 0)) then return

	avail = pc_check_quantities (check=all, sources=sources[cont_selected], /indices)
	if (any (avail lt 0)) then return

	tags = tag_names (all)
	num = n_elements (avail)
	for pos = 0, num-1 do begin
		avail_list[avail[pos]] = all.(avail[pos])
		if (has_tag (default, tags[avail[pos]])) then indices = [ indices, avail[pos] ]
	end
	active = where (indices ge 0, num_active)
	if (num_active ge 1) then indices = indices[active]
end


; Event handling of file dialog window
pro select_files_event, event

	common select_files_gui_common, b_var, b_add, b_ts, c_list, i_skip, i_step, f_load, f_comp, scal_x, scal_y, scal_z, c_cont, c_quant, c_over, d_slice, cut_co, cut_sl, sub_xs, sub_ys, sub_zs, sub_xe, sub_ye, sub_ze, sub_nx, sub_ny, sub_nz, i_last, nf
	common select_files_common, num_files, selected, num_selected, var_selected, add_selected, sources, sources_selected, num_cont, cont_selected, quant, quant_selected, quant_list, all_quant, quant_avail, over, over_selected, over_list, all_over, over_avail, cut_pos, max_pos, slice, skipping, stepping, data_dir, units, run_par, start_par, gb_per_file, cont_corr, subvol_corr, subvol_xs, subvol_ys, subvol_zs, subvol_xe, subvol_ye, subvol_ze, scaling_x, scaling_y, scaling_z, nx, ny, nz, nghost, skip_last, n_slice

	WIDGET_CONTROL, WIDGET_INFO (event.top, /CHILD)
	WIDGET_CONTROL, event.id, GET_UVALUE = eventval

	quit = -1

	SWITCH eventval of
	'ALL': begin
		WIDGET_CONTROL, c_list, SET_LIST_SELECT = indgen (num_files)
		num_selected = num_files
		pc_select_files_update
		WIDGET_CONTROL, nf, SET_VALUE = strtrim(string(num_selected),2) + ' / ' + strtrim(string(num_files),2)
		break
	end
	'NONE': begin
		WIDGET_CONTROL, c_list, SET_LIST_SELECT = -1
		num_selected = 0
		pc_select_files_update
		WIDGET_CONTROL, nf, SET_VALUE = strtrim(string(num_selected),2) + ' / ' + strtrim(string(num_files),2)
		break
	end
	'SKIP': begin
		skipping = event.value
		if (skipping lt 0) then skipping = 0
		if (skipping ge num_files - skip_last) then skipping = num_files - skip_last - 1
		if (stepping gt num_files - skipping - skip_last) then begin
			stepping = num_files - skipping
			WIDGET_CONTROL, i_step, SET_VALUE = stepping
		end
		WIDGET_CONTROL, i_skip, SET_VALUE = skipping
		break
	end
	'STEP': begin
		stepping = event.value
		if (stepping lt 1) then stepping = 1
		if (stepping gt num_files - skipping - skip_last) then stepping = num_files - skipping - skip_last
		WIDGET_CONTROL, i_step, SET_VALUE = stepping
		break
	end
	'SKIP_LAST': begin
		skip_last = event.value
		if (skip_last lt 0) then skip_last = 0
		if (skip_last ge num_files - skipping) then skip_last = num_files - skipping - 1
		if (stepping gt num_files - skipping - skip_last) then begin
			stepping = num_files - skipping - skip_last
			WIDGET_CONTROL, i_step, SET_VALUE = stepping
		end
		WIDGET_CONTROL, i_last, SET_VALUE = skip_last
		break 
	end
	'APPLY': begin
		selected = indgen (num_files)
		selected = selected[skipping:num_files - 1 - skip_last:stepping]
		WIDGET_CONTROL, c_list, SET_LIST_SELECT = selected
		WIDGET_CONTROL, nf, SET_VALUE = strtrim(string(num_selected),2) + ' / ' + strtrim(string(num_files),2)
	end
	'LIST': begin
		selected = WIDGET_INFO (c_list, /LIST_SELECT)
		if (any (selected ne -1)) then num_selected = n_elements (selected) else num_selected = 0
		pc_select_files_update
		WIDGET_CONTROL, nf, SET_VALUE = strtrim(string(num_selected),2) + ' / ' + strtrim(string(num_files),2)
		break
	end
	'ADD': begin
		add_selected = event.select
		pc_select_files_update
		break
	end
	'VAR': begin
		var_selected = event.select
		pc_select_files_update
		break
	end
	'SLICE': begin
		if (slice ne event.index) then begin
			slice = event.index
			max_pos = -1
			case n_slice of
			2: begin
				if (slice eq 1) and (nx gt 1) then max_pos = nx - 1
				if (slice eq 1) and (ny gt 1) then max_pos = ny - 1
				if (slice eq 1) and (nz gt 1) then max_pos = nz - 1
			end
			4: begin
				if slice eq 1 and (nx GT 1) and (ny GT 1) then max_pos = nx - 1
				if slice eq 2 and (nx GT 1) and (ny GT 1) then max_pos = ny - 1
				if slice eq 1 and (nx GT 1) and (nz GT 1) then max_pos = nx - 1
				if slice eq 2 and (nx GT 1) and (nz GT 1) then max_pos = nz - 1
				if slice eq 1 and (ny GT 1) and (nz GT 1) then max_pos = ny - 1
				if slice eq 2 and (ny GT 1) and (nz GT 1) then max_pos = nz - 1
			end
			5: begin
				if (slice eq 1) then max_pos = nx-1
				if (slice eq 2) then max_pos = ny-1
				if (slice eq 3) then max_pos = nz-1
			end
			endcase
			if (max_pos lt 1) then cut_pos = -1 else cut_pos = max_pos/2
			WIDGET_CONTROL, cut_co, SET_VALUE = cut_pos>0
			WIDGET_CONTROL, cut_sl, SET_SLIDER_MIN = 0<max_pos
			WIDGET_CONTROL, cut_sl, SET_SLIDER_MAX = max_pos>1
			WIDGET_CONTROL, cut_sl, SET_VALUE = cut_pos
			WIDGET_CONTROL, cut_co, SENSITIVE = (cut_pos ge 0)
			WIDGET_CONTROL, cut_sl, SENSITIVE = (cut_pos ge 0)
			WIDGET_CONTROL, sub_xs, SENSITIVE = ((slice eq n_slice-1) and (nx gt 1)) 
			WIDGET_CONTROL, sub_xe, SENSITIVE = ((slice eq n_slice-1) and (nx gt 1))
			WIDGET_CONTROL, sub_nx, SENSITIVE = ((slice eq n_slice-1) and (nx gt 1))
			WIDGET_CONTROL, sub_ys, SENSITIVE = ((slice eq n_slice-1) and (ny gt 1))
			WIDGET_CONTROL, sub_ye, SENSITIVE = ((slice eq n_slice-1) and (ny gt 1))
			WIDGET_CONTROL, sub_ny, SENSITIVE = ((slice eq n_slice-1) and (ny gt 1))
			WIDGET_CONTROL, sub_zs, SENSITIVE = ((slice eq n_slice-1) and (nz gt 1))
			WIDGET_CONTROL, sub_ze, SENSITIVE = ((slice eq n_slice-1) and (nz gt 1))
			WIDGET_CONTROL, sub_nz, SENSITIVE = ((slice eq n_slice-1) and (nz gt 1))

			pc_select_files_update
		end
		break
	end
	'SCAL_X': begin
		WIDGET_CONTROL, scal_x, GET_VALUE = scaling_x
		break
	end
	'SCAL_Y': begin
		WIDGET_CONTROL, scal_y, GET_VALUE = scaling_y
		break
	end
	'SCAL_Z': begin
		WIDGET_CONTROL, scal_z, GET_VALUE = scaling_z
		break
	end
	'SCAL_PLUS': begin
		scaling_x *= 2
		scaling_y *= 2
		scaling_z *= 2
		WIDGET_CONTROL, scal_x, SET_VALUE = scaling_x
		WIDGET_CONTROL, scal_y, SET_VALUE = scaling_y
		WIDGET_CONTROL, scal_z, SET_VALUE = scaling_z
		break
	end
	'SCAL_MINUS': begin
		scaling_x /= 2
		scaling_y /= 2
		scaling_z /= 2
		WIDGET_CONTROL, scal_x, SET_VALUE = scaling_x
		WIDGET_CONTROL, scal_y, SET_VALUE = scaling_y
		WIDGET_CONTROL, scal_z, SET_VALUE = scaling_z
		break
	end
	'CUT_CO': begin
		WIDGET_CONTROL, event.id, GET_VALUE = cut_pos
		if (cut_pos le 0) then cut_pos = 0
		if (cut_pos gt max_pos) then cut_pos = max_pos
		if (not any (slice eq [1, 2, 3])) then cut_pos = -1
		if (event.update) then WIDGET_CONTROL, cut_co, SET_VALUE = cut_pos
		WIDGET_CONTROL, cut_sl, SET_VALUE = cut_pos
		break
	end
	'CUT_SL': begin
		WIDGET_CONTROL, event.id, GET_VALUE = cut_pos
		WIDGET_CONTROL, cut_co, SET_VALUE = cut_pos
		break
	end
	'SUB_XS': begin
		subvol_nx = subvol_xe - subvol_xs + 1
		WIDGET_CONTROL, event.id, GET_VALUE = subvol_xs
		subvol_xs = (subvol_xs > 0) < (nx-1)
		subvol_xe = (subvol_xs + subvol_nx - 1) < (nx-1)
		subvol_nx = subvol_xe - subvol_xs + 1
		WIDGET_CONTROL, event.id, SET_VALUE = subvol_xs
		WIDGET_CONTROL, sub_xe, SET_VALUE = subvol_xe
		WIDGET_CONTROL, sub_nx, SET_VALUE = subvol_nx
		pc_select_files_update
		break
	end
	'SUB_YS': begin
		subvol_ny = subvol_ye - subvol_ys + 1
		WIDGET_CONTROL, event.id, GET_VALUE = subvol_ys
		subvol_ys = (subvol_ys > 0) < (ny-1)
		subvol_ye = (subvol_ys + subvol_ny - 1) < (ny-1)
		subvol_ny = subvol_ye - subvol_ys + 1
		WIDGET_CONTROL, event.id, SET_VALUE = subvol_ys
		WIDGET_CONTROL, sub_ye, SET_VALUE = subvol_ye
		WIDGET_CONTROL, sub_ny, SET_VALUE = subvol_ny
		pc_select_files_update
		break
	end
	'SUB_ZS': begin
		subvol_nz = subvol_ze - subvol_zs + 1
		WIDGET_CONTROL, event.id, GET_VALUE = subvol_zs
		subvol_zs = (subvol_zs > 0) < (nz-1)
		subvol_ze = (subvol_zs + subvol_nz - 1) < (nz-1)
		subvol_nz = subvol_ze - subvol_zs + 1
		WIDGET_CONTROL, event.id, SET_VALUE = subvol_zs
		WIDGET_CONTROL, sub_ze, SET_VALUE = subvol_ze
		WIDGET_CONTROL, sub_nz, SET_VALUE = subvol_nz
		pc_select_files_update
		break
	end
	'SUB_XE': begin
		subvol_nx = subvol_xe - subvol_xs + 1
		last_xe = subvol_xe
		WIDGET_CONTROL, event.id, GET_VALUE = subvol_xe
		if (subvol_xe gt last_xe) then subvol_nx = subvol_xe - subvol_xs + 1
		subvol_xe = (subvol_xe > 0) < (nx-1)
		subvol_nx = (subvol_xe - subvol_xs + 1) > 1
		subvol_xs = subvol_xe - subvol_nx + 1
		WIDGET_CONTROL, event.id, SET_VALUE = subvol_xe
		WIDGET_CONTROL, sub_xs, SET_VALUE = subvol_xs
		WIDGET_CONTROL, sub_nx, SET_VALUE = subvol_nx
		pc_select_files_update
		break
	end
	'SUB_YE': begin
		subvol_ny = subvol_ye - subvol_ys + 1
		last_ye = subvol_ye
		WIDGET_CONTROL, event.id, GET_VALUE = subvol_ye
		if (subvol_ye gt last_ye) then subvol_ny = subvol_ye - subvol_ys + 1
		subvol_ye = (subvol_ye > 0) < (ny-1)
		subvol_ny = (subvol_ye - subvol_ys + 1) > 1
		subvol_ys = subvol_ye - subvol_ny + 1
		WIDGET_CONTROL, event.id, SET_VALUE = subvol_ye
		WIDGET_CONTROL, sub_ys, SET_VALUE = subvol_ys
		WIDGET_CONTROL, sub_ny, SET_VALUE = subvol_ny
		pc_select_files_update
		break
	end
	'SUB_ZE': begin
		subvol_nz = subvol_ze - subvol_zs + 1
		last_ze = subvol_ze
		WIDGET_CONTROL, event.id, GET_VALUE = subvol_ze
		if (subvol_ze gt last_ze) then subvol_nz = subvol_ze - subvol_zs + 1
		subvol_ze = (subvol_ze > 0) < (nz-1)
		subvol_nz = (subvol_ze - subvol_zs + 1) > 1
		subvol_zs = subvol_ze - subvol_nz + 1
		WIDGET_CONTROL, event.id, SET_VALUE = subvol_ze
		WIDGET_CONTROL, sub_zs, SET_VALUE = subvol_zs
		WIDGET_CONTROL, sub_nz, SET_VALUE = subvol_nz
		pc_select_files_update
		break
	end
	'SUB_NX': begin
		WIDGET_CONTROL, event.id, GET_VALUE = subvol_nx
		subvol_nx = (subvol_nx > 1) < nx
		subvol_xe = (subvol_xs + subvol_nx - 1) < (nx-1)
		subvol_xs = subvol_xe - subvol_nx + 1
		WIDGET_CONTROL, event.id, SET_VALUE = subvol_nx
		WIDGET_CONTROL, sub_xs, SET_VALUE = subvol_xs
		pc_select_files_update
		WIDGET_CONTROL, sub_xe, SET_VALUE = subvol_xe
		break
	end
	'SUB_NY': begin
		WIDGET_CONTROL, event.id, GET_VALUE = subvol_ny
		subvol_ny = (subvol_ny > 1) < ny
		subvol_ye = (subvol_ys + subvol_ny - 1) < (ny-1)
		subvol_ys = subvol_ye - subvol_ny + 1
		WIDGET_CONTROL, event.id, SET_VALUE = subvol_ny
		WIDGET_CONTROL, sub_ys, SET_VALUE = subvol_ys
		WIDGET_CONTROL, sub_ye, SET_VALUE = subvol_ye
		pc_select_files_update
		break
	end
	'SUB_NZ': begin
		WIDGET_CONTROL, event.id, GET_VALUE = subvol_nz
		subvol_nz = (subvol_nz > 1) < nz
		subvol_ze = (subvol_zs + subvol_nz - 1) < (nz-1)
		subvol_zs = subvol_ze - subvol_nz + 1
		WIDGET_CONTROL, event.id, SET_VALUE = subvol_nz
		WIDGET_CONTROL, sub_zs, SET_VALUE = subvol_zs
		WIDGET_CONTROL, sub_ze, SET_VALUE = subvol_ze
		pc_select_files_update
		break
	end
	'QUANT': begin
		quant_selected = WIDGET_INFO (c_quant, /LIST_SELECT)
		pc_select_files_update
		break
	end
	'OVER': begin
		over_selected = WIDGET_INFO (c_over, /LIST_SELECT)
		pc_select_files_update
		break
	end
	'O_DEF':
	'Q_DEF':
	'CONT': begin
		cont_selected = WIDGET_INFO (c_cont, /LIST_SELECT)
		if (eventval ne 'O_DEF') then pc_select_files_update_list, quant_list, all_quant, quant_selected, default=quant, avail=quant_avail
		if (eventval ne 'Q_DEF') then pc_select_files_update_list, over_list, all_over, over_selected, default=over, avail=over_avail
		pc_select_files_update, quant=(eventval ne 'O_DEF'), over=(eventval ne 'Q_DEF')
		break
	end
	'Q_ALL': begin
		quant_selected = where (strmid (quant_avail, 0, 5) ne "[N/A]")
		pc_select_files_update, /quant
		break
	end
	'Q_NONE': begin
		quant_selected = -1
		pc_select_files_update, /quant
		break
	end
	'O_ALL': begin
		over_selected = where (strmid (over_avail, 0, 5) ne "[N/A]")
		pc_select_files_update, /over
		break
	end
	'O_NONE': begin
		over_selected = -1
		pc_select_files_update, /over
		break
	end
	'SHOW_TIME': begin
		WIDGET_CONTROL, b_ts, SENSITIVE = 0
		pc_show_ts, obj=ts, unit=units, start_param=start_par, run_param=run_par, datadir=data_dir
		WIDGET_CONTROL, b_ts, SENSITIVE = 1
		break
	end
	'OK': begin
		selected = WIDGET_INFO (c_list, /LIST_SELECT)
		cont_selected = WIDGET_INFO (c_cont, /LIST_SELECT)
		if (not finite (c_quant, /NaN)) then quant_selected = WIDGET_INFO (c_quant, /LIST_SELECT)
		if (not finite (c_over, /NaN)) then over_selected = WIDGET_INFO (c_over, /LIST_SELECT)
		WIDGET_CONTROL, scal_x, GET_VALUE = scaling_x
		WIDGET_CONTROL, scal_y, GET_VALUE = scaling_y
		WIDGET_CONTROL, scal_z, GET_VALUE = scaling_z
		quit = event.top
		break
	end
	'CANCEL': begin
		selected = -1
		add_selected = 0
		var_selected = 0
		quit = event.top
		break
	end
	endswitch

	WIDGET_CONTROL, WIDGET_INFO (event.top, /CHILD)

	if (quit ge 0) then WIDGET_CONTROL, quit, /DESTROY

	return
end


; Update file size and memory consumption display.
pro pc_select_files_update, quant_update=quant_update, over_update=over_update

	common select_files_gui_common, b_var, b_add, b_ts, c_list, i_skip, i_step, f_load, f_comp, scal_x, scal_y, scal_z, c_cont, c_quant, c_over, d_slice, cut_co, cut_sl, sub_xs, sub_ys, sub_zs, sub_xe, sub_ye, sub_ze, sub_nx, sub_ny, sub_nz, i_last, nf
	common select_files_common, num_files, selected, num_selected, var_selected, add_selected, sources, sources_selected, num_cont, cont_selected, quant, quant_selected, quant_list, all_quant, quant_avail, over, over_selected, over_list, all_over, over_avail, cut_pos, max_pos, slice, skipping, stepping, data_dir, units, run_par, start_par, gb_per_file, cont_corr, subvol_corr, subvol_xs, subvol_ys, subvol_zs, subvol_xe, subvol_ye, subvol_ze, scaling_x, scaling_y, scaling_z, nx, ny, nz, nghost, skip_last, n_slice

	num = 0
	for pos = 0, num_cont-1 do begin
		if (not any (pos eq cont_selected)) then continue
		if (any (strcmp (sources[pos], ["uu", "aa"], /fold_case))) then num += 3 else num++
	end
	cont_corr = num / float (num_cont)

	if ((slice ge 1) and (slice le 3)) then begin
		subvol_corr = (1.0+2*nghost)/(max_pos+1.0+2*nghost)
	end else if (slice eq 4) then begin
		subvol_nx = subvol_xe - subvol_xs + 1
		subvol_ny = subvol_ye - subvol_ys + 1
		subvol_nz = subvol_ze - subvol_zs + 1
		subvol_corr = product (([subvol_nx, subvol_ny, subvol_nz]+2*nghost)/([nx, ny, nz]+2.0*nghost))
	end else begin
		subvol_corr = 1
	end
	
	if (any (quant_selected ne -1)) then num_quant = n_elements (quant_selected) else num_quant = 0
	if (any (over_selected ne -1)) then num_over = n_elements (over_selected) else num_over = 0

	WIDGET_CONTROL, f_load, SET_VALUE = gb_per_file*cont_corr*subvol_corr*(num_selected+var_selected+add_selected)
	if (not finite (f_comp, /NaN)) then $
		WIDGET_CONTROL, f_comp, SET_VALUE = gb_per_file/num_cont*subvol_corr*(num_selected+var_selected+add_selected)*(num_quant+num_over*3)

	if (not finite (c_quant, /NaN) and keyword_set (quant_update)) then begin
		WIDGET_CONTROL, c_quant, SET_VALUE = quant_avail
		WIDGET_CONTROL, c_quant, SET_LIST_SELECT = quant_selected
	end
	if (not finite (c_over, /NaN) and keyword_set (over_update)) then begin
		WIDGET_CONTROL, c_over, SET_VALUE = over_avail
		WIDGET_CONTROL, c_over, SET_LIST_SELECT = over_selected
	end
end


; File selection dialog GUI.
pro pc_select_files, files=files, num_selected=num, pattern=pattern, varfile=varfile, addfile=addfile, datadir=datadir, allprocs=allprocs, reduced=reduced, procdir=procdir, unit=unit, dim=dim, start_param=start_param, run_param=run_param, quantities=quantities, overplots=overplots, varcontent=varcontent, var_list=var_list, cut_x=cut_x, cut_y=cut_y, cut_z=cut_z, xs=xs, ys=ys, zs=zs, xe=xe, ye=ye, ze=ze, min_display=min_display, max_display=max_display, hide_quantities=hide_quantities, hide_overplots=hide_overplots, scaling=scaling

	common select_files_gui_common, b_var, b_add, b_ts, c_list, i_skip, i_step, f_load, f_comp, scal_x, scal_y, scal_z, c_cont, c_quant, c_over, d_slice, cut_co, cut_sl, sub_xs, sub_ys, sub_zs, sub_xe, sub_ye, sub_ze, sub_nx, sub_ny, sub_nz, i_last, nf
	common select_files_common, num_files, selected, num_selected, var_selected, add_selected, sources, sources_selected, num_cont, cont_selected, quant, quant_selected, quant_list, all_quant, quant_avail, over, over_selected, over_list, all_over, over_avail, cut_pos, max_pos, slice, skipping, stepping, data_dir, units, run_par, start_par, gb_per_file, cont_corr, subvol_corr, subvol_xs, subvol_ys, subvol_zs, subvol_xe, subvol_ye, subvol_ze, scaling_x, scaling_y, scaling_z, nx, ny, nz, nghost, skip_last, n_slice

	; Default settings
	@pc_gui_settings
	if (size (datadir, /type) eq 0) then datadir = pc_get_datadir (datadir)
	default_varfile = 'var.dat'
	default_crashfile = 'crash.dat'
	if (file_test (datadir+'/allprocs/var.h5')) then begin
		default_varfile = 'var.h5'
		default_crashfile = 'crash.h5'
	end
	default, pattern, "VAR[0-9]*"
	default, varfile, default_varfile
	default, addfile, default_crashfile
	default, skipping, 0
	default, stepping, 10
	default, skip_last, 0
	default, cut_x, -1
	default, cut_y, -1
	default, cut_z, -1
	default, min_display, 25
	default, max_display, 25
	if (max_display lt 1) then max_display = 100
	if (min_display gt max_display) then min_display = max_display
	if (min_display lt 1) then min_display = 100 < max_display
	if (keyword_set (hide_quantities)) then hide_quant = hide_quantities else hide_quant = 0
	if (keyword_set (hide_overplots)) then hide_over = hide_overplots else hide_over = 0
	if (not keyword_set (scaling)) then scaling = 1

	; Load needed objects
	pc_find_config, varfile, datadir=datadir, procdir=procdir, dim=dim, allprocs=allprocs, reduced=reduced, f77=f77, start_param=start_param
	if (not keyword_set (varcontent)) then varcontent = pc_varcontent (datadir=datadir, dim=dim, param=start_param, /quiet)
	sources = varcontent
	all_quant = pc_check_quantities (sources=sources, /available)
	quant = all_quant
	if (keyword_set (quantities)) then begin
		quant = quantities
		if (size (quantities, /type) eq 7) then begin
			num = n_elements (quantities)
			for pos = 0, num-1 do begin
				if (pos eq 0) then begin
					quant = create_struct (quantities[pos], quantities[pos])
				end else begin
					quant = create_struct (quant, quantities[pos], quantities[pos])
				end
			end
		end
	end
	all_over = pc_check_quantities (sources=sources, /vectorfields)
	if (keyword_set (overplots)) then over = overplots else over = all_over

	; Fill common blocks
	if (datadir eq "") then datadir = "."
	data_dir = datadir
	if (n_elements (scaling) lt 3) then scaling = [ scaling, scaling, scaling ]
	scaling_x = scaling[0]
	scaling_y = scaling[1]
	scaling_z = scaling[2]
	nx = dim.nx
	ny = dim.ny
	nz = dim.nz
	nghost = max ([dim.nghostx, dim.nghosty, dim.nghostz])
	if (not keyword_set (unit)) then pc_units, obj=unit, datadir=datadir, param=start_param, dim=dim, /quiet
	if (not has_tag (unit, "default_length")) then unit = create_struct (unit, display_units)
	units = unit
	start_par = start_param
	if (keyword_set (run_param)) then run_par = run_param
	f_comp = !VALUES.D_NAN
	c_quant = !VALUES.D_NAN
	c_over = !VALUES.D_NAN

	; Get list of available snapshots
	files = file_search (procdir, pattern)
	num_files = n_elements (files)
	stepping = stepping < (num_files - skipping - skip_last)
	for pos = 0, num_files - 1 do begin
		cut = strpos (procdir, "/", /REVERSE_SEARCH)
		if (strmid (files[pos], cut, 1) eq "/") then cut += 1
		files[pos] = strmid (files[pos], cut)
	end
	sorted = files
	sorted[*] = ""
	target = 0
	for len = min (strlen (files)), max (strlen (files)) do begin
		indices = where (strlen (files) eq len, num_match)
		if (num_match ge 1) then begin
			sub = files[indices]
			reorder = sort (sub)
			for pos = 0, num_match-1 do begin
				sorted[target] = sub[reorder[pos]]
				target += 1
			end
		end
	end
	files = sorted
	sorted = ""

	; Get file size
	if (file_test (procdir+varfile)) then somefile = varfile else somefile = files[0]
	file_struct = file_info (procdir+somefile)
	multiplier = dim.nprocx * dim.nprocy * dim.nprocz
	if (allprocs eq 1) then multiplier = 1
	if (allprocs ge 2) then multiplier = dim.nprocz
	gb_per_file = (file_struct.size * multiplier) / 1024.0^3

	; Preselected snapshots and additional snapshots to be loaded
	selected = -1
	num_selected = 0
	if ((addfile ne "") and file_test (procdir+addfile)) then add_selected = 1 else add_selected = 0
	if ((varfile ne "") and file_test (procdir+varfile)) then var_selected = 1 else var_selected = 0

	; Pre-defined sub-volume settings
	if (n_elements (xs) gt 0) then subvol_xs = xs
	if (n_elements (ys) gt 0) then subvol_ys = ys
	if (n_elements (zs) gt 0) then subvol_zs = zs
	if (n_elements (xe) gt 0) then subvol_xe = xe
	if (n_elements (ye) gt 0) then subvol_ye = ye
	if (n_elements (ze) gt 0) then subvol_ze = ze
	default, subvol_xs, 0
	default, subvol_ys, 0
	default, subvol_zs, 0
	default, subvol_xe, nx-1
	default, subvol_ye, ny-1
	default, subvol_ze, nz-1
	subvol_xs = (subvol_xs < (nx-1)) > 0
	subvol_ys = (subvol_ys < (ny-1)) > 0
	subvol_zs = (subvol_zs < (nz-1)) > 0
	subvol_xe = subvol_xe < (nx-1) > subvol_xs
	subvol_ye = subvol_ye < (ny-1) > subvol_ys
	subvol_ze = subvol_ze < (nz-1) > subvol_zs

	; Pre-defined slice settings
	dimensionality = 0 + ((nx gt 1) + (ny gt 1) + (nz gt 1))
	slice = 0
	if ((subvol_xs gt 0) or (subvol_xe lt nx-1) or (subvol_ys gt 0) or (subvol_ye lt ny-1) or (subvol_zs gt 0) or (subvol_ze lt nz-1)) then slice = 4
	max_pos = -1
	cut_pos = -1
	if (dimensionality eq 3) then begin
		if (cut_x ge 0) then begin
			slice = 1
			max_pos = nx
			cut_pos = cut_x < max_pos
			subvol_xs = cut_pos
			subvol_xe = cut_pos
			subvol_ys = 0
			subvol_ye = ny-1
			subvol_zs = 0
			subvol_ze = nz-1
		end
		if (cut_y ge 0) then begin
			slice = 2
			max_pos = ny
			cut_pos = cut_y < max_pos
			subvol_xs = 0
			subvol_xe = nx-1
			subvol_ys = cut_pos
			subvol_ye = cut_pos
			subvol_zs = 0
			subvol_ze = nz-1
		end
		if (cut_z ge 0) then begin
			slice = 3
			max_pos = nz
			cut_pos = cut_z < max_pos
			subvol_xs = 0
			subvol_xe = nx-1
			subvol_ys = 0
			subvol_ye = ny-1
			subvol_zs = cut_pos
			subvol_ze = cut_pos
		end
	end

	subvol_nx = subvol_xe - subvol_xs + 1
	subvol_ny = subvol_ye - subvol_ys + 1
	subvol_nz = subvol_ze - subvol_zs + 1

	; Pre-defined varcontent settings
	content = varcontent.variable
	num_cont = n_elements (content)
	indices = where (content ne "UNKNOWN")
	if (any (indices ne -1)) then content = content[indices]
	num_content = n_elements (content)
	cont_selected = indgen (num_content)
	if (n_elements (var_list) gt 0) then begin
		for pos = 0, num_content-1 do begin
			if (not any (strcmp (sources[pos], var_list, /fold_case))) then cont_selected[pos] = -1
		end
		indices = where (cont_selected ge 0)
		if (any (indices ne -1)) then cont_selected = cont_selected[indices]
	end

	; Pre-defined quantities settings
	num_all_quant = n_tags (all_quant)
	quant_list = strarr (num_all_quant)
	num_quant = n_tags (quant)
	for pos = 0, num_all_quant-1 do quant_list[pos] = all_quant.(pos)
	pc_select_files_update_list, quant_list, all_quant, quant_selected, default=quant, avail=quant_avail

	; Pre-defined overplot settings
	num_all_over = n_tags (all_over)
	over_list = strarr (num_all_over)
	num_over = n_tags (over)
	for pos = 0, num_all_over-1 do over_list[pos] = all_over.(pos)
	pc_select_files_update_list, over_list, all_over, over_selected, default=over, avail=over_avail


	; Build GUI
	MOTHER	= WIDGET_BASE (title='PC file selector', EVENT_PRO=select_files_event)
	BASE	= WIDGET_BASE (MOTHER, /row)

	CTRL	= WIDGET_BASE (BASE, /col)

	tmp	= WIDGET_LABEL (CTRL, value='File Selection:', frame=0)
	EDIT	= WIDGET_BASE (CTRL, /row, /align_center)
	tmp	= WIDGET_BUTTON (EDIT, xsize=60, value='ALL', uvalue='ALL')
	tmp	= WIDGET_BUTTON (EDIT, xsize=60, value='NONE', uvalue='NONE')

	EDIT	= WIDGET_BASE (CTRL, /col, xsize=100, /align_center, frame=1)
	i_skip	= CW_FIELD (EDIT, title='skip first', /column, uvalue='SKIP', value=skipping, /integer, /return_events, xsize=8)
	i_step	= CW_FIELD (EDIT, title='stepping', /column, uvalue='STEP', value=stepping, /integer, /return_events, xsize=8)
	i_last  = CW_FIELD (EDIT, title='skip last', /column, uvalue='SKIP_LAST', value=skip_last, /integer, /return_events, xsize=8)
	tmp	= WIDGET_BUTTON (EDIT, xsize=80, value='APPLY', uvalue='APPLY')

	XTRA	= WIDGET_BASE (CTRL, /col, /align_center)
	tmp	= WIDGET_LABEL (XTRA, value='Analysis:', frame=0)
	BUT	= WIDGET_BASE (XTRA, /row, /align_center)
	b_ts	= WIDGET_BUTTON (BUT, value='show timeseries', uvalue='SHOW_TIME')

	tmp	= WIDGET_LABEL (XTRA, value='GB per snapshot:', frame=0)
	tmp	= WIDGET_LABEL (XTRA, value=strtrim (gb_per_file, 2), frame=1, xsize=100)
	f_load	= CW_FIELD (XTRA, title='Total GB to load:', /column, /float, /noedit)
	if (not hide_quant or not hide_over) then $
		f_comp	= CW_FIELD (XTRA, title='Total GB to compute:', /column, /float, /noedit)

	BUT	= WIDGET_BASE (XTRA, /row, /align_center, frame=1)
	tmp	= WIDGET_BUTTON (BUT, xsize=60, value='CANCEL', uvalue='CANCEL')
	tmp	= WIDGET_BUTTON (BUT, xsize=60, value='OK', uvalue='OK')

	SEL	= WIDGET_BASE (BASE, /col)

	tmp	= WIDGET_LABEL (SEL, value='Available snapshots:', frame=0)
	c_list	= WIDGET_LIST (SEL, value=files, uvalue='LIST', YSIZE=(num_files<max_display)>min_display, /multiple)
	tmp	= WIDGET_LABEL (SEL, value='Selceted snapshots:', frame=0)
	nf	= CW_FIELD (SEL, title='', value=strtrim(string(num_selected),2) + ' / ' + strtrim(string(num_files),2), xsize=13, /noedit)
	if (var_selected eq 1) then begin
		pc_read_var_time, t=var_time, varfile=varfile, datadir=datadir, allprocs=allprocs, reduced=reduced, /quiet
		b_var	= CW_BGROUP (SEL, varfile+' ('+strtrim (var_time*unit.time, 2)+' s)', /nonexcl, uvalue='VAR', set_value=1)
	end else if (varfile ne "") then begin
		b_var	= WIDGET_LABEL (SEL, value='"'+varfile+'" not found', frame=0)
	end else begin
		b_var	= WIDGET_LABEL (SEL, value='No varfile selected', frame=0)
	end

	if (addfile eq varfile) then begin
		b_add	= WIDGET_LABEL (SEL, value='Additional file is identical to "'+varfile+'"', frame=0)
	end else if (add_selected eq 1) then begin
		pc_read_var_time, t=add_time, varfile=addfile, datadir=datadir, allprocs=allprocs, reduced=reduced, /quiet
		b_add	= CW_BGROUP (SEL, addfile+' ('+strtrim (add_time*unit.time, 2)+' s)', /nonexcl, uvalue='ADD', set_value=1)
	end else begin
		b_add	= WIDGET_LABEL (SEL, value='No "'+addfile+'" found', frame=0)
	end

	VC	= WIDGET_BASE (BASE, /col)

	tmp	= WIDGET_LABEL (VC, value='Available content:', frame=0)
	c_cont	= WIDGET_LIST (VC, value=content, uvalue='CONT', YSIZE=num_content<max_display, /multiple)
	WIDGET_CONTROL, c_cont, SET_LIST_SELECT = cont_selected

	IO_scheme = ["distributed files", "collective files", "collect_xy files"]
	if (keyword_set (reduced)) then IO_scheme[allprocs] = "reduced files"
	tmp	= WIDGET_LABEL (VC, value='Load '+IO_scheme[allprocs]+":", frame=0)
	SEL	= WIDGET_BASE (VC, frame=1, /align_center, /col)

	case dimensionality of 
		1: begin
			if (allprocs eq 1) then dir_str = "/allprocs/" else dir_str = "/proc*/"
			if (keyword_set (reduced)) then dir_str = "/reduced/"
			tmp	= WIDGET_LABEL (SEL, value="From: "+datadir+dir_str, frame=0)
			tmp	= WIDGET_LABEL (SEL, value='full '+strtrim (dimensionality, 2)+'D data', frame=0)
			if nx gt 1 then begin
				load_list = ['full 1D data', '1D sub-volume (x)']
			endif	
			if ny gt 1 then begin
				load_list = ['full 1D data', '1D sub-volume (y)']
			endif
			if nz gt 1 then begin
				load_list = ['full 1D data', '1D sub-volume (z)']
			endif
			n_slice = n_elements(load_list)
			d_slice	= WIDGET_DROPLIST (SEL, value=load_list, /align_center, uvalue='SLICE')
			WIDGET_CONTROL, d_slice, SET_DROPLIST_SELECT = slice
			cut_co	= CW_FIELD (SEL, title='Slice position:', uvalue='CUT_CO', value=cut_pos>0, /integer, /return_events, xsize=8)
			WIDGET_CONTROL, cut_co, SENSITIVE = (cut_pos ge 0)
			cut_sl	= WIDGET_SLIDER (SEL, uvalue='CUT_SL', value=cut_pos, min=0<max_pos, max=max_pos>1, /drag, /suppress_value, sensitive=(cut_pos ge 0))
		   end
		2: begin
			if (allprocs eq 1) then dir_str = "/allprocs/" else dir_str = "/proc*/"
			if (keyword_set (reduced)) then dir_str = "/reduced/"
			if (nx GT 1) and (ny GT 1) then begin
				load_list = ['full 2D data', 'y-slice', 'x-slice', '2D sub-volume']
			endif 
			if (nx GT 1) and (nz GT 1) then begin
				load_list = ['full 2D data', 'z-slice', 'x-slice', '2D sub-volume']
			endif 
			if (ny GT 1) and (nz GT 1) then begin
				load_list = ['full 2D data', 'z-slice', 'y-slice', '2D sub-volume']
			endif 
			n_slice = n_elements(load_list)
			tmp	= WIDGET_LABEL (SEL, value="From: "+datadir+dir_str, frame=0)
			tmp	= WIDGET_LABEL (SEL, value='full '+strtrim (dimensionality, 2)+'D data', frame=0)	
			d_slice	= WIDGET_DROPLIST (SEL, value=load_list, /align_center, uvalue='SLICE')
			WIDGET_CONTROL, d_slice, SET_DROPLIST_SELECT = slice
			cut_co	= CW_FIELD (SEL, title='Slice position:', uvalue='CUT_CO', value=cut_pos>0, /integer, /return_events, xsize=8)
			WIDGET_CONTROL, cut_co, SENSITIVE = (cut_pos ge 0)
			cut_sl	= WIDGET_SLIDER (SEL, uvalue='CUT_SL', value=cut_pos, min=0<max_pos, max=max_pos>1, /drag, /suppress_value, sensitive=(cut_pos ge 0))

		   end
		3: begin
			load_list = ['full 3D data', 'yz-slice', 'xz-slice', 'xy-slice', '3D sub-volume']
			n_slice = n_elements(load_list)
			tmp	= WIDGET_LABEL (SEL, value="From: "+procdir, frame=0)
			d_slice	= WIDGET_DROPLIST (SEL, value=load_list, /align_center, uvalue='SLICE')
			WIDGET_CONTROL, d_slice, SET_DROPLIST_SELECT = slice
			cut_co	= CW_FIELD (SEL, title='Slice position:', uvalue='CUT_CO', value=cut_pos>0, /integer, /return_events, xsize=8)
			WIDGET_CONTROL, cut_co, SENSITIVE = (cut_pos ge 0)
			cut_sl	= WIDGET_SLIDER (SEL, uvalue='CUT_SL', value=cut_pos, min=0<max_pos, max=max_pos>1, /drag, /suppress_value, sensitive=(cut_pos ge 0))
		   end
	endcase
	tmp	= WIDGET_LABEL (SEL, value='Sub-volume (start, end, size):', frame=0)
	SUB	= WIDGET_BASE (SEL, frame=0, /align_center, /row)
	sub_xs	= CW_FIELD (SUB, title='X:', uvalue='SUB_XS', value=subvol_xs, /integer, /return_events, xsize=5)
	sub_xe	= CW_FIELD (SUB, title='', uvalue='SUB_XE', value=subvol_xe, /integer, /return_events, xsize=5)
	sub_nx	= CW_FIELD (SUB, title='=', uvalue='SUB_NX', value=subvol_nx, /integer, /return_events, xsize=5)
	SUB	= WIDGET_BASE (SEL, frame=0, /align_center, /row)
	sub_ys	= CW_FIELD (SUB, title='Y:', uvalue='SUB_YS', value=subvol_ys, /integer, /return_events, xsize=5)
	sub_ye	= CW_FIELD (SUB, title='', uvalue='SUB_YE', value=subvol_ye, /integer, /return_events, xsize=5)
	sub_ny	= CW_FIELD (SUB, title='=', uvalue='SUB_NY', value=subvol_ny, /integer, /return_events, xsize=5)
	SUB	= WIDGET_BASE (SEL, frame=0, /align_center, /row)
	sub_zs	= CW_FIELD (SUB, title='Z:', uvalue='SUB_ZS', value=subvol_zs, /integer, /return_events, xsize=5)
	sub_ze	= CW_FIELD (SUB, title='', uvalue='SUB_ZE', value=subvol_ze, /integer, /return_events, xsize=5)
	sub_nz	= CW_FIELD (SUB, title='=', uvalue='SUB_NZ', value=subvol_nz, /integer, /return_events, xsize=5)
	WIDGET_CONTROL, sub_xs, SENSITIVE = ((slice eq n_slice-1) and (nx gt 1)) 
	WIDGET_CONTROL, sub_xe, SENSITIVE = ((slice eq n_slice-1) and (nx gt 1))
	WIDGET_CONTROL, sub_nx, SENSITIVE = ((slice eq n_slice-1) and (nx gt 1))
	WIDGET_CONTROL, sub_ys, SENSITIVE = ((slice eq n_slice-1) and (ny gt 1))
	WIDGET_CONTROL, sub_ye, SENSITIVE = ((slice eq n_slice-1) and (ny gt 1))
	WIDGET_CONTROL, sub_ny, SENSITIVE = ((slice eq n_slice-1) and (ny gt 1))
	WIDGET_CONTROL, sub_zs, SENSITIVE = ((slice eq n_slice-1) and (nz gt 1))
	WIDGET_CONTROL, sub_ze, SENSITIVE = ((slice eq n_slice-1) and (nz gt 1))
	WIDGET_CONTROL, sub_nz, SENSITIVE = ((slice eq n_slice-1) and (nz gt 1))

	tmp	= WIDGET_LABEL (VC, value='Scaling factors (X,Y,Z):', frame=0)
	SCA	= WIDGET_BASE (VC, frame=1, /align_center, /row)
	scal_x	= CW_FIELD (SCA, title='', uvalue='SCAL_X', value=scaling_x, /float, /return_events, xsize=5)
	scal_y	= CW_FIELD (SCA, title='', uvalue='SCAL_Y', value=scaling_y, /float, /return_events, xsize=5)
	scal_z	= CW_FIELD (SCA, title='', uvalue='SCAL_Z', value=scaling_z, /float, /return_events, xsize=5)
	tmp	= WIDGET_BUTTON (SCA, xsize=16, value='+', uvalue='SCAL_PLUS')
	tmp	= WIDGET_BUTTON (SCA, xsize=16, value='-', uvalue='SCAL_MINUS')

	if (not hide_quant) then begin
		QU	= WIDGET_BASE (BASE, /col)
		tmp	= WIDGET_LABEL (QU, value='Available quantities:', frame=0)
		c_quant	= WIDGET_LIST (QU, value=quant_avail, uvalue='QUANT', YSIZE=(num_quant<max_display)>min_display, /multiple)
		SEL	= WIDGET_BASE (QU, /row, /align_center)
		tmp	= WIDGET_BUTTON (SEL, xsize=40, value='ALL', uvalue='Q_ALL')
		tmp	= WIDGET_BUTTON (SEL, xsize=60, value='DEFAULT', uvalue='Q_DEF')
		tmp	= WIDGET_BUTTON (SEL, xsize=40, value='NONE', uvalue='Q_NONE')
	end

	if (not hide_over) then begin
		OV	= WIDGET_BASE (BASE, /col)
		tmp	= WIDGET_LABEL (OV, value='Available overplots:', frame=0)
		c_over	= WIDGET_LIST (OV, value=over_avail, uvalue='OVER', YSIZE=(num_over<max_display)>min_display, /multiple)
		SEL	= WIDGET_BASE (OV, /row, /align_center)
		tmp	= WIDGET_BUTTON (SEL, xsize=40, value='ALL', uvalue='O_ALL')
		tmp	= WIDGET_BUTTON (SEL, xsize=60, value='DEFAULT', uvalue='O_DEF')
		tmp	= WIDGET_BUTTON (SEL, xsize=40, value='NONE', uvalue='O_NONE')
	end

	pc_select_files_update, /quant, /over

	WIDGET_CONTROL, MOTHER, /REALIZE
	wimg = !d.window

	WIDGET_CONTROL, BASE

	XMANAGER, "select_files", MOTHER


	; Build list of selected snapshots
	if (any (selected ne -1)) then files = files[selected] else files = [ "" ]
	if (var_selected) then files = [ files, varfile ]
	if (add_selected) then files = [ files, addfile ]
	indices = where (files ne "")
	if (any (indices ne -1)) then begin
		files = files[indices]
		num_selected = n_elements (files)
	end else begin
		files = -1
		num_selected = 0
	end

	; Build list of selected variables
	if (any (cont_selected ge 0)) then begin
		content = varcontent.idlvar
		content = content[where (content ne "dummy")]
		var_list = content[cont_selected]
	end else begin
		var_list = -1
	end

	; Build indices of selected sub-volume
	cut_x = -1
	cut_y = -1
	cut_z = -1
	xs = 0
	ys = 0
	zs = 0
	xe = nx-1
	ye = ny-1
	ze = nz-1
	; selects the right xs, xe, ys, ye, zs, ze values in 1D and 2D case
	case dimensionality of
		1: begin
			if (slice eq 1) then begin
				if (nx gt 1) then begin
					xs = subvol_xs
					xe = subvol_xe
				end
				if (ny gt 1) then begin
					ys = subvol_ys
					ye = subvol_ye
				end
				if (nz gt 1) then begin
					zs = subvol_zs
					ze = subvol_ze
				end
			end
		end
		
		2: begin
			if (nx eq 1) then begin
				if (slice eq 1) then begin
					cut_y = cut_pos
					ys = cut_y
					ye = cut_y
				end
				
				if (slice eq 2) then begin
					cut_z = cut_pos
					zs = cut_z
					ze = cut_z
				end
				
				if (slice eq 3) then begin
					ys = subvol_ys
					zs = subvol_zs
					ye = subvol_ye
					ze = subvol_ze
				end
			end
			if (ny eq 1) then begin
				if (slice eq 1) then begin
					cut_x = cut_pos
					xs = cut_x
					xe = cut_x
				end
				
				if (slice eq 2) then begin
					cut_z = cut_pos
					zs = cut_z
					ze = cut_z
				end
				if (slice eq 3) then begin
					xs = subvol_xs
					zs = subvol_zs
					xe = subvol_xe
					ze = subvol_ze
				end
			end
			if (nz eq 1) then begin
				if (slice eq 1) then begin
					cut_x = cut_pos
					xs = cut_x
					xe = cut_x
				end
				
				if (slice eq 2) then begin
					cut_y = cut_pos
					ys = cut_y
					ye = cut_y
				end
				if (slice eq 3) then begin
					xs = subvol_xs
					ys = subvol_ys
					xe = subvol_xe
					ye = subvol_ye
				end
			end
		end
		
		3: begin
			if (slice eq 1) then begin
				cut_x = cut_pos
				xs = cut_x
				xe = cut_x
			end
			if (slice eq 2) then begin
				cut_y = cut_pos
				ys = cut_y
				ye = cut_y
			end
			if (slice eq 3) then begin	
				cut_z = cut_pos
				zs = cut_z
				ze = cut_z
			end
			if (slice eq 4) then begin
				xs = subvol_xs
				ys = subvol_ys
				zs = subvol_zs
				xe = subvol_xe
				ye = subvol_ye
				ze = subvol_ze
			end
		end
	endcase
	; Build scaling factors array
	scaling = [ scaling_x, scaling_y, scaling_z ]

	; Build list of selected quantities
	if (any (quant_selected ge 0)) then begin
		num = n_elements (quant_selected)
		tags = tag_names (all_quant)
		for pos=0, num-1 do begin
			tag = tags[quant_selected[pos]]
			if (pos eq 0) then begin
				quantities = create_struct (tag, all_quant.(quant_selected[pos]))
			end else begin
				quantities = create_struct (quantities, tag, all_quant.(quant_selected[pos]))
			end
		end
	end else begin
		undefine, quantities
		print, "No physical quantities have been selected for computation."
	end

	; Build list of selected overplots
	if (any (over_selected ge 0)) then begin
		num = n_elements (over_selected)
		tags = tag_names (all_over)
		for pos=0, num-1 do begin
			tag = tags[over_selected[pos]]
			if (pos eq 0) then begin
				overplots = create_struct (tag, all_over.(over_selected[pos]))
			end else begin
				overplots = create_struct (overplots, tag, all_over.(over_selected[pos]))
			end
		end
	end else begin
		overplots = { none:'none' }
	end

	if (not keyword_set (quiet)) then begin
		; Print summary
		print, ""
		print, "Selected snapshots:"
		if (num_selected le 0) then begin
			print, "none"
		end else begin
			print, files
			print, "This corresponds to ", strtrim (num_selected * gb_per_file * cont_corr * subvol_corr, 2), " GB = ", strtrim (num_selected, 2), " files"
		end
		print, ""
	end

	num = num_selected
end
