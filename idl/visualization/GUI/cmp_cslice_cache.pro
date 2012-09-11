;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   cmp_cslice_cache.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  $Id$
;;;
;;;  Description:
;;;   Fast and powerful to use tool to view and compare slices of 3D data
;;;  To do:
;;;   Add more comments


; Check value range and extend it, if necessary (for sliders or plotting)
function cslice_get_range, data

	common slider_common, bin_x, bin_y, bin_z, num_x, num_y, num_z, pos_b, pos_t, pos_over, val_min, val_max, val_range, over_max, dimensionality, frozen
	common settings_common, px, py, pz, cut, log_plot, abs_scale, show_cross, show_cuts, sub_aver, selected_cube, selected_overplot, selected_snapshot, af_x, af_y, af_z, destretch

	if (abs_scale) then begin
		min = val_min
		max = val_max
	end else begin
		cslice_get_minmax_value, data, min, max
	end

	range = get_val_range ([min, max])

	return, [range[0], range[1]]
end


; Get values for minimum and maximum of the selected data
pro cslice_get_minmax_value, data, min, max

	tmp_minmax = minmax (data)

	min = min (tmp_minmax)
	max = max (tmp_minmax)
end


; Event handling of visualisation window
pro cslice_event, event

	common varset_common, set, overplot, oversets, unit, coord, varsets, varfiles, datadir, sources, param, run_param, var_list
	common event_common, button_pressed_yz, button_pressed_xz, button_pressed_xy
	common cslice_common, cube, field, num_cubes, num_overs, num_snapshots
	common slider_common, bin_x, bin_y, bin_z, num_x, num_y, num_z, pos_b, pos_t, pos_over, val_min, val_max, val_range, over_max, dimensionality, frozen
	common gui_common, wimg_yz, wimg_xz, wimg_xy, wcut_x, wcut_y, wcut_z, co_t, sl_t, co_x, co_y, co_z, sl_x, sl_y, sl_z, b_load, b_log, b_sub, b_abs, b_cro, b_des, aver, timeser, vars, over, snap, prev, next, play, image, sl_min, sl_max, sl_over, min_max, freeze, jump_min, jump_max
	common settings_common, px, py, pz, cut, log_plot, abs_scale, show_cross, show_cuts, sub_aver, selected_cube, selected_overplot, selected_snapshot, af_x, af_y, af_z, destretch

	WIDGET_CONTROL, WIDGET_INFO (event.top, /CHILD)

	; SETTINGS:
	; filename for saving the settings
	settings_file = 'GUI_settings.xdr'

	; time in seconds to wait after showing each frame for 1D, 2D, and 3D movies (0=fastest possible)
	min_wait_time = [ 0.2, 0.1, 0.05 ]

	quit = -1
	DRAW_IMAGE_1=0  &  DRAW_IMAGE_2=0  &  DRAW_IMAGE_3=0

	WIDGET_CONTROL, event.id, GET_UVALUE = eventval


	SWITCH eventval of
	'LOG_PLOT': begin
		log_plot = event.select
		DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
		break
	end
	'ABS_SCALE': begin
		abs_scale = event.select
		DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
		break
	end
	'SUB_AVER': begin
		pos_b[selected_cube,sub_aver] = val_min
		pos_t[selected_cube,sub_aver] = val_max
		sub_aver = event.select
		cslice_prepare_cube, -1
		DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
		break
	end
	'SHOW_CROSS': begin
		show_cross = event.select
		DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
		break
	end
	'DESTRETCH': begin
		destretch[where (coord.lequidist eq 0)] = event.select
		cslice_prepare_overplot
		DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
		break
	end
	'SHOW_AVER': begin
		cslice_draw_averages, selected_snapshot
		break
	end
	'SHOW_TIME': begin
		pc_show_ts, unit=unit, param=param, run_param=run_param, datadir=datadir
		break
	end
	'COT':
	'SLT': begin
		WIDGET_CONTROL, event.id, GET_VALUE = pos
		if (eventval eq 'COT') then begin
			diff = abs (pos - varfiles[*].time)
			pos = (where (diff eq min (diff)))[0]
			if (event.update) then WIDGET_CONTROL, co_t, SET_VALUE = varfiles[pos].time
			WIDGET_CONTROL, sl_t, SET_VALUE = num_snapshots-1-pos
		end else begin
			pos = num_snapshots - 1 - pos
			WIDGET_CONTROL, co_t, SET_VALUE = varfiles[pos].time
		end
		if (pos ne selected_snapshot) then begin
			cslice_prepare_set, pos
			cslice_prepare_cube, -1
			cslice_prepare_overplot
			WIDGET_CONTROL, snap, SET_DROPLIST_SELECT = selected_snapshot
			WIDGET_CONTROL, prev, SENSITIVE = (selected_snapshot lt num_snapshots - 1)
			WIDGET_CONTROL, next, SENSITIVE = (selected_snapshot gt 0)
			DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
		end
		break
	end
	'COX': begin
		WIDGET_CONTROL, event.id, GET_VALUE = pos
		px = pc_find_index (pos, coord.x, num_x)
		if (event.update) then WIDGET_CONTROL, co_x, SET_VALUE = coord.x[px]
		WIDGET_CONTROL, sl_x, SET_VALUE = px
		DRAW_IMAGE_1 = 1
		break
	end
	'COY': begin
		WIDGET_CONTROL, event.id, GET_VALUE = pos
		py = pc_find_index (pos, coord.y, num_y)
		if (event.update) then WIDGET_CONTROL, co_y, SET_VALUE = coord.y[py]
		WIDGET_CONTROL, sl_y, SET_VALUE = py
		DRAW_IMAGE_2 = 1
		break
	end
	'COZ': begin
		WIDGET_CONTROL, event.id, GET_VALUE = pos
		pz = pc_find_index (pos, coord.z, num_z)
		if (event.update) then WIDGET_CONTROL, co_z, SET_VALUE = coord.z[pz]
		WIDGET_CONTROL, sl_z, SET_VALUE = pz
		DRAW_IMAGE_3 = 1
		break
	end
	'SLX': begin
		WIDGET_CONTROL, event.id, GET_VALUE = pos
		px = pos
		WIDGET_CONTROL, co_x, SET_VALUE = coord.x[px]
		DRAW_IMAGE_1 = 1
		break
	end
	'SLY': begin
		WIDGET_CONTROL, event.id, GET_VALUE = pos
		py = pos
		WIDGET_CONTROL, co_y, SET_VALUE = coord.y[py]
		DRAW_IMAGE_2 = 1
		break
	end
	'SLZ': begin
		WIDGET_CONTROL, event.id, GET_VALUE = pos
		pz = pos
		WIDGET_CONTROL, co_z, SET_VALUE = coord.z[pz]
		DRAW_IMAGE_3 = 1
		break
	end
	'DRAW_YZ': begin
		if (event.press) then button_pressed_yz = 1
		if (button_pressed_yz) then begin
			last_py = py
			last_pz = pz
			py = event.x / bin_y > 0 < (num_y-1)
			pz = event.y / bin_z > 0 < (num_z-1)
			if (destretch[1]) then py = pc_find_index (py * (coord.y[num_y-1] - coord.y[0]) / (num_y-1) + coord.y[0], coord.y, num_y)
			if (destretch[2]) then pz = pc_find_index (pz * (coord.z[num_z-1] - coord.z[0]) / (num_z-1) + coord.z[0], coord.z, num_z)
			if ((py ne last_py) or (pz ne last_pz)) then begin
				WIDGET_CONTROL, sl_y, SET_VALUE = py
				WIDGET_CONTROL, sl_z, SET_VALUE = pz
				WIDGET_CONTROL, co_y, SET_VALUE = coord.y[py]
				WIDGET_CONTROL, co_z, SET_VALUE = coord.z[pz]
				DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
			end
		end
		if (event.release) then button_pressed_yz = 0
		break
	end
	'DRAW_XZ': begin
		if (event.press) then button_pressed_xz = 1
		if (button_pressed_xz) then begin
			last_px = px
			last_pz = pz
			px = event.x / bin_x > 0 < (num_x-1)
			pz = event.y / bin_z > 0 < (num_z-1)
			if (destretch[0]) then px = pc_find_index (px * (coord.x[num_x-1] - coord.x[0]) / (num_x-1) + coord.x[0], coord.x, num_x)
			if (destretch[2]) then pz = pc_find_index (pz * (coord.z[num_z-1] - coord.z[0]) / (num_z-1) + coord.z[0], coord.z, num_z)
			if ((px ne last_px) or (pz ne last_pz)) then begin
				WIDGET_CONTROL, sl_x, SET_VALUE = px
				WIDGET_CONTROL, sl_z, SET_VALUE = pz
				WIDGET_CONTROL, co_x, SET_VALUE = coord.x[px]
				WIDGET_CONTROL, co_z, SET_VALUE = coord.z[pz]
				DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
			end
		end
		if (event.release) then button_pressed_xz = 0
		break
	end
	'DRAW_XY': begin
		if (event.press) then button_pressed_xy = 1
		if (button_pressed_xy) then begin
			last_px = px
			last_py = py
			px = event.x / bin_x > 0 < (num_x-1)
			py = event.y / bin_y > 0 < (num_y-1)
			if (destretch[0]) then px = pc_find_index (px * (coord.x[num_x-1] - coord.x[0]) / (num_x-1) + coord.x[0], coord.x, num_x)
			if (destretch[1]) then py = pc_find_index (py * (coord.y[num_y-1] - coord.y[0]) / (num_y-1) + coord.y[0], coord.y, num_y)
			if ((px ne last_px) or (py ne last_py)) then begin
				WIDGET_CONTROL, sl_x, SET_VALUE = px
				WIDGET_CONTROL, sl_y, SET_VALUE = py
				WIDGET_CONTROL, co_x, SET_VALUE = coord.x[px]
				WIDGET_CONTROL, co_y, SET_VALUE = coord.y[py]
				DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
			end
		end
		if (event.release) then button_pressed_xy = 0
		break
	end
	'SCALE_BOT': begin
		WIDGET_CONTROL, sl_min, GET_VALUE = val_min
		if (val_min gt val_max) then begin
			val_min = val_max
			WIDGET_CONTROL, sl_min, SET_VALUE = val_min
		end
		pos_b[selected_cube,sub_aver] = val_min
		DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
		break
	end
	'SCALE_TOP': begin
		WIDGET_CONTROL, sl_max, GET_VALUE = val_max
		if (val_max lt val_min) then begin
			val_max = val_min
			WIDGET_CONTROL, sl_max, SET_VALUE = val_max
		end
		pos_t[selected_cube,sub_aver] = val_max
		DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
		break
	end
	'SCALE_OVER': begin
		WIDGET_CONTROL, sl_over, GET_VALUE = tmp_over
		if (tmp_over gt 1.0) then tmp_over = 1.0 + 1.0 / ((2.0 - tmp_over) > 1.e-8)
		pos_over[selected_overplot] = tmp_over > 1.e-21
		DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
		break
	end
	'MIN_MAX': begin
		cslice_get_minmax_value, cube, val_min, val_max
		WIDGET_CONTROL, sl_min, SET_VALUE = val_min
		WIDGET_CONTROL, sl_max, SET_VALUE = val_max
		pos_b[selected_cube,sub_aver] = val_min
		pos_t[selected_cube,sub_aver] = val_max
		DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
		break
	end
	'JUMP_MIN': begin
		pos_min = min (cube, location)
		dims = size (cube, /dimensions)
		pos_min = array_indices (dims, location, /dimensions)
		px = pos_min[0]
		py = pos_min[1]
		pz = pos_min[2]
		WIDGET_CONTROL, sl_x, SET_VALUE = px
		WIDGET_CONTROL, sl_y, SET_VALUE = py
		WIDGET_CONTROL, sl_z, SET_VALUE = pz
		WIDGET_CONTROL, co_x, SET_VALUE = coord.x[px]
		WIDGET_CONTROL, co_y, SET_VALUE = coord.y[py]
		WIDGET_CONTROL, co_z, SET_VALUE = coord.z[pz]
		DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
		break
	end
	'JUMP_MAX': begin
		pos_max = max (cube, location)
		dims = size (cube, /dimensions)
		pos_max = array_indices (dims, location, /dimensions)
		px = pos_max[0]
		py = pos_max[1]
		pz = pos_max[2]
		WIDGET_CONTROL, sl_x, SET_VALUE = px
		WIDGET_CONTROL, sl_y, SET_VALUE = py
		WIDGET_CONTROL, sl_z, SET_VALUE = pz
		WIDGET_CONTROL, co_x, SET_VALUE = coord.x[px]
		WIDGET_CONTROL, co_y, SET_VALUE = coord.y[py]
		WIDGET_CONTROL, co_z, SET_VALUE = coord.z[pz]
		DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
		break
	end
	'FREEZE': begin
		WIDGET_CONTROL, sl_min, SENSITIVE = 0
		WIDGET_CONTROL, sl_max, SENSITIVE = 0
		WIDGET_CONTROL, sl_over, SENSITIVE = 0
		WIDGET_CONTROL, freeze, set_value='RELEASE RANGE', set_uvalue='RELEASE'
		pos_b[selected_cube,sub_aver] = val_min
		pos_t[selected_cube,sub_aver] = val_max
		frozen = 1
		break
	end
	'RELEASE': begin
		WIDGET_CONTROL, sl_min, SENSITIVE = 1
		WIDGET_CONTROL, sl_max, SENSITIVE = 1
		WIDGET_CONTROL, sl_over, SENSITIVE = (selected_overplot gt 0)
		WIDGET_CONTROL, freeze, set_value='FREEZE RANGE', set_uvalue='FREEZE'
		frozen = 0
		cslice_prepare_cube, -1
		DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
		break
	end
	'VAR': begin
		if (selected_cube ne event.index) then begin
			cslice_prepare_cube, event.index
			DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
		end
		break
	end
	'SNAP': begin
		if (selected_snapshot ne event.index) then begin
			cslice_prepare_set, event.index
			cslice_prepare_cube, -1
			cslice_prepare_overplot
			DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
			WIDGET_CONTROL, prev, SENSITIVE = (selected_snapshot lt num_snapshots - 1)
			WIDGET_CONTROL, next, SENSITIVE = (selected_snapshot gt 0)

			window, 0, xsize=8, ysize=8, retain=2
			!P.MULTI = [0, 1, 1]
			wdelete
		end
		break
	end
	'NEXT': begin
		selected_snapshot -= 1
		if (selected_snapshot le 0) then begin
			WIDGET_CONTROL, next, SENSITIVE = 0
			selected_snapshot = 0
		end
		WIDGET_CONTROL, prev, SENSITIVE = 1
		WIDGET_CONTROL, snap, SET_DROPLIST_SELECT = selected_snapshot
		cslice_prepare_set, selected_snapshot
		cslice_prepare_cube, -1
		cslice_prepare_overplot
		DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
		break
	end
	'PREV': begin
		selected_snapshot += 1
		if (selected_snapshot ge num_snapshots - 1) then begin
			WIDGET_CONTROL, prev, SENSITIVE = 0
			selected_snapshot = num_snapshots - 1
		end
		WIDGET_CONTROL, next, SENSITIVE = 1
		WIDGET_CONTROL, snap, SET_DROPLIST_SELECT = selected_snapshot
		cslice_prepare_set, selected_snapshot
		cslice_prepare_cube, -1
		cslice_prepare_overplot
		DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
		break
	end
	'OVER': begin
		if (selected_overplot ne event.index) then begin
			selected_overplot = event.index
			cslice_prepare_overplot
			DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
		end
		break
	end
	'RESET': begin
		cslice_reset_GUI
		DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
		break
	end
	'LOAD': begin
		if (not file_test (settings_file, /read)) then begin
			WIDGET_CONTROL, b_load, SENSITIVE = 0
			break
		end
		restore, settings_file
		selected_var = (where (tag_names (varsets) eq var_names[selected_var]))[0]
		if (selected_var ge 0) then selected_cube = selected_var
		num = n_elements (var_names)
		for pos = 0, num-1 do begin
			insert = (where (tag_names (varsets) eq var_names[pos]))[0]
			if (insert ge 0) then begin
				pos_b[insert] = pos_bot[pos]
				pos_t[insert] = pos_top[pos]
			end
		end
		selected_over = (where (tag_names (overplot) eq var_names[selected_over]))[0]
		if (selected_over ge 0) then selected_overplot = selected_over
		num = n_elements (over_names)
		for pos = 0, num-1 do begin
			insert = (where (tag_names (overplot) eq over_names[pos]))[0]
			if (insert ge 0) then pos_over[insert] = pos_overplot[pos]
		end
		if (px gt (num_x - 1)) then px = num_x - 1
		if (py gt (num_y - 1)) then py = num_y - 1
		if (pz gt (num_z - 1)) then pz = num_z - 1
		cslice_prepare_cube, -1
		cslice_prepare_overplot
		DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
		WIDGET_CONTROL, b_log, SET_VALUE = log_plot
		WIDGET_CONTROL, b_abs, SET_VALUE = abs_scale
		WIDGET_CONTROL, b_sub, SET_VALUE = sub_aver
		WIDGET_CONTROL, b_cro, SET_VALUE = show_cross
		if (any (coord.lequidist eq 0)) then WIDGET_CONTROL, b_des, SET_VALUE = any (destretch)
		WIDGET_CONTROL, co_x, SET_VALUE = varfiles[selected_snapshot].time
		WIDGET_CONTROL, sl_x, SET_VALUE = selected_snapshot
		WIDGET_CONTROL, co_x, SET_VALUE = coord.x[px]
		WIDGET_CONTROL, co_y, SET_VALUE = coord.y[py]
		WIDGET_CONTROL, co_z, SET_VALUE = coord.z[pz]
		WIDGET_CONTROL, sl_x, SET_VALUE = px
		WIDGET_CONTROL, sl_y, SET_VALUE = py
		WIDGET_CONTROL, sl_z, SET_VALUE = pz
		WIDGET_CONTROL, sl_min, SET_VALUE = pos_b[selected_cube,sub_aver]
		WIDGET_CONTROL, sl_max, SET_VALUE = pos_t[selected_cube,sub_aver]
		WIDGET_CONTROL, sl_over, SET_VALUE = [ pos_over[selected_overplot], 0.0, 2.0 ]
		WIDGET_CONTROL, vars, SET_DROPLIST_SELECT = selected_cube
		WIDGET_CONTROL, over, SET_DROPLIST_SELECT = selected_overplot
		WIDGET_CONTROL, snap, SET_DROPLIST_SELECT = selected_snapshot
		break
	end
	'SAVE': begin
		var_names = tag_names (varsets)
		over_names = tag_names (overplot)
		pos_bot = pos_b
		pos_top = pos_t
		pos_overplot = pos_over
		pos_bot[selected_cube,sub_aver] = val_min
		pos_top[selected_cube,sub_aver] = val_max
		selected_var = selected_cube
		selected_over = selected_overplot
		save, filename=settings_file, var_names, over_names, selected_var, selected_over, log_plot, abs_scale, sub_aver, show_cross, destretch, px, py, pz, pos_bot, pos_top, pos_overplot
		WIDGET_CONTROL, b_load, SENSITIVE = 1
		break
	end
	'PLAY':
	'IMAGE': begin
		WIDGET_CONTROL, vars, SENSITIVE = 0
		WIDGET_CONTROL, over, SENSITIVE = 0
		WIDGET_CONTROL, snap, SENSITIVE = 0
		WIDGET_CONTROL, prev, SENSITIVE = 0
		WIDGET_CONTROL, next, SENSITIVE = 0
		WIDGET_CONTROL, co_t, SENSITIVE = 1
		WIDGET_CONTROL, sl_t, SENSITIVE = 0
		WIDGET_CONTROL, play, SENSITIVE = 0
		WIDGET_CONTROL, image, SENSITIVE = 0
		WIDGET_CONTROL, aver, SENSITIVE = 0
		WIDGET_CONTROL, timeser, SENSITIVE = 0
		WIDGET_CONTROL, min_max, SENSITIVE = 0
		WIDGET_CONTROL, freeze, SENSITIVE = 0
		WIDGET_CONTROL, sl_min, SENSITIVE = 0
		WIDGET_CONTROL, sl_max, SENSITIVE = 0
		WIDGET_CONTROL, sl_over, SENSITIVE = 0
		WIDGET_CONTROL, jump_min, SENSITIVE = 0
		WIDGET_CONTROL, jump_max, SENSITIVE = 0
		orig_frozen = frozen
		frozen = 1
		if (num_snapshots gt 1) then begin
			if (eventval eq "IMAGE") then cslice_save_movie, num_snapshots, /allocate
			previous_snapshot = selected_snapshot
			for i = num_snapshots-1, 0, -1 do begin
				cslice_prepare_set, i
				cslice_prepare_cube, -1
				cslice_prepare_overplot
				cslice_draw, 1, 1, 1
				if (eventval eq "IMAGE") then begin
					cslice_save_images, "PNG", movie_frame=i
				end else begin
					wait, min_wait_time[dimensionality-1]
				end
			end
			if (previous_snapshot ne 0) then begin
				cslice_prepare_set, previous_snapshot
				cslice_prepare_cube, -1
				cslice_prepare_overplot
			end
		end else begin
			if (eventval eq "IMAGE") then cslice_save_images, "PNG", /slices
		end
		frozen = orig_frozen
		WIDGET_CONTROL, vars, SENSITIVE = (num_cubes ge 2)
		WIDGET_CONTROL, over, SENSITIVE = (num_overs ge 2)
		WIDGET_CONTROL, snap, SENSITIVE = (num_snapshots ge 2)
		WIDGET_CONTROL, prev, SENSITIVE = (selected_snapshot lt num_snapshots - 1)
		WIDGET_CONTROL, next, SENSITIVE = (selected_snapshot gt 0)
		WIDGET_CONTROL, co_t, SENSITIVE = (num_snapshots ge 2)
		WIDGET_CONTROL, sl_t, SENSITIVE = (num_snapshots ge 2)
		WIDGET_CONTROL, play, SENSITIVE = (num_snapshots ge 2)
		WIDGET_CONTROL, image, SENSITIVE = 1
		WIDGET_CONTROL, aver, SENSITIVE = 1
		WIDGET_CONTROL, timeser, SENSITIVE = 1
		WIDGET_CONTROL, min_max, SENSITIVE = 1
		WIDGET_CONTROL, freeze, SENSITIVE = 1
		if (not frozen) then begin
			WIDGET_CONTROL, sl_min, SENSITIVE = 1
			WIDGET_CONTROL, sl_max, SENSITIVE = 1
			WIDGET_CONTROL, sl_over, SENSITIVE = (selected_overplot gt 0)
		end
		WIDGET_CONTROL, jump_min, SENSITIVE = 1
		WIDGET_CONTROL, jump_max, SENSITIVE = 1

		if (show_cuts) then begin
			window, 0, xsize=8, ysize=8, retain=2
			!P.MULTI = [0, 1, 1]
			wdelete
		end
		break
	end
	'QUIT': begin
		quit = event.top
		break
	end
	endswitch

	cslice_draw, DRAW_IMAGE_1, DRAW_IMAGE_2, DRAW_IMAGE_3

	WIDGET_CONTROL, WIDGET_INFO (event.top, /CHILD)

	if (quit ge 0) then WIDGET_CONTROL, quit, /DESTROY

	return
end


; Draws the slices into the window
pro cslice_draw, DRAW_IMAGE_1, DRAW_IMAGE_2, DRAW_IMAGE_3

	common varset_common, set, overplot, oversets, unit, coord, varsets, varfiles, datadir, sources, param, run_param, var_list
	common cslice_common, cube, field, num_cubes, num_overs, num_snapshots
	common overplot_common, overplot_contour, field_x_y, field_x_z, field_y_x, field_y_z, field_z_x, field_z_y, field_x_indices, field_y_indices, field_z_indices, vector_distance, vector_length
	common slider_common, bin_x, bin_y, bin_z, num_x, num_y, num_z, pos_b, pos_t, pos_over, val_min, val_max, val_range, over_max, dimensionality, frozen
	common gui_common, wimg_yz, wimg_xz, wimg_xy, wcut_x, wcut_y, wcut_z, co_t, sl_t, co_x, co_y, co_z, sl_x, sl_y, sl_z, b_load, b_log, b_sub, b_abs, b_cro, b_des, aver, timeser, vars, over, snap, prev, next, play, image, sl_min, sl_max, sl_over, min_max, freeze, jump_min, jump_max
	common settings_common, px, py, pz, cut, log_plot, abs_scale, show_cross, show_cuts, sub_aver, selected_cube, selected_overplot, selected_snapshot, af_x, af_y, af_z, destretch

	; stepping of crosshairs
	step = 4

	; number of levels for contour plot
	default, nlevels, 50

	; default target is the screen (=0)
	default, target, 0

	!P.MULTI = [0, 1, 1]

	ox = floor (bin_x / 2.0)
	oy = floor (bin_y / 2.0)
	oz = floor (bin_z / 2.0)
	if (destretch[0]) then cx = pc_find_index (coord.x[px], pc_destretch (coord.x, coord.x, dim=1), num_x) else cx = px
	if (destretch[1]) then cy = pc_find_index (coord.y[py], pc_destretch (coord.y, coord.y, dim=1), num_y) else cy = py
	if (destretch[2]) then cz = pc_find_index (coord.z[pz], pc_destretch (coord.z, coord.z, dim=1), num_z) else cz = pz

	num_over_x = n_elements (field_x_indices)
	num_over_y = n_elements (field_y_indices)
	num_over_z = n_elements (field_z_indices)

	over_length = pos_over[selected_overplot] > 1.e-42
	local_max = over_max

	if (DRAW_IMAGE_1 or DRAW_IMAGE_2 or DRAW_IMAGE_3) then begin
		data = reform ((cube[px,*,*] > val_min) < val_max, num_y, num_z)
		if (destretch[1]) then data = pc_destretch (data, coord.y, dim=1)
		if (destretch[2]) then data = pc_destretch (data, coord.z, dim=2)
		if ((bin_y ne 1) or (bin_z ne 1)) then data = congrid (data, fix (num_y*bin_y), fix (num_z*bin_z), cubic=0)
		val_range = cslice_get_range (data)
		wset, wimg_yz
		if (show_cross) then begin
			if (cy gt af_y) then for i = fix ((cy-af_y)*bin_y), 0, -step do data[i:i+1,cz*bin_z+oz] = val_range
			if (cy lt num_y-1-af_y) then for i = fix ((cy+af_y)*bin_y), fix (num_y*bin_y)-2, step do data[i:i+1,cz*bin_z+oz] = val_range
			if (cz gt af_z) then for i = fix ((cz-af_z)*bin_z), 0, -step do data[cy*bin_y+oy,i:i+1] = val_range
			if (cz lt num_z-1-af_z) then for i = fix ((cz+af_z)*bin_z), fix (num_z*bin_z)-2, step do data[cy*bin_y+oy,i:i+1] = val_range
			if ((cy le af_y) and (cy ge num_y-1-af_y) and (cz le af_z) and (cz ge num_z-1-af_z)) then data[0:1,0] = val_range
		end else if (abs_scale) then begin
			if (min (data) gt val_min) then data[0,0] = val_min
			if (max (data) lt val_max) then data[1,0] = val_max
		end
		if (log_plot) then data = alog10 (data)
		tvscl, data
		if (selected_overplot gt 0) then begin
			if (overplot_contour eq 1) then begin
				if ((n_elements (field_y_indices) gt 1) and (n_elements (field_z_indices) gt 1)) then begin
					contour, reform (field_x_y[cx,*,*], num_over_y, num_over_z), field_y_indices, field_z_indices, nlevels=nlevels, xs=4, ys=4, color=200, /noerase, pos=[0.0,0.0,1.0,1.0]
				endif
			end else begin
				if (abs_scale) then local_max = max (abs ([field_y_x[cx,*,*],field_z_x[cx,*,*]]))
				vec_len = vector_length * over_length * local_max / over_max
				velovect, reform (field_y_x[cx,*,*], num_over_y, num_over_z), reform (field_z_x[cx,*,*], num_over_y, num_over_z), field_y_indices, field_z_indices, length=vec_len, xr=[0.0,1.0], yr=[0.0,1.0], xs=4, ys=4, color=200, /noerase, pos=[0.0,0.0,1.0,1.0]
			end
		end
		if (show_cuts and (DRAW_IMAGE_1 or DRAW_IMAGE_3)) then begin
			wset, wcut_y
			data = reform (cube[px,*,pz])
			if (any (destretch)) then begin
				plot, coord.y, data, xrange=minmax (coord.y), yrange=val_range, xstyle=1, ystyle=1, xmargin=[0,0], ymargin=[0,0], ylog=log_plot
				oplot, coord.y, data, psym=3, color=200
			end else begin
				plot, data, xrange=[0,num_y-1], yrange=val_range, xstyle=1, ystyle=1, xmargin=[0,0], ymargin=[0,0], ylog=log_plot
				oplot, data, psym=3, color=200
			end
			axis, 0, 0, xaxis=1, xstyle=1, ystyle=1
			axis, 0, 0, yaxis=1, xstyle=1, ystyle=1
		end
	end

	if (DRAW_IMAGE_1 or DRAW_IMAGE_2 or DRAW_IMAGE_3) then begin
		data = reform ((cube[*,py,*] > val_min) < val_max, num_x, num_z)
		if (destretch[0]) then data = pc_destretch (data, coord.x, dim=1)
		if (destretch[2]) then data = pc_destretch (data, coord.z, dim=2)
		if ((bin_x ne 1) or (bin_z ne 1)) then data = congrid (data, fix (num_x*bin_x), fix (num_z*bin_z), cubic=0)
		val_range = cslice_get_range (data)
		wset, wimg_xz
		if (show_cross) then begin
			if (cx gt af_x) then for i = fix ((cx-af_x)*bin_x), 0, -step do data[i:i+1,cz*bin_z+oz] = val_range
			if (cx lt num_x-1-af_x) then for i = fix ((cx+af_x)*bin_x), fix (num_x*bin_x)-2, step do data[i:i+1,cz*bin_z+oz] = val_range
			if (cz gt af_z) then for i = fix ((cz-af_z)*bin_z), 0, -step do data[cx*bin_x+ox,i:i+1] = val_range
			if (cz lt num_z-1-af_z) then for i = fix ((cz+af_z)*bin_z), fix (num_z*bin_z)-2, step do data[cx*bin_x+ox,i:i+1] = val_range
			if ((cx le af_x) and (cx ge num_x-1-af_x) and (cz le af_z) and (cz ge num_z-1-af_z)) then data[0:1,0] = val_range
		end else if (abs_scale) then begin
			if (min (data) gt val_min) then data[0,0] = val_min
			if (max (data) lt val_max) then data[1,0] = val_max
		end
		if (log_plot) then data = alog10 (data)
		tvscl, data
		if (selected_overplot gt 0) then begin
			if (overplot_contour eq 1) then begin
				if ((n_elements (field_x_indices) gt 1) and (n_elements (field_z_indices) gt 1)) then begin
					contour, reform (field_y_x[*,cy,*], num_over_x, num_over_z), field_x_indices, field_z_indices, nlevels=nlevels, xs=4, ys=4, color=200, /noerase, pos=[0.0,0.0,1.0,1.0]
				endif
			end else begin
				if (abs_scale) then local_max = max (abs ([field_x_y[*,cy,*],field_z_y[*,cy,*]]))
				vec_len = vector_length * over_length * local_max / over_max
				velovect, reform (field_x_y[*,cy,*], num_over_x, num_over_z), reform (field_z_y[*,cy,*], num_over_x, num_over_z), field_x_indices, field_z_indices, length=vec_len, xr=[0.0,1.0], yr=[0.0,1.0], xs=4, ys=4, color=200, /noerase, pos=[0.0,0.0,1.0,1.0]
			end
		end
		if (show_cuts and (DRAW_IMAGE_2 or DRAW_IMAGE_3)) then begin
			wset, wcut_x
			data = reform (cube[*,py,pz])
			if (any (destretch)) then begin
				plot, coord.x, data, xrange=minmax (coord.x), yrange=val_range, xstyle=1, ystyle=1, xmargin=[0,0], ymargin=[0,0], ylog=log_plot
				oplot, coord.x, data, psym=3, color=200
			end else begin
				plot, data, xrange=[0,num_x-1], yrange=val_range, xstyle=1, ystyle=1, xmargin=[0,0], ymargin=[0,0], ylog=log_plot
				oplot, data, psym=3, color=200
			end
			axis, 0, 0, xaxis=1, xstyle=1, ystyle=1
			axis, 0, 0, yaxis=1, xstyle=1, ystyle=1
		end
	end

	if (DRAW_IMAGE_1 or DRAW_IMAGE_2 or DRAW_IMAGE_3) then begin
		data = reform ((cube[*,*,pz] > val_min) < val_max, num_x, num_y)
		if (destretch[0]) then data = pc_destretch (data, coord.x, dim=1)
		if (destretch[1]) then data = pc_destretch (data, coord.y, dim=2)
		if ((bin_x ne 1) or (bin_y ne 1)) then data = congrid (data, fix (num_x*bin_x), fix (num_y*bin_y), cubic=0)
		val_range = cslice_get_range (data)
		wset, wimg_xy
		if (show_cross) then begin
			if (cx gt af_x) then for i = fix ((cx-af_x)*bin_x), 0, -step do data[i:i+1,cy*bin_y+oy] = val_range
			if (cx lt num_x-1-af_x) then for i = fix ((cx+af_x)*bin_x), fix (num_x*bin_x)-2, step do data[i:i+1,cy*bin_y+oy] = val_range
			if (cy gt af_y) then for i = fix ((cy-af_y)*bin_y), 0, -step do data[cx*bin_x+ox,i:i+1] = val_range
			if (cy lt num_y-1-af_y) then for i = fix ((cy+af_y)*bin_y), fix (num_y*bin_y)-2, step do data[cx*bin_x+ox,i:i+1] = val_range
			if ((cx le af_x) and (cx ge num_x-1-af_x) and (cy le af_y) and (cy ge num_y-1-af_y)) then data[0:1,0] = val_range
		end else if (abs_scale) then begin
			if (min (data) gt val_min) then data[0,0] = val_min
			if (max (data) lt val_max) then data[1,0] = val_max
		end
		if (log_plot) then data = alog10 (data)
		tvscl, data
		if (selected_overplot gt 0) then begin
			if (overplot_contour eq 1) then begin
				if ((n_elements (field_x_indices) gt 1) and (n_elements (field_y_indices) gt 1)) then begin
					contour, reform (field_z_x[*,*,cz], num_over_x, num_over_y), field_x_indices, field_y_indices, nlevels=nlevels, xs=4, ys=4, color=200, /noerase, pos=[0.0,0.0,1.0,1.0]
				endif
			end else begin
				if (abs_scale) then local_max = max (abs ([field_x_z[*,*,cz],field_y_z[*,*,cz]]))
				vec_len = vector_length * over_length * local_max / over_max
				velovect, reform (field_x_z[*,*,cz], num_over_x, num_over_y), reform (field_y_z[*,*,cz], num_over_x, num_over_y), field_x_indices, field_y_indices, length=vec_len, xr=[0.0,1.0], yr=[0.0,1.0], xs=4, ys=4, color=200, /noerase, pos=[0.0,0.0,1.0,1.0]
			end
		end
		if (show_cuts and (DRAW_IMAGE_1 or DRAW_IMAGE_2)) then begin
			wset, wcut_z
			data = reform (cube[px,py,*])
			if (any (destretch)) then begin
				plot, coord.z, data, xrange=minmax (coord.z), yrange=val_range, xstyle=1, ystyle=1, xmargin=[0,0], ymargin=[0,0], ylog=log_plot
				oplot, coord.z, data, psym=3, color=200
			end else begin
				plot, data, xrange=[0,num_z-1], yrange=val_range, xstyle=1, ystyle=1, xmargin=[0,0], ymargin=[0,0], ylog=log_plot
				oplot, data, psym=3, color=200
			end
			axis, 0, 0, xaxis=1, xstyle=1, ystyle=1
			axis, 0, 0, yaxis=1, xstyle=1, ystyle=1
		end
	end
end


; Saves images of the slices with the given format
pro cslice_save_images, img_type, slices=slices, movie_frame=movie_frame

	common varset_common, set, overplot, oversets, unit, coord, varsets, varfiles, datadir, sources, param, run_param, var_list
	common gui_common, wimg_yz, wimg_xz, wimg_xy, wcut_x, wcut_y, wcut_z, co_t, sl_t, co_x, co_y, co_z, sl_x, sl_y, sl_z, b_load, b_log, b_sub, b_abs, b_cro, b_des, aver, timeser, vars, over, snap, prev, next, play, image, sl_min, sl_max, sl_over, min_max, freeze, jump_min, jump_max
	common settings_common, px, py, pz, cut, log_plot, abs_scale, show_cross, show_cuts, sub_aver, selected_cube, selected_overplot, selected_snapshot, af_x, af_y, af_z, destretch

	quantity = (tag_names (set))[selected_cube]
	prefix = varfiles[selected_snapshot].title + "_" + quantity
	if (selected_overplot gt 0) then prefix += "_overplot-" + (tag_names (overplot))[selected_overplot]
	suffix = "." + strlowcase (img_type)

	pc_save_image, prefix+"_xy"+suffix, window=wimg_xy
	pc_save_image, prefix+"_xz"+suffix, window=wimg_xz
	pc_save_image, prefix+"_yz"+suffix, window=wimg_yz
	pc_save_image, prefix+"_x"+suffix, window=wcut_x
	pc_save_image, prefix+"_y"+suffix, window=wcut_y
	pc_save_image, prefix+"_z"+suffix, window=wcut_z

	if (keyword_set (slices)) then cslice_save_slices
	if (n_elements (movie_frame)) then cslice_save_movie, movie_frame
end


; Saves the slices into a XDR file
pro cslice_save_slices

	common cslice_common, cube, field, num_cubes, num_overs, num_snapshots
	common varset_common, set, overplot, oversets, unit, coord, varsets, varfiles, datadir, sources, param, run_param, var_list
	common slider_common, bin_x, bin_y, bin_z, num_x, num_y, num_z, pos_b, pos_t, pos_over, val_min, val_max, val_range, over_max, dimensionality, frozen
	common settings_common, px, py, pz, cut, log_plot, abs_scale, show_cross, show_cuts, sub_aver, selected_cube, selected_overplot, selected_snapshot, af_x, af_y, af_z, destretch

        quantity = (tag_names (set))[selected_cube]
	prefix = varfiles[selected_snapshot].title + "_" + quantity

	cut_xy = reform (cube[*,*,pz], num_x, num_y)
	cut_xz = reform (cube[*,py,*], num_x, num_z)
	cut_yz = reform (cube[px,*,*], num_y, num_z)
	cut_x = reform (cube[*,py,pz])
	cut_y = reform (cube[px,*,pz])
	cut_z = reform (cube[px,py,*])
	x = coord.x * unit.default_length
	y = coord.y * unit.default_length
	z = coord.z * unit.default_length
	dx = coord.dx
	dy = coord.dy
	dz = coord.dz
	pos_x = px
	pos_y = py
	pos_z = pz
	time = varfiles[selected_snapshot].time

	save, filename=prefix+"_cuts.xdr", time, quantity, cut_xy, cut_xz, cut_yz, cut_x, cut_y, cut_z, pos_x, pos_y, pos_z, x, y, z, dx, dy, dz
end


; Saves a movie of slices into a XDR file
pro cslice_save_movie, frame, allocate=allocate

	common cslice_movie_common, cut_xy, cut_xz, cut_yz, cut_x, cut_y, cut_z, time
	common cslice_common, cube, field, num_cubes, num_overs, num_snapshots
	common varset_common, set, overplot, oversets, unit, coord, varsets, varfiles, datadir, sources, param, run_param, var_list
	common slider_common, bin_x, bin_y, bin_z, num_x, num_y, num_z, pos_b, pos_t, pos_over, val_min, val_max, val_range, over_max, dimensionality, frozen
	common settings_common, px, py, pz, cut, log_plot, abs_scale, show_cross, show_cuts, sub_aver, selected_cube, selected_overplot, selected_snapshot, af_x, af_y, af_z, destretch

	if (keyword_set (allocate)) then begin
		cut_xy = dblarr (num_x, num_y, frame)
		cut_xz = dblarr (num_x, num_z, frame)
		cut_yz = dblarr (num_y, num_z, frame)
		cut_x = dblarr (num_x, frame)
		cut_y = dblarr (num_y, frame)
		cut_z = dblarr (num_z, frame)
		time = dblarr (frame)
		time[*] = !Values.D_NaN
		return
	end

	cut_xy[*,*,frame] = reform (cube[*,*,pz], num_x, num_y)
	cut_xz[*,*,frame] = reform (cube[*,py,*], num_x, num_z)
	cut_yz[*,*,frame] = reform (cube[px,*,*], num_y, num_z)
	cut_x[*,frame] = reform (cube[*,py,pz])
	cut_y[*,frame] = reform (cube[px,*,pz])
	cut_z[*,frame] = reform (cube[px,py,*])
	time[frame] = varfiles[selected_snapshot].time
	if (any (finite (time, /NaN))) then return

        quantity = (tag_names (set))[selected_cube]
	x = coord.x * unit.default_length
	y = coord.y * unit.default_length
	z = coord.z * unit.default_length
	dx = coord.dx
	dy = coord.dy
	dz = coord.dz
	pos_x = px
	pos_y = py
	pos_z = pz

	save, filename=quantity+"_movie.xdr", time, quantity, cut_xy, cut_xz, cut_yz, cut_x, cut_y, cut_z, pos_x, pos_y, pos_z, x, y, z, dx, dy, dz

	cut_xy = 0
	cut_xz = 0
	cut_yz = 0
	cut_x = 0
	cut_y = 0
	cut_z = 0
	time = 0
end


; Draws horizontally averaged vertical profiles into a second window
pro cslice_draw_averages, number

	common varset_common, set, overplot, oversets, unit, coord, varsets, varfiles, datadir, sources, param, run_param, var_list
	common settings_common, px, py, pz, cut, log_plot, abs_scale, show_cross, show_cuts, sub_aver, selected_cube, selected_overplot, selected_snapshot, af_x, af_y, af_z, destretch
	common slider_common, bin_x, bin_y, bin_z, num_x, num_y, num_z, pos_b, pos_t, pos_over, val_min, val_max, val_range, over_max, dimensionality, frozen

	vert_label = "height"
	if (unit.default_length_str) then begin
		vert_label += ' ['+unit.default_length_str+']'
	end else if (strupcase (param.unit_system) eq "SI") then begin
		vert_label += ' [m]'
	end else if (strupcase (param.unit_system) eq "CGS") then begin
		vert_label += ' [cm]'
	end else begin
		vert_label += ' [code units]'
	end

	prefix = varfiles[selected_snapshot].title + "_" + (tag_names (set))[selected_cube]
	time = strtrim (varfiles[selected_snapshot].time * unit.time/unit.default_time, 2) + " " + unit.default_time_str
	pc_vert_profile, reform (varsets[number].(selected_cube)[cut], num_x, num_y, num_z), coord=coord.z, title=set.(selected_cube), log=log_plot, horiz_label='['+param.unit_system+']', vert_label=vert_label, file_label=prefix, time=time
end


; Prepares a data set for visualisation
pro cslice_prepare_set, i

	common varset_common, set, overplot, oversets, unit, coord, varsets, varfiles, datadir, sources, param, run_param, var_list
	common cslice_common, cube, field, num_cubes, num_overs, num_snapshots
	common gui_common, wimg_yz, wimg_xz, wimg_xy, wcut_x, wcut_y, wcut_z, co_t, sl_t, co_x, co_y, co_z, sl_x, sl_y, sl_z, b_load, b_log, b_sub, b_abs, b_cro, b_des, aver, timeser, vars, over, snap, prev, next, play, image, sl_min, sl_max, sl_over, min_max, freeze, jump_min, jump_max
	common settings_common, px, py, pz, cut, log_plot, abs_scale, show_cross, show_cuts, sub_aver, selected_cube, selected_overplot, selected_snapshot, af_x, af_y, af_z, destretch

	selected_snapshot = i

	if ((selected_snapshot ge 0) and not finite (co_t, /NaN) and not finite (sl_t, /NaN)) then begin
		WIDGET_CONTROL, co_t, SET_VALUE = varfiles[selected_snapshot].time
		WIDGET_CONTROL, sl_t, SET_VALUE = num_snapshots-1-selected_snapshot
	end

	pc_gui_precalc, i

	num_cubes = n_tags (set)

	return
end


; Prepares a cube for visualisation
pro cslice_prepare_cube, cube_index

	common varset_common, set, overplot, oversets, unit, coord, varsets, varfiles, datadir, sources, param, run_param, var_list
	common cslice_common, cube, field, num_cubes, num_overs, num_snapshots
	common slider_common, bin_x, bin_y, bin_z, num_x, num_y, num_z, pos_b, pos_t, pos_over, val_min, val_max, val_range, over_max, dimensionality, frozen
	common gui_common, wimg_yz, wimg_xz, wimg_xy, wcut_x, wcut_y, wcut_z, co_t, sl_t, co_x, co_y, co_z, sl_x, sl_y, sl_z, b_load, b_log, b_sub, b_abs, b_cro, b_des, aver, timeser, vars, over, snap, prev, next, play, image, sl_min, sl_max, sl_over, min_max, freeze, jump_min, jump_max
	common settings_common, px, py, pz, cut, log_plot, abs_scale, show_cross, show_cuts, sub_aver, selected_cube, selected_overplot, selected_snapshot, af_x, af_y, af_z, destretch

	; SETTINGS:
	; fraction of box width to keep free of crosshairs at center
	af_fraction = 1.0 / 8.0
	; minimum size of crosshairs
	af_minimum = 6
	af_maximum = 32

	; get selected cube from set
	if (cube_index ge 0) then selected_cube = cube_index
	cube = reform (varsets[selected_snapshot].(selected_cube)[cut], num_x, num_y, num_z)

	; substract horizontal averages
	if (sub_aver) then for z=0, num_z-1 do cube[*,*,z] -= mean (cube [*,*,z])

	; find minimum and maximum values
	cslice_get_minmax_value, cube, tmp_min, tmp_max

	if (frozen) then begin
		; update slider to intermediate min/max values
		WIDGET_CONTROL, sl_min, SET_VALUE = [ val_min, (tmp_min < val_min), tmp_max ]
		WIDGET_CONTROL, sl_max, SET_VALUE = [ val_max, tmp_min, (tmp_max > val_max) ]
	end else begin
		; find minimum and maximum values
		val_min = tmp_min
		val_max = tmp_max
		val_range = cslice_get_range ([tmp_min, tmp_max])

		; set default slider positions (min/max)
		if (finite (pos_b[selected_cube,sub_aver], /NaN)) then pos_b[selected_cube,sub_aver] = val_min
		if (finite (pos_t[selected_cube,sub_aver], /NaN)) then pos_t[selected_cube,sub_aver] = val_max

		; adjust slider positions to fit inside value range
		if (pos_b[selected_cube,sub_aver] lt val_min) then pos_b[selected_cube,sub_aver] = val_min
		if (pos_t[selected_cube,sub_aver] gt val_max) then pos_t[selected_cube,sub_aver] = val_max

		; update slider
		if (sl_min ne 0) then WIDGET_CONTROL, sl_min, SET_VALUE = [ pos_b[selected_cube,sub_aver], val_range ]
		if (sl_max ne 0) then WIDGET_CONTROL, sl_max, SET_VALUE = [ pos_t[selected_cube,sub_aver], val_range ]

		; get desired min/max values from sliders
		val_min = pos_b[selected_cube,sub_aver]
		val_max = pos_t[selected_cube,sub_aver]
	end

	; setup crosshairs parameters
	af_x = (round (num_x * af_fraction) > af_minimum) < af_maximum
	af_y = (round (num_y * af_fraction) > af_minimum) < af_maximum
	af_z = (round (num_z * af_fraction) > af_minimum) < af_maximum

	return
end


; Prepares an overplot for visualisation
pro cslice_prepare_overplot

	common varset_common, set, overplot, oversets, unit, coord, varsets, varfiles, datadir, sources, param, run_param, var_list
	common cslice_common, cube, field, num_cubes, num_overs, num_snapshots
	common overplot_common, overplot_contour, field_x_y, field_x_z, field_y_x, field_y_z, field_z_x, field_z_y, field_x_indices, field_y_indices, field_z_indices, vector_distance, vector_length
	common slider_common, bin_x, bin_y, bin_z, num_x, num_y, num_z, pos_b, pos_t, pos_over, val_min, val_max, val_range, over_max, dimensionality, frozen
	common gui_common, wimg_yz, wimg_xz, wimg_xy, wcut_x, wcut_y, wcut_z, co_t, sl_t, co_x, co_y, co_z, sl_x, sl_y, sl_z, b_load, b_log, b_sub, b_abs, b_cro, b_des, aver, timeser, vars, over, snap, prev, next, play, image, sl_min, sl_max, sl_over, min_max, freeze, jump_min, jump_max
	common settings_common, px, py, pz, cut, log_plot, abs_scale, show_cross, show_cuts, sub_aver, selected_cube, selected_overplot, selected_snapshot, af_x, af_y, af_z, destretch

	; SETTINGS:
	; distance of vector footpoint locations
	vector_distance = 8
	; maximum length of vectors
	vector_length = vector_distance * sqrt (2)
	; default plot routine: 0=velovect (1=contour)
	overplot_contour = 0

	if (strpos (strlowcase ((tag_names (overplot))[selected_overplot]), "_contour") gt 0) then overplot_contour = 1
	if ((selected_overplot le 0) or overplot_contour) then begin
		WIDGET_CONTROL, sl_over, SET_VALUE = [ 1.0, 0.0, 2.0 ]
		WIDGET_CONTROL, sl_over, SENSITIVE = 0
		if (selected_overplot le 0) then return
	end

	field_x = reform ((oversets[selected_snapshot].(selected_overplot-1)[*,*,*,0])[cut], num_x, num_y, num_z)
	field_y = reform ((oversets[selected_snapshot].(selected_overplot-1)[*,*,*,1])[cut], num_x, num_y, num_z)
	field_z = reform ((oversets[selected_snapshot].(selected_overplot-1)[*,*,*,2])[cut], num_x, num_y, num_z)
	if (destretch[0]) then begin
		field_x = pc_destretch (field_x, coord.x, dim=1)
		field_y = pc_destretch (field_y, coord.x, dim=1)
		field_z = pc_destretch (field_z, coord.x, dim=1)
	end
	if (destretch[1]) then begin
		field_x = pc_destretch (field_x, coord.y, dim=2)
		field_y = pc_destretch (field_y, coord.y, dim=2)
		field_z = pc_destretch (field_z, coord.y, dim=2)
	end
	if (destretch[2]) then begin
		field_x = pc_destretch (field_x, coord.z, dim=3)
		field_y = pc_destretch (field_y, coord.z, dim=3)
		field_z = pc_destretch (field_z, coord.z, dim=3)
	end

	; find minimum and maximum values and set slider position
	if (not overplot_contour) then over_max = max ([field_x, field_y, field_z]) > 1.e-42

	if (overplot_contour eq 1) then begin
		; setup contour plot
		field_x_y = field_x
		field_x_z = 0.0
		field_y_x = field_y
		field_y_z = 0.0
		field_z_x = field_z
		field_z_y = 0.0

		; setup field indices
		field_x_indices = (findgen (num_x) + 0.25) / num_x
		field_y_indices = (findgen (num_y) + 0.25) / num_y
		field_z_indices = (findgen (num_z) + 0.25) / num_z
	end else begin
		; setup vector field
		field_x_y = congrid (field_x, num_x*bin_x/vector_distance, num_y, num_z*bin_z/vector_distance, /center)
		field_x_z = congrid (field_x, num_x*bin_x/vector_distance, num_y*bin_y/vector_distance, num_z, /center)
		field_y_x = congrid (field_y, num_x, num_y*bin_y/vector_distance, num_z*bin_z/vector_distance, /center)
		field_y_z = congrid (field_y, num_x*bin_x/vector_distance, num_y*bin_y/vector_distance, num_z, /center)
		field_z_x = congrid (field_z, num_x, num_y*bin_y/vector_distance, num_z*bin_z/vector_distance, /center)
		field_z_y = congrid (field_z, num_x*bin_x/vector_distance, num_y, num_z*bin_z/vector_distance, /center)

		; setup field indices
		field_x_indices = (findgen (num_x*bin_x/vector_distance) + 0.5) / (num_x*bin_x/vector_distance)
		field_y_indices = (findgen (num_y*bin_y/vector_distance) + 0.5) / (num_y*bin_y/vector_distance)
		field_z_indices = (findgen (num_z*bin_z/vector_distance) + 0.5) / (num_z*bin_z/vector_distance)

		WIDGET_CONTROL, sl_over, SENSITIVE = (selected_overplot gt 0)
		WIDGET_CONTROL, sl_over, SET_VALUE = [ pos_over[selected_overplot], 0.0, 2.0 ]
	end

	return
end


; Resets everything and redraws the window
pro cslice_reset_GUI

	common varset_common, set, overplot, oversets, unit, coord, varsets, varfiles, datadir, sources, param, run_param, var_list
	common cslice_common, cube, field, num_cubes, num_overs, num_snapshots
	common slider_common, bin_x, bin_y, bin_z, num_x, num_y, num_z, pos_b, pos_t, pos_over, val_min, val_max, val_range, over_max, dimensionality, frozen
	common gui_common, wimg_yz, wimg_xz, wimg_xy, wcut_x, wcut_y, wcut_z, co_t, sl_t, co_x, co_y, co_z, sl_x, sl_y, sl_z, b_load, b_log, b_sub, b_abs, b_cro, b_des, aver, timeser, vars, over, snap, prev, next, play, image, sl_min, sl_max, sl_over, min_max, freeze, jump_min, jump_max
	common settings_common, px, py, pz, cut, log_plot, abs_scale, show_cross, show_cuts, sub_aver, selected_cube, selected_overplot, selected_snapshot, af_x, af_y, af_z, destretch

	selected_cube = 0
	selected_overplot = 0
	selected_snapshot = 0
	log_plot = 0
	sub_aver = 0
	abs_scale = 1
	show_cross = 1
	destretch = [0, 0, 0]
	frozen = 0
	px = num_x / 2
	py = num_y / 2
	pz = num_z / 2
	pos_b = replicate (!VALUES.D_NAN, num_cubes, 2, 2)
	pos_t = replicate (!VALUES.D_NAN, num_cubes, 2, 2)
	pos_over = replicate (1.0, num_overs)
	over_max = 0.0

	cslice_prepare_cube, -1
	cslice_prepare_overplot

	WIDGET_CONTROL, b_log, SET_VALUE = log_plot
	WIDGET_CONTROL, b_abs, SET_VALUE = abs_scale
	WIDGET_CONTROL, b_sub, SET_VALUE = sub_aver
	WIDGET_CONTROL, b_cro, SET_VALUE = show_cross
	if (any (coord.lequidist eq 0)) then WIDGET_CONTROL, b_des, SET_VALUE = any (destretch)
	WIDGET_CONTROL, co_t, SET_VALUE = varfiles[selected_snapshot].time
	WIDGET_CONTROL, sl_t, SET_VALUE = num_snapshots-1-selected_snapshot
	WIDGET_CONTROL, co_x, SET_VALUE = coord.x[px]
	WIDGET_CONTROL, co_y, SET_VALUE = coord.y[py]
	WIDGET_CONTROL, co_z, SET_VALUE = coord.z[pz]
	WIDGET_CONTROL, sl_x, SET_VALUE = px
	WIDGET_CONTROL, sl_y, SET_VALUE = py
	WIDGET_CONTROL, sl_z, SET_VALUE = pz
	WIDGET_CONTROL, sl_min, SET_VALUE = pos_b[selected_cube,sub_aver]
	WIDGET_CONTROL, sl_max, SET_VALUE = pos_t[selected_cube,sub_aver]
	WIDGET_CONTROL, sl_over, SET_VALUE = [1.0, 0.0, 2.0]
	WIDGET_CONTROL, vars, SET_DROPLIST_SELECT = selected_cube
	WIDGET_CONTROL, over, SET_DROPLIST_SELECT = selected_overplot
	WIDGET_CONTROL, snap, SET_DROPLIST_SELECT = selected_snapshot

	WIDGET_CONTROL, prev, SENSITIVE = (num_snapshots ge 2)
	WIDGET_CONTROL, next, SENSITIVE = 0
	WIDGET_CONTROL, min_max, SENSITIVE = 1
	WIDGET_CONTROL, freeze, SENSITIVE = 1
	WIDGET_CONTROL, sl_over, SENSITIVE = 0
end


; Sophisticated interface with caching of VAR-files
pro cmp_cslice_cache, set_names, set_content=set_content, set_files=set_files, limits=limits, units=units, coords=coords, scaling=scaling, overplots=overplots

	common varset_common, set, overplot, oversets, unit, coord, varsets, varfiles, datadir, sources, param, run_param, var_list
	common cslice_common, cube, field, num_cubes, num_overs, num_snapshots
	common event_common, button_pressed_yz, button_pressed_xz, button_pressed_xy
	common slider_common, bin_x, bin_y, bin_z, num_x, num_y, num_z, pos_b, pos_t, pos_over, val_min, val_max, val_range, over_max, dimensionality, frozen
	common gui_common, wimg_yz, wimg_xz, wimg_xy, wcut_x, wcut_y, wcut_z, co_t, sl_t, co_x, co_y, co_z, sl_x, sl_y, sl_z, b_load, b_log, b_sub, b_abs, b_cro, b_des, aver, timeser, vars, over, snap, prev, next, play, image, sl_min, sl_max, sl_over, min_max, freeze, jump_min, jump_max
	common settings_common, px, py, pz, cut, log_plot, abs_scale, show_cross, show_cuts, sub_aver, selected_cube, selected_overplot, selected_snapshot, af_x, af_y, af_z, destretch

	; DEFAULT SETTINGS:
	show_cuts = 1
	log_plot = 0
	sub_aver = 0
	abs_scale = 1
	show_cross = 1
	destretch = [0, 0, 0]
	frozen = 0
	selected_cube = 0
	selected_overplot = 0
	selected_snapshot = 0
	af_x = 0
	af_y = 0
	af_z = 0
	min_size = 8


	set = set_names
	if (keyword_set (set_content)) then varsets = set_content
	if (keyword_set (set_files)) then varfiles = set_files
	if (not keyword_set (overplots)) then begin
		overplots = {none:'none'}
	end else begin
		if (not any (strlowcase (tag_names (overplots)) eq 'none')) then begin
			overplots = create_struct ({none:'none'}, overplots)
		end
	end
	overplot = overplots

	if (keyword_set (units)) then unit = units
	if (not keyword_set (unit) ) then begin
		print, "WARNING: setting units to unity."
		unit = { length:1, default_length:1, default_length_str:'-', velocity:1, default_velocity:1, default_velocity_str:'-', time:1, default_time:1, default_time_str:'-', temperature:1, default_temperature:1, default_temperature_str:'-', density:1, default_density:1, default_density_str:'-', mass:1, default_mass:1, default_mass_str:'-', magnetic_field:1, default_magnetic_field:1, default_magnetic_field_str:'-', current_density:1, default_current_density:1, default_current_density_str:'-' }
	end

	if (not keyword_set (scaling)) then scaling = 1
	if (n_elements (scaling) eq 1) then scaling = [ scaling, scaling, scaling ]
	if (keyword_set (coords)) then coord = coords

	; setup dimensions
	dims = size (varsets[0].(0))
	dimensionality = dims[0]
	num_x = dims[1]
	if (dimensionality ge 2) then num_y = dims[2] else num_y = 1
	if (dimensionality ge 3) then num_z = dims[3] else num_z = 1

	; setup limits, if necessary
	if (not keyword_set (limits)) then begin
		nghost_x = (dims[1] - num_x) / 2
		if (dimensionality ge 2) then nghost_y = (dims[2] - num_y) / 2 else nghost_y = 0
		if (dimensionality ge 3) then nghost_z = (dims[3] - num_z) / 2 else nghost_z = 0
		limits = reform (nghost_x+spread(indgen(num_x),[1,2],[num_y,num_z]) + dims[1]*(nghost_y+spread(indgen(num_y),[0,2],[num_x,num_z])) + dims[1]*dims[2]*(nghost_z+spread(indgen(num_z),[0,1],[num_x,num_y])), num_x, num_y, num_z)
	end

	cut = reform (limits, num_x, num_y, num_z)

	; setup coordinate system, if necessary
	if (n_elements (coord) eq 0) then begin
		print, "WARNING: setting the pixel size to default length."
		nghost_x = (dims[1] - num_x) / 2
		if (dimensionality ge 2) then nghost_y = (dims[2] - num_y) / 2 else nghost_y = 0
		if (dimensionality ge 3) then nghost_z = (dims[3] - num_z) / 2 else nghost_z = 0
		lequidist = indgen (dimensionality) ne -1
		lperi = indgen (dimensionality) eq -1
		ldegenerated = indgen (dimensionality) eq -1
		coord = { x:findgen(num_x)*unit.length/unit.default_length, y:findgen(num_y)*unit.length/unit.default_length, z:findgen(num_z)*unit.length/unit.default_length, nx:num_x, ny:num_y, nz:num_z, l1:nghost_x, l2:nghost_x+num_x-1, m1:nghost_y, m2:nghost_y+num_y-1, n1:nghost_z, n2:nghost_z+num_z-1, lequidist:lequidist, lperi:lperi, ldegenerated:ldegenerated }
	end

	sl_min = 0
	sl_max = 0
	sl_over = 0
	min_max = 0
	freeze = 0

	num_snapshots = n_elements (varfiles)
	snaps = varfiles[*].title + " (" + strtrim (varfiles[*].time, 2) + " " + unit.default_time_str + ")"

	num_cubes = n_tags (set)
	tags = strarr (num_cubes)
	for pos=0, num_cubes-1 do tags[pos] = set.(pos)

	num_overs = n_tags (overplot)
	overs = strarr (num_overs)
	for pos=0, num_overs-1 do overs[pos] = overplot.(pos)

	pos_b = replicate (!VALUES.D_NAN, num_cubes, 2, 2)
	pos_t = replicate (!VALUES.D_NAN, num_cubes, 2, 2)
	pos_over = replicate (1.0, num_overs)

	b_des = !VALUES.D_NAN
	co_t = !VALUES.D_NAN
	sl_t = !VALUES.D_NAN


	cslice_prepare_set, 0
	cslice_prepare_cube, -1


	if (num_x*scaling[0] lt min_size) then scaling[0] = min_size / num_x
	if (num_y*scaling[1] lt min_size) then scaling[1] = min_size / num_y
	if (num_z*scaling[2] lt min_size) then scaling[2] = min_size / num_z

	bin_x = scaling[0]
	bin_y = scaling[1]
	bin_z = scaling[2]


	px = num_x / 2
	py = num_y / 2
	pz = num_z / 2


	vars_active = (num_cubes ge 2)
	over_active = (num_overs ge 2)
	snap_active = (num_snapshots ge 2)
	prev_active = (num_snapshots ge 2)
	next_active = 0
	if (file_test ('GUI_settings.xdr', /read)) then load_active = 1 else load_active = 0

	coord_x_active = (num_x gt 1)
	coord_y_active = (num_y gt 1)
	coord_z_active = (num_z gt 1)

	dimensionality = coord_x_active + coord_y_active + coord_z_active
	if (dimensionality eq 0) then begin
		print, "Are you sure, you want to visualize 0D data?"
		print, "(If yes, just type '.continue' and press the return key.)"
		stop
	end

	button_pressed_yz = 0
	button_pressed_xz = 0
	button_pressed_xy = 0

	sl_size	= ((2*num_x*bin_x+num_y*bin_y)/2.5 > (400+max([num_x*bin_x,num_y*bin_y,num_z*bin_z]))/2) < 500

	MOTHER	= WIDGET_BASE (title='compare cube-slices')
	BASE	= WIDGET_BASE (MOTHER, /col)
	CTRL	= WIDGET_BASE (BASE, /row)

	bcol	= WIDGET_BASE (CTRL, /col)
	tmp	= WIDGET_BUTTON (bcol, value='RESET', uvalue='RESET', xsize=100)
	b_load	= WIDGET_BUTTON (bcol, value='LOAD SETTINGS', uvalue='LOAD', xsize=100, sensitive=load_active)
	tmp	= WIDGET_BUTTON (bcol, value='SAVE SETTINGS', uvalue='SAVE', xsize=100)
	play	= WIDGET_BUTTON (bcol, value='PLAY', uvalue='PLAY', xsize=100, sensitive=snap_active)
	if (snap_active) then save_str='SAVE MOVIE' else save_str='SAVE IMAGE'
	image	= WIDGET_BUTTON (bcol, value=save_str, uvalue='IMAGE', xsize=100)
	tmp	= WIDGET_BUTTON (bcol, value='QUIT', uvalue='QUIT', xsize=100)

	bcol	= WIDGET_BASE (CTRL, /col)
	bcot	= WIDGET_BASE (bcol, /row)
	snap	= WIDGET_DROPLIST (bcot, value=snaps, uvalue='SNAP', sensitive=snap_active, EVENT_PRO=cslice_event, title='snapshot:')
	prev	= WIDGET_BUTTON (bcot, value='-', uvalue='PREV', sensitive=prev_active, EVENT_PRO=cslice_event)
	next	= WIDGET_BUTTON (bcot, value='+', uvalue='NEXT', sensitive=next_active, EVENT_PRO=cslice_event)
	bcot	= WIDGET_BASE (bcol, /row)
	vars	= WIDGET_DROPLIST (bcot, value=tags, uvalue='VAR', sensitive=vars_active, EVENT_PRO=cslice_event, title='quantity:')
	bcot	= WIDGET_BASE (bcol, /col, frame=1)
	over	= WIDGET_DROPLIST (bcot, value=overs, uvalue='OVER', sensitive=over_active, EVENT_PRO=cslice_event, title='overplot:')
	bsubrow	= WIDGET_BASE (bcot, /row)
	tmp	= WIDGET_LABEL (bsubrow, value='vector length:', frame=0)
	sl_over	= CW_FSLIDER (bsubrow, uvalue='SCALE_OVER', /double, /edit, /suppress_value, min=0.0, max=2.0, /drag, value=1.0, xsize=150)
	WIDGET_CONTROL, sl_over, SENSITIVE = (selected_overplot gt 0)

	bcol	= WIDGET_BASE (CTRL, /col)
	timeser	= WIDGET_BUTTON (bcol, value='timeseries analysis', uvalue='SHOW_TIME')
	aver	= WIDGET_BUTTON (bcol, value='vertical profile', uvalue='SHOW_AVER')
	b_log	= CW_BGROUP (bcol, 'logarithmic plot', /nonexcl, uvalue='LOG_PLOT', set_value=log_plot)
	b_sub	= CW_BGROUP (bcol, 'substract profile', /nonexcl, uvalue='SUB_AVER', set_value=sub_aver)
	b_abs	= CW_BGROUP (bcol, 'absolute scaling', /nonexcl, uvalue='ABS_SCALE', set_value=abs_scale)
	b_cro	= CW_BGROUP (bcol, 'show crosshairs', /nonexcl, uvalue='SHOW_CROSS', set_value=show_cross)
	if (any (coord.lequidist eq 0)) then $
		b_des = CW_BGROUP (bcol, 'destretch grid', /nonexcl, uvalue='DESTRETCH', set_value=0)

	scol	= WIDGET_BASE (CTRL, /col)

	if (unit.default_time_str) then begin
		title_add = ' ['+unit.default_time_str+']:'
	end else begin
		title_add = ':'
	end
	scot	= WIDGET_BASE (scol, /row, /base_align_center)
	co_t	= CW_FIELD (scot, title='time'+title_add, uvalue='COT', value=varfiles[selected_snapshot].time, noedit=1-snap_active, /floating, /return_events, xsize=10)
	sl_t	= WIDGET_SLIDER (scot, uvalue='SLT', /suppress_value, min=0, max=(num_snapshots-1)>1, /drag, value=num_snapshots-1-selected_snapshot, xsize=(((num_snapshots+1)*8>64)+10)<256, sensitive=snap_active)
	if (unit.default_length_str) then begin
		title_add = ' ['+unit.default_length_str+']:'
		co_int = 0
	end else begin
		title_add = ':'
		co_int = 1
	end
	scot	= WIDGET_BASE (scol, /row, /base_align_center)
	co_x	= CW_FIELD (scot, title='X'+title_add, uvalue='COX', value=coord.x[px], noedit=1-coord_x_active, integer=co_int, floating=(1-co_int), /return_events, xsize=12)
	sl_x	= WIDGET_SLIDER (scot, uvalue='SLX', value=px, min=0, max=(num_x-1)>1, xsize=((num_x*bin_x>128)+10)<512, /drag, sensitive=coord_x_active)
	scot	= WIDGET_BASE (scol, /row, /base_align_center)
	co_y	= CW_FIELD (scot, title='Y'+title_add, uvalue='COY', value=coord.y[py], noedit=1-coord_y_active, integer=co_int, floating=(1-co_int), /return_events, xsize=12)
	sl_y	= WIDGET_SLIDER (scot, uvalue='SLY', value=py, min=0, max=(num_y-1)>1, xsize=((num_y*bin_y>128)+10)<512, /drag, sensitive=coord_y_active)
	scot	= WIDGET_BASE (scol, /row, /base_align_center)
	co_z	= CW_FIELD (scot, title='Z'+title_add, uvalue='COZ', value=coord.z[pz], noedit=1-coord_z_active, integer=co_int, floating=(1-co_int), /return_events, xsize=12)
	sl_z	= WIDGET_SLIDER (scot, uvalue='SLZ', value=pz, min=0, max=(num_z-1)>1, xsize=((num_z*bin_z>128)+10)<512, /drag, sensitive=coord_z_active)

	DISP	= WIDGET_BASE (BASE, /row)
	MID	= WIDGET_BASE (BASE, /col)
	bcot	= WIDGET_BASE (MID, /row)

	sl_min	= CW_FSLIDER (bcot, title='lower value (black level)', uvalue='SCALE_BOT', /double, /edit, min=val_range[0], max=val_range[1], /drag, value=val_min, xsize=sl_size)
	bcol	= WIDGET_BASE (bcot, /col)
	min_max	= WIDGET_BUTTON (bcol, value='<= min SET max =>', uvalue='MIN_MAX', xsize=120)
	freeze	= WIDGET_BUTTON (bcol, value='FREEZE RANGE', uvalue='FREEZE', xsize=120)
	sl_max	= CW_FSLIDER (bcot, title='upper value (white level)', uvalue='SCALE_TOP', /double, /edit, min=val_range[0], max=val_range[1], /drag, value=val_max, xsize=sl_size)

	bsubcol	= WIDGET_BASE (bcot, /col)
	tmp	= WIDGET_LABEL (bsubcol, value='Jump to:', frame=0)
	bsubrow	= WIDGET_BASE (bsubcol, /row)
	jump_min= WIDGET_BUTTON (bsubrow, value='MIN', uvalue='JUMP_MIN')
	jump_max= WIDGET_BUTTON (bsubrow, value='MAX', uvalue='JUMP_MAX')

	cut_height = min ([num_x*bin_x,num_y*bin_y,num_z*bin_z]) > 256
	CUTS	= WIDGET_BASE (BASE, /row)

	dimg_yz	= WIDGET_DRAW (DISP, UVALUE='DRAW_YZ', xsize=num_y*bin_y, ysize=num_z*bin_z, retain=2, /button_events, /motion_events)
	dimg_xz	= WIDGET_DRAW (DISP, UVALUE='DRAW_XZ', xsize=num_x*bin_x, ysize=num_z*bin_z, retain=2, /button_events, /motion_events)
	dimg_xy	= WIDGET_DRAW (DISP, UVALUE='DRAW_XY', xsize=num_x*bin_x, ysize=num_y*bin_y, retain=2, /button_events, /motion_events)

	dcut_y	= WIDGET_DRAW (CUTS, xsize=num_y*bin_y, ysize=cut_height, retain=2)
	dcut_x	= WIDGET_DRAW (CUTS, xsize=num_x*bin_x, ysize=cut_height, retain=2)
	dcut_z	= WIDGET_DRAW (CUTS, xsize=num_z*bin_z>128, ysize=cut_height, retain=2)

	WIDGET_CONTROL, BASE

	WIDGET_CONTROL, MOTHER, /REALIZE
	WIDGET_CONTROL, dimg_yz, GET_VALUE = wimg_yz
	WIDGET_CONTROL, dimg_xz, GET_VALUE = wimg_xz
	WIDGET_CONTROL, dimg_xy, GET_VALUE = wimg_xy
	WIDGET_CONTROL, dcut_y, GET_VALUE = wcut_y
	WIDGET_CONTROL, dcut_x, GET_VALUE = wcut_x
	WIDGET_CONTROL, dcut_z, GET_VALUE = wcut_z

	XMANAGER, "cslice", MOTHER, /no_block

	cslice_reset_GUI
	cslice_draw, 1, 1, 1
end

