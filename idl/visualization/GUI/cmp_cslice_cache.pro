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
function get_range, data

	common slider_common, bin_x, bin_y, bin_z, num_x, num_y, num_z, pos_b, pos_t, csmin, csmax, dimensionality
	common settings_common, px, py, pz, cut, abs_scale, show_cross, show_cuts, sub_aver, selected_cube, selected_overplot, selected_snapshot, af_x, af_y, af_z

	if (abs_scale) then begin
		min = csmin
		max = csmax
	end else begin
		min = min (data)
		max = max (data)
	end

	if (min eq max) then begin
		; extend value range a little, if necessary (must have min < max)
		; a uniform value should appear as a 50% saturation gray
		if (min eq 0.0) then begin
			min = -1d-42
			max = 1d-42
		endif else begin
			min *= 0.99999
			max *= 1.00001
		endelse
	endif

	return, [min, max]
end


; Event handling of visualisation window
pro cslice_event, event

	common event_common, button_pressed_yz, button_pressed_xz, button_pressed_xy
	common cslice_common, cube, field, num_cubes, num_overs, num_snapshots
	common slider_common, bin_x, bin_y, bin_z, num_x, num_y, num_z, pos_b, pos_t, csmin, csmax, dimensionality
	common gui_common, wimg_yz, wimg_xz, wimg_xy, wcut_x, wcut_y, wcut_z, sl_x, sl_y, sl_z, b_abs, b_sub, b_cro, aver, vars, over, snap, prev, next, play, scal_b, scal_t
	common settings_common, px, py, pz, cut, abs_scale, show_cross, show_cuts, sub_aver, selected_cube, selected_overplot, selected_snapshot, af_x, af_y, af_z

	WIDGET_CONTROL, WIDGET_INFO (event.top, /CHILD)

	; SETTINGS:
	; filename for saving the settings
	settings_file = 'GUI_settings.xdr'

	; time in seconds to wait after showing each frame for 1D, 2D, and 3D movies (0=fastest possible)
	min_wait_time = [ 0.2, 0.1, 0.05 ]

	quit = -1
	DRAW_IMAGE_1=0  &  DRAW_IMAGE_2=0  &  DRAW_IMAGE_3=0

	WIDGET_CONTROL, event.id, GET_UVALUE = eventval


	CASE eventval of
	'SCL':  begin 
		abs_scale = event.select
		DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
	end
	'SUB_AVER':  begin 
		pos_b[selected_cube,sub_aver] = csmin
		pos_t[selected_cube,sub_aver] = csmax
		sub_aver = event.select
		prepare_cube, -1
		DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
	end
	'CHS':  begin 
		show_cross = event.select
		DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
	end
	'SHOW_AVER':  begin 
		draw_averages, selected_snapshot
	end
	'SLX':  begin
		WIDGET_CONTROL, event.id, GET_VALUE = pos
		px = pos
		DRAW_IMAGE_1 = 1
	end
	'SLY':  begin
		WIDGET_CONTROL, event.id, GET_VALUE = pos
		py = pos
		DRAW_IMAGE_2 = 1
	end
	'SLZ':  begin
		WIDGET_CONTROL, event.id, GET_VALUE = pos
		pz = pos
		DRAW_IMAGE_3 = 1
	end
	'DRAW_YZ':  begin
		if (event.press) then button_pressed_yz = 1
		if (button_pressed_yz) then begin
			last_py = py
			last_pz = pz
			py = event.x / bin_y > 0 < (num_y-1)
			pz = event.y / bin_z > 0 < (num_z-1)
			if ((py ne last_py) or (pz ne last_pz)) then begin
				WIDGET_CONTROL, sl_y, SET_VALUE = py
				WIDGET_CONTROL, sl_z, SET_VALUE = pz
				DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
			end
		endif
		if (event.release) then button_pressed_yz = 0
	end
	'DRAW_XZ':  begin
		if (event.press) then button_pressed_xz = 1
		if (button_pressed_xz) then begin
			last_px = px
			last_pz = pz
			px = event.x / bin_x > 0 < (num_x-1)
			pz = event.y / bin_z > 0 < (num_z-1)
			if ((px ne last_px) or (pz ne last_pz)) then begin
				WIDGET_CONTROL, sl_x, SET_VALUE = px
				WIDGET_CONTROL, sl_z, SET_VALUE = pz
				DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
			end
		endif
		if (event.release) then button_pressed_xz = 0
	end
	'DRAW_XY':  begin
		if (event.press) then button_pressed_xy = 1
		if (button_pressed_xy) then begin
			last_px = px
			last_py = py
			px = event.x / bin_x > 0 < (num_x-1)
			py = event.y / bin_y > 0 < (num_y-1)
			if ((px ne last_px) or (py ne last_py)) then begin
				WIDGET_CONTROL, sl_x, SET_VALUE = px
				WIDGET_CONTROL, sl_y, SET_VALUE = py
				DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
			end
		endif
		if (event.release) then button_pressed_xy = 0
	end
	'IMG_SCLB': begin
		WIDGET_CONTROL, scal_b, GET_VALUE = csmin
		if (csmin gt csmax) then begin
			csmin = csmax
			WIDGET_CONTROL, scal_b, SET_VALUE = csmin
		end
		DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
	end
	'IMG_SCLT': begin
		WIDGET_CONTROL, scal_t, GET_VALUE = csmax
		if (csmax lt csmin) then begin
			csmax = csmin
			WIDGET_CONTROL, scal_t, SET_VALUE = csmax
		end
		DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
	end
	'VAR': begin
		last = selected_cube
		selected_cube = event.index
		if (last ne event.index) then begin
			prepare_cube, last
			DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
		end
	end
	'SNAP': begin
		last = selected_cube
		if (selected_snapshot ne event.index) then begin
			prepare_set, event.index
			prepare_cube, last
			DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
			if (selected_snapshot lt num_snapshots - 1) then prev_active = 1 else prev_active = 0
			if (selected_snapshot gt 0) then next_active = 1 else next_active = 0
			WIDGET_CONTROL, prev, SENSITIVE = prev_active
			WIDGET_CONTROL, next, SENSITIVE = next_active

			window, 0, xsize=8, ysize=8, retain=2
			!P.MULTI = [0, 1, 1]
			wdelete
		end
	end
	'NEXT': begin
		selected_snapshot -= 1
		WIDGET_CONTROL, snap, SET_DROPLIST_SELECT = selected_snapshot
		prepare_set, selected_snapshot
		prepare_cube, selected_cube
		WIDGET_CONTROL, prev, SENSITIVE = 1
		if (selected_snapshot le 0) then WIDGET_CONTROL, next, SENSITIVE = 0
		DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
	end
	'PREV': begin
		selected_snapshot += 1
		WIDGET_CONTROL, snap, SET_DROPLIST_SELECT = selected_snapshot
		prepare_set, selected_snapshot
		prepare_cube, selected_cube
		if (selected_snapshot ge num_snapshots - 1) then WIDGET_CONTROL, prev, SENSITIVE = 0
		WIDGET_CONTROL, next, SENSITIVE = 1
		DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
	end
	'OVER': begin
		if (selected_overplot ne event.index) then begin
			selected_overplot = event.index
			prepare_overplot
			DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
		end
	end
	'RESET': begin
		reset_GUI
		DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
	end
	'LOAD': begin
		if (file_test (settings_file, /read)) then begin
			restore, settings_file

			if (px gt (num_x - 1)) then px = num_x - 1
			if (py gt (num_y - 1)) then py = num_y - 1
			if (pz gt (num_z - 1)) then pz = num_z - 1

			prepare_cube, -1
			prepare_overplot
			DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1

			WIDGET_CONTROL, b_abs, SET_VALUE = abs_scale
			WIDGET_CONTROL, b_sub, SET_VALUE = sub_aver
			WIDGET_CONTROL, b_cro, SET_VALUE = show_cross
			WIDGET_CONTROL, sl_x, SET_VALUE = px
			WIDGET_CONTROL, sl_y, SET_VALUE = py
			WIDGET_CONTROL, sl_z, SET_VALUE = pz
			WIDGET_CONTROL, scal_b, SET_VALUE = pos_b[selected_cube,sub_aver]
			WIDGET_CONTROL, scal_t, SET_VALUE = pos_t[selected_cube,sub_aver]
			WIDGET_CONTROL, vars, SET_DROPLIST_SELECT = selected_cube
			WIDGET_CONTROL, over, SET_DROPLIST_SELECT = selected_overplot
			WIDGET_CONTROL, snap, SET_DROPLIST_SELECT = selected_snapshot
		end
	end
	'SAVE': begin
		pos_b[selected_cube,sub_aver] = csmin
		pos_t[selected_cube,sub_aver] = csmax
		save, filename=settings_file, num_cubes, selected_cube, selected_overplot, abs_scale, sub_aver, show_cross, px, py, pz, pos_b, pos_t
	end
	'PLAY': begin
		WIDGET_CONTROL, vars, SENSITIVE = 0
		WIDGET_CONTROL, over, SENSITIVE = 0
		WIDGET_CONTROL, snap, SENSITIVE = 0
		WIDGET_CONTROL, prev, SENSITIVE = 0
		WIDGET_CONTROL, next, SENSITIVE = 0
		WIDGET_CONTROL, play, SENSITIVE = 0
		WIDGET_CONTROL, aver, SENSITIVE = 0
		previous_snapshot = selected_snapshot
		if (num_snapshots gt 1) then begin
			for i = num_snapshots-1, 1, -1 do begin
				prepare_set, i
				prepare_cube, -1, 0
				draw_images, 1, 1, 1
				wait, min_wait_time[dimensionality-1]
			end
		end
		prepare_set, previous_snapshot
		prepare_cube, -1, 0
		DRAW_IMAGE_1=1  &  DRAW_IMAGE_2=1  &  DRAW_IMAGE_3=1
		if (num_cubes ge 2) then vars_active = 1 else vars_active = 0
		if (num_overs ge 2) then over_active = 1 else over_active = 0
		if (num_snapshots ge 2) then snap_active = 1 else snap_active = 0
		if (selected_snapshot lt num_snapshots - 1) then prev_active = 1 else prev_active = 0
		if (selected_snapshot gt 0) then next_active = 1 else next_active = 0
		WIDGET_CONTROL, vars, SENSITIVE = vars_active
		WIDGET_CONTROL, over, SENSITIVE = over_active
		WIDGET_CONTROL, snap, SENSITIVE = snap_active
		WIDGET_CONTROL, prev, SENSITIVE = prev_active
		WIDGET_CONTROL, next, SENSITIVE = next_active
		WIDGET_CONTROL, play, SENSITIVE = 1
		WIDGET_CONTROL, aver, SENSITIVE = 1

		if (show_cuts) then begin
			window, 0, xsize=8, ysize=8, retain=2
			!P.MULTI = [0, 1, 1]
			wdelete
		end
	end
	'QUIT': begin
		quit = event.top
	end
	endcase

	draw_images, DRAW_IMAGE_1, DRAW_IMAGE_2, DRAW_IMAGE_3

	WIDGET_CONTROL, WIDGET_INFO (event.top, /CHILD)

	if (quit ge 0) then  WIDGET_CONTROL, quit, /DESTROY

	return
end


; Draws the slices into the window
pro draw_images, DRAW_IMAGE_1, DRAW_IMAGE_2, DRAW_IMAGE_3

	common cslice_common, cube, field, num_cubes, num_overs, num_snapshots
	common overplot_common, overplot_contour, field_x_y, field_x_z, field_y_x, field_y_z, field_z_x, field_z_y, field_x_indices, field_y_indices, field_z_indices, vector_distance, vector_length, field_x_max, field_y_max, field_z_max
	common slider_common, bin_x, bin_y, bin_z, num_x, num_y, num_z, pos_b, pos_t, csmin, csmax, dimensionality
	common gui_common, wimg_yz, wimg_xz, wimg_xy, wcut_x, wcut_y, wcut_z, sl_x, sl_y, sl_z, b_abs, b_sub, b_cro, aver, vars, over, snap, prev, next, play, scal_b, scal_t
	common settings_common, px, py, pz, cut, abs_scale, show_cross, show_cuts, sub_aver, selected_cube, selected_overplot, selected_snapshot, af_x, af_y, af_z

	; stepping of crosshairs
	step = 4

	; number of levels for contour plot
	default, nlevels, 50

	!P.MULTI = [0, 1, 1]

	ox = round (bin_x / 2.0) - 1
	oy = round (bin_y / 2.0) - 1
	oz = round (bin_z / 2.0) - 1

	num_over_x = n_elements (field_x_indices)
	num_over_y = n_elements (field_y_indices)
	num_over_z = n_elements (field_z_indices)

	if (DRAW_IMAGE_1 or DRAW_IMAGE_2 or DRAW_IMAGE_3) then begin
		ii = (reform (cube[px,*,*], num_y, num_z) > csmin) < csmax
		if ((bin_y ne 1) or (bin_z ne 1)) then ii = congrid (ii, fix (num_y*bin_y), fix (num_z*bin_z), cubic = 0)
		range = get_range (ii)
		wset, wimg_yz
		if (show_cross) then begin
			if (py gt af_y) then for i = fix ((py-af_y)*bin_y), 0, -step do ii[i:i+1, pz*bin_z+oz] = range
			if (py lt num_y-1-af_y) then for i = fix ((py+af_y)*bin_y), fix (num_y*bin_y)-2, step do ii[i:i+1, pz*bin_z+oz] = range
			if (pz gt af_z) then for i = fix ((pz-af_z)*bin_z), 0, -step do ii[py*bin_y+oy, i:i+1] = range
			if (pz lt num_z-1-af_z) then for i = fix ((pz+af_z)*bin_z), fix (num_z*bin_z)-2, step do ii[py*bin_y+oy, i:i+1] = range
			if ((py le af_y) and (py lt num_y-1-af_y) and (pz le af_z) and (pz lt num_z-1-af_z)) then ii[0:1, 0] = range
		end $
		else if (abs_scale) then ii[0:1, 0] = range
		tvscl, ii
		if (selected_overplot gt 0) then begin
			if (overplot_contour eq 1) then begin
				contour, reform (field_x_y[px, *, *], num_over_y, num_over_z), field_y_indices, field_z_indices, nlevels=nlevels, xs=4, ys=4, color=200, /noerase, pos=[0.0,0.0,1.0,1.0]
			end else begin
				velovect, reform (field_y_x[px, *, *], num_over_y, num_over_z), reform (field_z_x[px, *, *], num_over_y, num_over_z), field_y_indices, field_z_indices, length=vector_length, xr=[0.0,1.0], yr=[0.0,1.0], xs=4, ys=4, color=200, /noerase, pos=[0.0,0.0,1.0,1.0]
			end
		end
		if (show_cuts and (DRAW_IMAGE_1 or DRAW_IMAGE_3)) then begin
			wset, wcut_x
			plot, cube[px,*,pz], xrange=[0,num_y], yrange=range, xstyle=1, ystyle=1, xmargin=[0,0], ymargin=[0,0]
			axis, 0, 0, xaxis=1, xstyle=1, ystyle=1
			axis, 0, 0, yaxis=1, xstyle=1, ystyle=1
		end
	end

	if (DRAW_IMAGE_1 or DRAW_IMAGE_2 or DRAW_IMAGE_3) then begin
		ii = (reform (cube[*, py, *], num_x, num_z) > csmin) < csmax
		if ((bin_x ne 1) or (bin_z ne 1)) then ii = congrid (ii, fix (num_x*bin_x), fix (num_z*bin_z), cubic = 0)
		range = get_range (ii)
		wset, wimg_xz
		if (show_cross) then begin
			if (px gt af_x) then for i = fix ((px-af_x)*bin_x), 0, -step do ii[i:i+1, pz*bin_z+oz] = range
			if (px lt num_x-1-af_x) then for i = fix ((px+af_x)*bin_x), fix (num_x*bin_x)-2, step do ii[i:i+1, pz*bin_z+oz] = range
			if (pz gt af_z) then for i = fix ((pz-af_z)*bin_z), 0, -step do ii[px*bin_x+ox, i:i+1] = range
			if (pz lt num_z-1-af_z) then for i = fix ((pz+af_z)*bin_z), fix (num_z*bin_z)-2, step do ii[px*bin_x+ox, i:i+1] = range
			if ((px le af_x) and (px ge num_x-1-af_x) and (pz le af_z) and (pz ge num_z-1-af_z)) then ii[0:1, 0] = range
		end $
		else if (abs_scale) then ii[0:1, 0] = range
		tvscl, ii
		if (selected_overplot gt 0) then begin
			if (overplot_contour eq 1) then begin
				contour, reform (field_y_x[*, py, *], num_over_x, num_over_z), field_x_indices, field_z_indices, nlevels=nlevels, xs=4, ys=4, color=200, /noerase, pos=[0.0,0.0,1.0,1.0]
			end else begin
				velovect, reform (field_x_y[*, py, *], num_over_x, num_over_z), reform (field_z_y[*, py, *], num_over_x, num_over_z), field_x_indices, field_z_indices, length=vector_length, xr=[0.0,1.0], yr=[0.0,1.0], xs=4, ys=4, color=200, /noerase, pos=[0.0,0.0,1.0,1.0]
			end
		end
		if (show_cuts and (DRAW_IMAGE_2 or DRAW_IMAGE_3)) then begin
			wset, wcut_y
			plot, cube[*,py,pz], xrange=[0,num_x], yrange=range, xstyle=1, ystyle=1, xmargin=[0,0], ymargin=[0,0]
			axis, 0, 0, xaxis=1, xstyle=1, ystyle=1
			axis, 0, 0, yaxis=1, xstyle=1, ystyle=1
		end
	end

	if (DRAW_IMAGE_1 or DRAW_IMAGE_2 or DRAW_IMAGE_3) then begin
		ii = (reform (cube[*, *, pz], num_x, num_y) > csmin) < csmax
		if ((bin_x ne 1) or (bin_y ne 1)) then ii = congrid (ii, fix (num_x*bin_x), fix (num_y*bin_y), cubic = 0)
		range = get_range (ii)
		wset, wimg_xy
		if (show_cross) then begin
			if (px gt af_x) then for i = fix ((px-af_x)*bin_x), 0, -step do ii[i:i+1, py*bin_y+oy] = range
			if (px lt num_x-1-af_x) then for i = fix ((px+af_x)*bin_x), fix (num_x*bin_x)-2, step do ii[i:i+1, py*bin_y+oy] = range
			if (py gt af_y) then for i = fix ((py-af_y)*bin_y), 0, -step do ii[px*bin_x+ox, i:i+1] = range
			if (py lt num_y-1-af_y) then for i = fix ((py+af_y)*bin_y), fix (num_y*bin_y)-2, step do ii[px*bin_x+ox, i:i+1] = range
			if ((px le af_x) and (px ge num_x-1-af_x) and (py le af_y) and (py ge num_y-1-af_y)) then ii[0:1, 0] = range
		end $
		else if (abs_scale) then ii[0:1, 0] = range
		tvscl, ii
		if (selected_overplot gt 0) then begin
			if (overplot_contour eq 1) then begin
				contour, reform (field_z_x[*, *, pz], num_over_x, num_over_y), field_x_indices, field_y_indices, nlevels=nlevels, xs=4, ys=4, color=200, /noerase, pos=[0.0,0.0,1.0,1.0]
			end else begin
				velovect, reform (field_x_z[*, *, pz], num_over_x, num_over_y), reform (field_y_z[*, *, pz], num_over_x, num_over_y), field_x_indices, field_y_indices, length=vector_length, xr=[0.0,1.0], yr=[0.0,1.0], xs=4, ys=4, color=200, /noerase, pos=[0.0,0.0,1.0,1.0]
			end
		end
		if (show_cuts and (DRAW_IMAGE_1 or DRAW_IMAGE_2)) then begin
			wset, wcut_z
			plot, cube[px,py,*], xrange=[0,num_z], yrange=range, xstyle=1, ystyle=1, xmargin=[0,0], ymargin=[0,0]
			axis, 0, 0, xaxis=1, xstyle=1, ystyle=1
			axis, 0, 0, yaxis=1, xstyle=1, ystyle=1
		end
	end
end


; Draws horizontally averaged vertical profiles into a second window
pro draw_averages, number

	common varset_common, set, overplot, oversets, unit, coord, varsets, varfiles, sources
	common settings_common, px, py, pz, cut, abs_scale, show_cross, show_cuts, sub_aver, selected_cube, selected_overplot, selected_snapshot, af_x, af_y, af_z

	tags = tag_names (varsets[number])

	if (tags eq ['cube']) then begin
		window, 13, xsize=500, ysize=400, title = 'vertical profile analysis', retain=2
		!P.MULTI = [0, 1, 1]
		vert_prof, varsets[number].cube[cut], coord=coord.z, title = 'horizontal averages'
	end else begin
		window, 13, xsize=1000, ysize=800, title = 'vertical profile analysis', retain=2
		!P.MULTI = [0, 2, 2]
		max_subplots = 4
		num_subplots = 0
		if (any (strcmp (tags, 'ln_rho', /fold_case)) and (num_subplots lt max_subplots)) then begin
			num_subplots += 1
			vert_prof, exp (varsets[number].ln_rho[cut]), coord=coord.z, title = 'density ['+unit.default_density_str+']', log=1
		end else if (any (strcmp (tags, 'log_rho', /fold_case)) and (num_subplots lt max_subplots)) then begin
			num_subplots += 1
			vert_prof, 10.0^(varsets[number].log_rho[cut]), coord=coord.z, title = 'density ['+unit.default_density_str+']', log=1
		end else if (any (strcmp (tags, 'rho', /fold_case)) and (num_subplots lt max_subplots)) then begin
			num_subplots += 1
			vert_prof, varsets[number].rho[cut], coord=coord.z, title = 'density ['+unit.default_density_str+']'
		end
		if (any (strcmp (tags, 'u_abs', /fold_case)) and (num_subplots lt max_subplots)) then begin
			num_subplots += 1
			vert_prof, varsets[number].u_abs[cut], coord=coord.z, title = 'absolute velocity ['+unit.default_velocity_str+']'
		end
		if (any (strcmp (tags, 'rho_u_z', /fold_case)) and (num_subplots lt max_subplots)) then begin
			num_subplots += 1
			vert_prof, varsets[number].rho_u_z[cut], coord=coord.z, title = 'vertical impulse density ['+unit.default_density_str+' * '+unit.default_velocity_str+']'
		end
		if (any (strcmp (tags, 'Temp', /fold_case)) and (num_subplots lt max_subplots)) then begin
			num_subplots += 1
			vert_prof, varsets[number].Temp[cut], coord=coord.z, title = 'temperature [K]', log=1
		end
		if (any (strcmp (tags, 'j', /fold_case)) and (num_subplots lt max_subplots)) then begin
			num_subplots += 1
			vert_prof, varsets[number].j[cut], coord=coord.z, title = 'current density', log=1
		end
	end
end


; Prepares a data set for visualisation
pro prepare_set, i

	common varset_common, set, overplot, oversets, unit, coord, varsets, varfiles, sources
	common cslice_common, cube, field, num_cubes, num_overs, num_snapshots
	common settings_common, px, py, pz, cut, abs_scale, show_cross, show_cuts, sub_aver, selected_cube, selected_overplot, selected_snapshot, af_x, af_y, af_z

	selected_snapshot = i

	precalc, i

	num_cubes = n_tags (set)

	return
end


; Prepares a cube for visualisation
pro prepare_cube, last_index, update_slider

	common varset_common, set, overplot, oversets, unit, coord, varsets, varfiles, sources
	common cslice_common, cube, field, num_cubes, num_overs, num_snapshots
	common slider_common, bin_x, bin_y, bin_z, num_x, num_y, num_z, pos_b, pos_t, csmin, csmax, dimensionality
	common gui_common, wimg_yz, wimg_xz, wimg_xy, wcut_x, wcut_y, wcut_z, sl_x, sl_y, sl_z, b_abs, b_sub, b_cro, aver, vars, over, snap, prev, next, play, scal_b, scal_t
	common settings_common, px, py, pz, cut, abs_scale, show_cross, show_cuts, sub_aver, selected_cube, selected_overplot, selected_snapshot, af_x, af_y, af_z

	; SETTINGS:
	; fraction of box width to keep free of crosshairs at center
	af_fraction = 1.0 / 8.0
	; minimum size of crosshairs
	af_minimum = 6

	default, update_slider, 1

	; get selected cube from set
	tag = set.(selected_cube)
	res = execute ("cube = reform (varsets[selected_snapshot]."+tag+"[cut], num_x, num_y, num_z)")
	if (not res) then begin
		print, "Could not select dataset!"
		stop
	end

	; substract horizontal averages
	if (sub_aver) then for z=0, num_z-1 do cube[*,*,z] -= mean (cube [*,*,z])

	if (update_slider) then begin
		; find minimum and maximum values
		csmin = min (cube)
		csmax = max (cube)
		range = get_range (cube)

		if (last_index ge 0) then begin
			; save slider positions of previous cube
			WIDGET_CONTROL, scal_b, GET_VALUE = b
			pos_b(last_index, sub_aver) = b
			WIDGET_CONTROL, scal_t, GET_VALUE = t
			pos_t(last_index, sub_aver) = t
		end

		; set default slider positions (min/max)
		if (not finite (pos_b[selected_cube,sub_aver], /NaN)) then pos_b[selected_cube,sub_aver] = csmin
		if (not finite (pos_t[selected_cube,sub_aver], /NaN)) then pos_t[selected_cube,sub_aver] = csmax

		; adjust slider positions to fit inside value range
		if (pos_b[selected_cube,sub_aver] lt csmin) then pos_b[selected_cube,sub_aver] = csmin
		if (pos_t[selected_cube,sub_aver] gt csmax) then pos_t[selected_cube,sub_aver] = csmax

		; update slider
		if (scal_b ne 0) then WIDGET_CONTROL, scal_b, SET_VALUE = [ pos_b[selected_cube,sub_aver], range ]
		if (scal_t ne 0) then WIDGET_CONTROL, scal_t, SET_VALUE = [ pos_t[selected_cube,sub_aver], range ]

		; get desired min/max values from sliders
		csmin = pos_b[selected_cube,sub_aver]
		csmax = pos_t[selected_cube,sub_aver]
	end

	; determine dimesions
	num_x = (size (cube))[1]
	num_y = (size (cube))[2]
	num_z = (size (cube))[3]

	; setup crosshairs parameters
	af_x = round (num_x * af_fraction)
	af_y = round (num_y * af_fraction)
	af_z = round (num_z * af_fraction)
	if (af_x < af_minimum) then af_x = af_minimum
	if (af_y < af_minimum) then af_y = af_minimum
	if (af_z < af_minimum) then af_z = af_minimum

	return
end


; Prepares an overplot for visualisation
pro prepare_overplot

	common varset_common, set, overplot, oversets, unit, coord, varsets, varfiles, sources
	common cslice_common, cube, field, num_cubes, num_overs, num_snapshots
	common overplot_common, overplot_contour, field_x_y, field_x_z, field_y_x, field_y_z, field_z_x, field_z_y, field_x_indices, field_y_indices, field_z_indices, vector_distance, vector_length, field_x_max, field_y_max, field_z_max
	common slider_common, bin_x, bin_y, bin_z, num_x, num_y, num_z, pos_b, pos_t, csmin, csmax, dimensionality
	common gui_common, wimg_yz, wimg_xz, wimg_xy, wcut_x, wcut_y, wcut_z, sl_x, sl_y, sl_z, b_abs, b_sub, b_cro, aver, vars, over, snap, prev, next, play, scal_b, scal_t
	common settings_common, px, py, pz, cut, abs_scale, show_cross, show_cuts, sub_aver, selected_cube, selected_overplot, selected_snapshot, af_x, af_y, af_z

	; SETTINGS:
	; distance of vector footpoint locations
	vector_distance = 8
	; maximum length of vectors
	vector_length = vector_distance * 0.75
	; default plot routine: 0=velovect (1=contour)
	overplot_contour = 0

	if (selected_overplot le 0) then return

	; get selected overplot from set
	tag = overplot.(selected_overplot)
	res_x = execute ("field_x = oversets[selected_snapshot]."+tag+"[*,*,*,0]")
	res_y = execute ("field_y = oversets[selected_snapshot]."+tag+"[*,*,*,1]")
	res_z = execute ("field_z = oversets[selected_snapshot]."+tag+"[*,*,*,2]")
	if ((not res_x) or (not res_y) or (not res_z)) then begin
		print, "Could not select overplot dataset!"
		stop
	end

	size_field_x = size (field_x[cut])
	size_field_y = size (field_y[cut])
	size_field_z = size (field_z[cut])
	field_x = reform (field_x[cut], size_field_x[1], size_field_x[2], size_field_x[3])
	field_y = reform (field_y[cut], size_field_y[1], size_field_y[2], size_field_y[3])
	field_z = reform (field_z[cut], size_field_z[1], size_field_z[2], size_field_z[3])

	if (strpos (tag, "_velovect") gt 0) then overplot_contour = 0
	if (strpos (tag, "_contour") gt 0) then overplot_contour = 1

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

		; setup maximum values of x, y, and z overplots
		field_x_max = max (field_x)
		field_y_max = max (field_y)
		field_z_max = max (field_z)
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

		; setup vector lengthes for x, y, and z overplots
		field_x_max = max (field_x)
		field_y_max = max (field_y)
		field_z_max = max (field_z)

		; normalize maximum value to 1.0
;		field_x_y /= field_x_max
;		field_x_z /= field_x_max
;		field_y_x /= field_y_max
;		field_y_z /= field_y_max
;		field_z_x /= field_z_max
;		field_z_y /= field_z_max
	end

	return
end


; Resets everything and redraws the window
pro reset_GUI

	common cslice_common, cube, field, num_cubes, num_overs, num_snapshots
	common slider_common, bin_x, bin_y, bin_z, num_x, num_y, num_z, pos_b, pos_t, csmin, csmax, dimensionality
	common gui_common, wimg_yz, wimg_xz, wimg_xy, wcut_x, wcut_y, wcut_z, sl_x, sl_y, sl_z, b_abs, b_sub, b_cro, aver, vars, over, snap, prev, next, play, scal_b, scal_t
	common settings_common, px, py, pz, cut, abs_scale, show_cross, show_cuts, sub_aver, selected_cube, selected_overplot, selected_snapshot, af_x, af_y, af_z

	selected_cube = 0
	selected_overplot = 0
	abs_scale = 1
	sub_aver = 0
	show_cross = 1
	px = num_x / 2
	py = num_y / 2
	pz = num_z / 2
	pos_b = replicate (-1.0, num_cubes, 2)
	pos_t = replicate (-1.0, num_cubes, 2)

	prepare_cube, -1

	WIDGET_CONTROL, b_abs, SET_VALUE = abs_scale
	WIDGET_CONTROL, b_sub, SET_VALUE = sub_aver
	WIDGET_CONTROL, b_cro, SET_VALUE = show_cross
	WIDGET_CONTROL, sl_x, SET_VALUE = px
	WIDGET_CONTROL, sl_y, SET_VALUE = py
	WIDGET_CONTROL, sl_z, SET_VALUE = pz
	WIDGET_CONTROL, scal_b, SET_VALUE = pos_b[selected_cube,sub_aver]
	WIDGET_CONTROL, scal_t, SET_VALUE = pos_t[selected_cube,sub_aver]
	WIDGET_CONTROL, vars, SET_DROPLIST_SELECT = selected_cube
	WIDGET_CONTROL, over, SET_DROPLIST_SELECT = selected_overplot
	WIDGET_CONTROL, snap, SET_DROPLIST_SELECT = selected_snapshot

	if (num_snapshots ge 2) then prev_active = 1 else prev_active = 0
	next_active = 0
	WIDGET_CONTROL, prev, SENSITIVE = prev_active
	WIDGET_CONTROL, next, SENSITIVE = next_active
end


; Sophisticated interface with caching of VAR-files
pro cmp_cslice_cache, set_names, limits, units=units, coords=coords, scaling=scaling, overplots=overplots

	common varset_common, set, overplot, oversets, unit, coord, varsets, varfiles, sources
	common cslice_common, cube, field, num_cubes, num_overs, num_snapshots
	common event_common, button_pressed_yz, button_pressed_xz, button_pressed_xy
	common slider_common, bin_x, bin_y, bin_z, num_x, num_y, num_z, pos_b, pos_t, csmin, csmax, dimensionality
	common gui_common, wimg_yz, wimg_xz, wimg_xy, wcut_x, wcut_y, wcut_z, sl_x, sl_y, sl_z, b_abs, b_sub, b_cro, aver, vars, over, snap, prev, next, play, scal_b, scal_t
	common settings_common, px, py, pz, cut, abs_scale, show_cross, show_cuts, sub_aver, selected_cube, selected_overplot, selected_snapshot, af_x, af_y, af_z

	; DEFAULT SETTINGS:
	abs_scale = 1
	show_cross = 1
	show_cuts = 1
	sub_aver = 0
	selected_cube = 0
	selected_overplot = 0
	selected_snapshot = 0
	af_x = 0
	af_y = 0
	af_z = 0
	min_size = 8


	set = set_names
	if (n_elements (overplots) eq 0) then overplots = {none:'none'} else overplots = create_struct ({none:'none'}, overplots)
	overplot = overplots

	if (n_elements (units) ge 1) then unit = units
	if (n_elements (unit) eq 0) then begin
		print, "WARNING: setting units to unity."
		unit = { velocity:1, temperature:1, length:1, time:1, density:1 }
	end

	if (n_elements (scaling) eq 0) then scaling = 1
	if (n_elements (scaling) eq 1) then scaling = [ scaling, scaling, scaling ]

	; setup limits, if necessary
	if (n_elements (limits) eq 0) then begin
		dims = size (varsets.(0))
		nghost_x = (dims[1] - (size (coord.x))[1]) / 2
		nghost_y = (dims[2] - (size (coord.y))[1]) / 2
		nghost_z = (dims[3] - (size (coord.z))[1]) / 2
		num_x = dims[1] - 2*nghost_x
		num_y = dims[2] - 2*nghost_y
		num_z = dims[3] - 2*nghost_z
		limits = reform (nghost_x+spread(indgen(num_x),[1,2],[num_y,num_z]) + dims[1]*(nghost_y+spread(indgen(num_y),[0,2],[num_x,num_z])) + dims[1]*dims[2]*(nghost_z+spread(indgen(num_z),[0,1],[num_x,num_y])), num_x, num_y, num_z)
	end

	cut = reform (limits, (size (limits))[1], (size (limits))[2], (size (limits))[3])

	wimg_yz = !d.window
	wimg_xz = !d.window
	wimg_xy = !d.window

	wcut_x = !d.window
	wcut_y = !d.window
	wcut_z = !d.window

	scal_b = 0
	scal_t = 0

	num_snapshots = n_elements (varfiles)
	snaps = varfiles[*].title

	tags = tag_names (set)
	num_cubes = n_tags (set)

	overs = tag_names (overplot)
	num_overs = n_tags (overplot)

	pos_b = replicate (!VALUES.D_NAN, num_cubes, 2)
	pos_t = replicate (!VALUES.D_NAN, num_cubes, 2)

	
	prepare_set, 0
	prepare_cube, -1


	if (num_x*scaling[0] lt min_size) then scaling[0] = min_size / num_x
	if (num_y*scaling[1] lt min_size) then scaling[1] = min_size / num_y
	if (num_z*scaling[2] lt min_size) then scaling[2] = min_size / num_z

	bin_x = scaling[0]
	bin_y = scaling[1]
	bin_z = scaling[2]

	if (n_elements (coords) ge 1) then coord = coords
	if (n_elements (coord) eq 0) then begin
		print, "WARNING: setting the pixel size to unit length."
		coord = { x:findgen(num_x)*unit.length, y:findgen(num_y)*unit.length, z:findgen(num_z)*unit.length }
	end


	px = num_x / 2
	py = num_y / 2
	pz = num_z / 2


	if (num_cubes ge 2) then vars_active = 1 else vars_active = 0
	if (num_overs ge 2) then over_active = 1 else over_active = 0
	if (num_snapshots ge 2) then snap_active = 1 else snap_active = 0
	if (num_snapshots ge 2) then prev_active = 1 else prev_active = 0
	next_active = 0

	if (num_x gt 1) then sl_x_active = 1 else sl_x_active = 0
	if (num_y gt 1) then sl_y_active = 1 else sl_y_active = 0
	if (num_z gt 1) then sl_z_active = 1 else sl_z_active = 0

	dimensionality = sl_x_active + sl_y_active + sl_z_active
	if (dimensionality eq 0) then begin
		print, "Are you sure, you want to visualize 0D data?"
		print, "(If yes, just type '.continue' and press the return key.)"
		stop
	end

	button_pressed_yz = 0
	button_pressed_xz = 0
	button_pressed_xy = 0

	MOTHER	= WIDGET_BASE (title='compare cube-slices')
	BASE    = WIDGET_BASE (MOTHER, /col)
	TOP     = WIDGET_BASE (BASE, /row)
	scol    = WIDGET_BASE (top, /col)
	scot    = WIDGET_BASE (scol, /col)
	sl_x    = WIDGET_SLIDER (scot, uvalue='SLX', value=px, min=0, max=(num_x-1)>1, xsize=(num_x*bin_x>128)+10, /drag, sensitive=sl_x_active)
	scot    = WIDGET_BASE (scol, /col)
	sl_y    = WIDGET_SLIDER (scot, uvalue='SLY', value=py, min=0, max=(num_y-1)>1, xsize=(num_y*bin_y>128)+10, /drag, sensitive=sl_y_active)
	scot    = WIDGET_BASE (scol, /col)
	sl_z    = WIDGET_SLIDER (scot, uvalue='SLZ', value=pz, min=0, max=(num_z-1)>1, xsize=(num_z*bin_z>128)+10, /drag, sensitive=sl_z_active)
	bcol    = WIDGET_BASE (top, /col)
	b_abs   = CW_BGROUP (bcol, 'absolute scaling', /nonexcl, uvalue='SCL', set_value=abs_scale)
	b_sub   = CW_BGROUP (bcol, 'substract averages', /nonexcl, uvalue='SUB_AVER', set_value=sub_aver)
	b_cro   = CW_BGROUP (bcol, 'show crosshairs', /nonexcl, uvalue='CHS', set_value=show_cross)
	aver    = WIDGET_BUTTON (bcol, value='vertical profile', uvalue='SHOW_AVER')
	bcol    = WIDGET_BASE (top, /col)
	bcot    = WIDGET_BASE (bcol, /row)
	vars    = WIDGET_DROPLIST (bcot, value=tags, uvalue='VAR', sensitive=vars_active, EVENT_PRO=cslice_event, title='data set')
	bcot    = WIDGET_BASE (bcol, /row)
	over    = WIDGET_DROPLIST (bcot, value=overs, uvalue='OVER', sensitive=over_active, EVENT_PRO=cslice_event, title='overplot')
	bcot    = WIDGET_BASE (bcol, /row)
	snap    = WIDGET_DROPLIST (bcot, value=snaps, uvalue='SNAP', sensitive=snap_active, EVENT_PRO=cslice_event, title='time step')
	prev    = WIDGET_BUTTON (bcot, value='-', uvalue='PREV', sensitive=prev_active, EVENT_PRO=cslice_event)
	next    = WIDGET_BUTTON (bcot, value='+', uvalue='NEXT', sensitive=next_active, EVENT_PRO=cslice_event)
	bcol    = WIDGET_BASE (top, /col)
	tmp	= WIDGET_BUTTON (bcol, value='RESET', uvalue='RESET', xsize=100)
	tmp	= WIDGET_BUTTON (bcol, value='LOAD', uvalue='LOAD', xsize=100)
	tmp	= WIDGET_BUTTON (bcol, value='SAVE', uvalue='SAVE', xsize=100)
	play	= WIDGET_BUTTON (bcol, value='PLAY', uvalue='PLAY', xsize=100, sensitive=snap_active)
	tmp	= WIDGET_BUTTON (bcol, value='QUIT', uvalue='QUIT', xsize=100)
	drow    = WIDGET_BASE (BASE, /row)
	tmp     = WIDGET_DRAW (drow, UVALUE='DRAW_YZ', xsize=num_y*bin_y, ysize=num_z*bin_z, /button_events, /motion_events)
	WIDGET_CONTROL, tmp, /REALIZE
	wimg_yz = !d.window
	tmp     = WIDGET_DRAW (drow, UVALUE='DRAW_XZ', xsize=num_x*bin_x, ysize=num_z*bin_z, /button_events, /motion_events)
	WIDGET_CONTROL, tmp, /REALIZE
	wimg_xz = !d.window
	tmp     = WIDGET_DRAW (drow, UVALUE='DRAW_XY', xsize=num_x*bin_x, ysize=num_y*bin_y, /button_events, /motion_events)
	WIDGET_CONTROL, tmp, /REALIZE
	wimg_xy = !d.window
	MID     = WIDGET_BASE (BASE, /col)
	bcot    = WIDGET_BASE (MID, /row)

	scal_b = CW_FSLIDER (bcot, title='lower value (black level)', uvalue='IMG_SCLB', /double, /edit, min=csmin, max=csmax, drag=1, value=csmin, xsize=((2*num_x*bin_x+num_y*bin_y)/2>(500+max([num_x*bin_x,num_y*bin_y,num_z*bin_z]))/2)<500 )
	scal_t = CW_FSLIDER (bcot, title='upper value (white level)', uvalue='IMG_SCLT', /double, /edit, min=csmin, max=csmax, drag=1, value=csmax, xsize=((2*num_x*bin_x+num_y*bin_y)/2>(500+max([num_x*bin_x,num_y*bin_y,num_z*bin_z]))/2)<500 )

	WIDGET_CONTROL, MOTHER, /REALIZE
	wimg = !d.window

	cut_height = min([num_x*bin_x,num_y*bin_y,num_z*bin_z]) > 256
	LOW     = WIDGET_BASE (BASE, /row)
	tmp     = WIDGET_DRAW (LOW, UVALUE='CUT1', xsize=num_y*bin_y, ysize=cut_height)
	WIDGET_CONTROL, tmp, /REALIZE
	wcut_x  = !d.window
	tmp     = WIDGET_DRAW (LOW, UVALUE='CUT2', xsize=num_x*bin_x, ysize=cut_height)
	WIDGET_CONTROL, tmp, /REALIZE
	wcut_y  = !d.window
	tmp     = WIDGET_DRAW (LOW, UVALUE='CUT3', xsize=num_z*bin_z>128, ysize=cut_height)
	WIDGET_CONTROL, tmp, /REALIZE
	wcut_z  = !d.window

	WIDGET_CONTROL, BASE

	XMANAGER, "cslice", MOTHER, /no_block

	reset_GUI
	draw_images, 1, 1, 1
end

