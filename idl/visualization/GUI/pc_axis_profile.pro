;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_axis_profile.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;;  $Id$
;;;
;;;  Description:
;;;    Plots a vertical profile of the given 3D quantity.
;;;


; Event handling of vertical profile window
pro pc_axis_profile_event, event

	common axis_prof_common, coords, num, num_coord, axis, t, prof_name, prof_mean, prof_min, prof_max, a_range, a_min, a_max, a_label, v_range, v_min, v_max, v_label, file_name
	common axis_prof_GUI_common, win, l_plot, l_line, plot_style, line_style, b_zero, b_line, b_log, show_zero, show_line, log_plot, sa_set, sv_set, sv_fr, v_coupled, sa_max, sa_min, sv_max, sv_min

	WIDGET_CONTROL, WIDGET_INFO (event.top, /CHILD)
	WIDGET_CONTROL, event.id, GET_UVALUE = eventval

	quit = -1
	DRAW_PROF = 0

	SWITCH eventval of
	'A_MIN': begin
		WIDGET_CONTROL, sa_min, GET_VALUE = val_min
		if (val_min gt a_range[1]) then begin
			val_min = a_range[1]
			WIDGET_CONTROL, sa_min, SET_VALUE = val_min
		end
		a_range[0] = val_min
		DRAW_PROF = 1
		break
	end
	'A_MAX': begin
		WIDGET_CONTROL, sa_max, GET_VALUE = val_max
		if (val_max lt a_range[0]) then begin
			val_max = a_range[0]
			WIDGET_CONTROL, sa_max, SET_VALUE = val_max
		end
		a_range[1] = val_max
		DRAW_PROF = 1
		break
	end
	'A_CENTER': begin
		a_width = a_range[1] - a_range[0]
		a_range = ([-a_width / 2.0, a_width / 2.0] < a_max) > a_min
		range_min = min (abs (a_range))
		a_range = (a_range < range_min) > (-range_min)
		a_width = a_range[1] - a_range[0]
		a_range = [-a_width / 2.0, a_width / 2.0]
		WIDGET_CONTROL, sa_min, SET_VALUE = a_range[0]
		WIDGET_CONTROL, sa_max, SET_VALUE = a_range[1]
		DRAW_PROF = 1
		break
	end
	'A_MINMAX': begin
		a_range = [a_min, a_max]
		WIDGET_CONTROL, sa_min, SET_VALUE = a_range[0]
		WIDGET_CONTROL, sa_max, SET_VALUE = a_range[1]
		DRAW_PROF = 1
		break
	end
	'V_MINMAX': begin
		v_range = [v_min, v_max]
		WIDGET_CONTROL, sv_min, SET_VALUE = v_range[0]
		WIDGET_CONTROL, sv_max, SET_VALUE = v_range[1]
		DRAW_PROF = 1
		break
	end
	'V_MIN': begin
		WIDGET_CONTROL, sv_min, GET_VALUE = val_min
		val_max = v_range[1]
		if (v_coupled lt 0) then begin
			if (val_min gt val_max) then begin
				val_min = val_max
				WIDGET_CONTROL, sv_min, SET_VALUE = val_min
			end
		end else begin
			if (val_min gt v_max - v_coupled) then begin
				val_min = v_max - v_coupled
				WIDGET_CONTROL, sv_min, SET_VALUE = val_min
			end
			val_max = val_min + v_coupled
			WIDGET_CONTROL, sv_max, SET_VALUE = val_max
			v_range[1] = val_max
		end
		v_range[0] = val_min
		DRAW_PROF = 1
		break
	end
	'V_MAX': begin
		WIDGET_CONTROL, sv_max, GET_VALUE = val_max
		val_min = v_range[0]
		if (v_coupled lt 0) then begin
			if (val_max lt val_min) then begin
				val_max = val_min
				WIDGET_CONTROL, sv_max, SET_VALUE = val_max
			end
		end else begin
			if (val_max lt v_min + v_coupled) then begin
				val_max = v_min + v_coupled
				WIDGET_CONTROL, sv_max, SET_VALUE = val_max
			end
			val_min = val_max - v_coupled
			WIDGET_CONTROL, sv_min, SET_VALUE = val_min
			v_range[0] = val_min
		end
		v_range[1] = val_max
		DRAW_PROF = 1
		break
	end
	'STYLE': begin
		if (plot_style ne event.index) then begin
			plot_style = event.index
			DRAW_PROF = 1
		end
		break
	end
	'LINESTYLE': begin
		if (line_style ne event.index) then begin
			line_style = event.index
			DRAW_PROF = 1
		end
		break
	end
	'ZERO': begin
		show_zero = event.select
		DRAW_PROF = 1
		break
	end
	'LINE': begin
		show_line = event.select
		DRAW_PROF = 1
		break
	end
	'LOG_PLOT': begin
		log_plot = event.select
		if (log_plot) then begin
			show_zero = 0
			WIDGET_CONTROL, b_zero, SET_VALUE = show_zero
		end
		DRAW_PROF = 1
		break
	end
	'RESET': begin
		pc_axis_profile_reset
		DRAW_PROF = 1
		break
	end
	'IMAGE': begin
		WIDGET_CONTROL, event.id, SENSITIVE = 0
		pc_save_image, file_name+".png", window=win
		save, coords, num, num_coord, axis, t, prof_name, prof_mean, prof_min, prof_max, a_range, a_min, a_max, a_label, v_range, v_min, v_max, v_label, file_name, file=file_name+".xdr"
		WIDGET_CONTROL, event.id, SENSITIVE = 1
		break
	end
	'V_COUPLE': begin
		WIDGET_CONTROL, sv_fr, set_value='<= RELEASE =>', set_uvalue='V_RELEASE'
		v_coupled = v_range[1] - v_range[0]
		break
	end
	'V_RELEASE': begin
		WIDGET_CONTROL, sv_fr, set_value='<= COUPLE =>', set_uvalue='V_COUPLE'
		v_coupled = -1
		break
	end
	'QUIT': begin
		quit = event.top
		break
	end
	endswitch

	if (DRAW_PROF) then pc_axis_profile_draw

	WIDGET_CONTROL, sa_set, sensitive=((a_range[0] gt a_min) or (a_range[1] lt a_max))
	WIDGET_CONTROL, sv_set, sensitive=(((v_range[0] gt v_min) or (v_range[1] lt v_max)) and (v_coupled le -1.0))
	WIDGET_CONTROL, sv_fr, sensitive=((v_range[0] gt v_min) or (v_range[1] lt v_max))

	WIDGET_CONTROL, WIDGET_INFO (event.top, /CHILD)

	if (quit ge 0) then WIDGET_CONTROL, quit, /DESTROY

	return
end


; Draw the timeseries plots
pro pc_axis_profile_draw

	common axis_prof_common, coords, num, num_coord, axis, t, prof_name, prof_mean, prof_min, prof_max, a_range, a_min, a_max, a_label, v_range, v_min, v_max, v_label, file_name
	common axis_prof_GUI_common, win, l_plot, l_line, plot_style, line_style, b_zero, b_line, b_log, show_zero, show_line, log_plot, sa_set, sv_set, sv_fr, v_coupled, sa_max, sa_min, sv_max, sv_min

	wset, win

	; increase charsize
	normal_charsize = !P.CHARSIZE
	if (normal_charsize le 0.0) then normal_charsize = 1.0
	!P.CHARSIZE = 1.25 * normal_charsize

	; setup horizontal plot range
	range = get_val_range (a_range)
	if (show_zero and (range[0] gt 0.0)) then range[0] = 0.0
	if (show_zero and (range[1] lt 0.0)) then range[1] = 0.0

	start_pos = 0
	; if the coordinates contain ghost cells, but not the data, center the plot:
	if (num_coord gt num) then start_pos = (num_coord - num) / 2
	end_pos = start_pos + num - 1

	; plot profile mean value
	if (plot_style ge 3) then psym = 3 else psym = 0
	plot, prof_mean, coords[start_pos:end_pos], xlog=log_plot, xs=3, ys=1, xr=range, yr=get_val_range (v_range), psym=psym, title=prof_name, xtitle=a_label, ytitle=v_label
	if (show_line) then oplot, [0.0, 0.0], minmax (coords), linestyle=1, color=20020
	if (plot_style le 2) then oplot, prof_mean, coords[start_pos:end_pos]
	if (plot_style gt 0) then begin
		psym = 3
		color = 200
		if ((plot_style eq 2) or (plot_style eq 4)) then psym = 2
		if (plot_style ge 3) then color = -1
		oplot, prof_mean, coords[start_pos:end_pos], psym=psym, color=color
	end

	; plot profile minimum and maximum values
	if (line_style gt 0) then begin
		sym_type = 0
		if (line_style eq 4) then sym_type = 1
		if (num eq 1) then color = 200 else color = -1
		line_type = 0
		if (line_style le 3) then line_type = 3 - line_style
		oplot, prof_min, coords[start_pos:end_pos], linestyle=line_type, psym=sym_type, color=color
		oplot, prof_max, coords[start_pos:end_pos], linestyle=line_type, psym=sym_type, color=color
	end

	; reset charsize
	!P.CHARSIZE = normal_charsize
end


; Reset to defaults
pro pc_axis_profile_reset

	common axis_prof_common, coords, num, num_coord, axis, t, prof_name, prof_mean, prof_min, prof_max, a_range, a_min, a_max, a_label, v_range, v_min, v_max, v_label, file_name
	common axis_prof_GUI_common, win, l_plot, l_line, plot_style, line_style, b_zero, b_line, b_log, show_zero, show_line, log_plot, sa_set, sv_set, sv_fr, v_coupled, sa_max, sa_min, sv_max, sv_min

	; initial GUI settings
	a_range = [a_min, a_max]
	v_range = [v_min, v_max]
	v_coupled = -1
	if (num gt 1) then plot_style = 1 else plot_style = 4
	if (num gt 1) then line_style = 1 else line_style = 4
	show_zero = 0
	show_line = 1
	log_plot = 0

	if (win ge 0) then begin
		WIDGET_CONTROL, l_plot, SET_DROPLIST_SELECT = plot_style
		WIDGET_CONTROL, l_line, SET_DROPLIST_SELECT = line_style
		WIDGET_CONTROL, b_zero, SET_VALUE = show_zero
		WIDGET_CONTROL, b_line, SET_VALUE = show_line
		WIDGET_CONTROL, b_log, SET_VALUE = log_plot
		WIDGET_CONTROL, b_log, sensitive = (a_min gt 0.0)
		WIDGET_CONTROL, sa_min, SET_VALUE = a_range[0]
		WIDGET_CONTROL, sa_max, SET_VALUE = a_range[1]
		WIDGET_CONTROL, sv_min, SET_VALUE = v_range[0]
		WIDGET_CONTROL, sv_max, SET_VALUE = v_range[1]
		pc_axis_profile_draw
	end
end


; Calculate and draw a vertical profile of the given 3D data
;
; axis_dir:    direction of the axis for profiling (0-2 or 'X', 'Y', 'Z')
; data:        3D data cube (can be including ghost cells)
; coord:       coordinates for the vertical position of data in the cube,
;              asumed to be in the center of the data cube (eg. without ghost cells).
;              If omitted, the index numbers are used as coordinates.
; title:       title string for the plot
; min:         initial minimum value for data display
; max:         initial maximum value for data display
; log:         set this to use a logarithmic scale for data display
; horiz_label: label string for the axis (horizontal)
; vert_label:  label string for the value (vertical)
; file_label:  label string for filenames (special characters will be filtered)
; time:        timestamp of the displayed data
;
pro pc_axis_profile, axis_dir, data, coord=coord, title=title, horiz_label=horiz_label, vert_label=vert_label, min=min, max=max, log=log, file_label=file_label, time=time

	common axis_prof_common, coords, num, num_coord, axis, t, prof_name, prof_mean, prof_min, prof_max, a_range, a_min, a_max, a_label, v_range, v_min, v_max, v_label, file_name
	common axis_prof_GUI_common, win, l_plot, l_line, plot_style, line_style, b_zero, b_line, b_log, show_zero, show_line, log_plot, sa_set, sv_set, sv_fr, v_coupled, sa_max, sa_min, sv_max, sv_min

	if (size (axis_dir, /type) eq 7) then begin
		if (strupcase (axis_dir) eq 'X') then axis = 0
		if (strupcase (axis_dir) eq 'Y') then axis = 1
		if (strupcase (axis_dir) eq 'Z') then axis = 2
		message, "pc_axis_profile: unknown axis direction '"+axis_dir+"'."
	end else begin
		axis = round (axis_dir)
	end

	num_dimensions = size (data, /n_dimensions)
	if (axis ge num_dimensions) then message, "pc_axis_profile: axis direction unavailable."

	num = (size (data, /dimensions))[axis]
	if (n_elements (coord) eq 0) then coord = findgen (num)
	num_coord = size (coord, /n_elements)
	coords = reform (coord, num_coord)

	t = "N/A"
	if (keyword_set (time)) then t = time

	prof_name = ""
	if (keyword_set (title)) then prof_name = title

	if (keyword_set (file_label)) then file_name = file_label else file_name = prof_name
	if (file_name eq "") then file_name = "axis_profile"
	pos = stregex (file_name, '[ \/]+', length=len)
	while (pos gt 0) do begin
		file_name = strmid (file_name, 0, pos) + "_" + strmid (file_name, pos+len)
		pos = stregex (file_name, '[ \/]+', length=len)
	end
	file_label = file_name

	start_pos = 0
	; if the data contains ghost cells, but not the coordinates, center the plot:
	if (num gt num_coord) then start_pos = (num - num_coord) / 2

	prof_mean = dblarr (num)
	prof_min = dblarr (num)
	prof_max = dblarr (num)

	for pos = 0, num - 1 do begin
		if (axis eq 0) then begin
			prof_mean[pos] = mean (data[start_pos + pos,*,*], /double)
			tmp = minmax (data[start_pos + pos,*,*])
		end else if (axis eq 1) then begin
			prof_mean[pos] = mean (data[*,start_pos + pos,*], /double)
			tmp = minmax (data[*,start_pos + pos,*])
		end else if (axis eq 2) then begin
			prof_mean[pos] = mean (data[*,*,start_pos + pos], /double)
			tmp = minmax (data[*,*,start_pos + pos])
		end
		prof_min[pos] = tmp[0]
		prof_max[pos] = tmp[1]
	end

	a_min = min (prof_min)
	a_max = max (prof_max)
	a_range = [a_min, a_max]
	if (keyword_set (horiz_label)) then a_label = horiz_label else a_label = ""

	; extend vertical profile plot range by half a grid distance
	if (num_coord le 1) then begin
		v_range = [coord[0]-1, coord[0]+1]
	end else begin
		v_range = [2*coord[0]-coord[1], 2*coord[num_coord-1]-coord[num_coord-2]]
	end
	v_min = v_range[0]
	v_max = v_range[1]
	if (keyword_set (vert_label)) then v_label = vert_label else v_label = ""

	; GUI default settings
	win = -1
	plot_width = 600
	plot_height = 500
	sl_width = (plot_width - 100) / 2
	pc_axis_profile_reset

	if (prof_name eq "") then add = "" else add = " of "
	MOTHER	= WIDGET_BASE (title="PC vertical profile analysis"+add+prof_name)
	APP	= WIDGET_BASE (MOTHER, /col)

	BASE	= WIDGET_BASE (APP, /row)

	CTRL	= WIDGET_BASE (BASE, /row)

	BUT	= WIDGET_BASE (CTRL, /col)
	tmp	= WIDGET_BUTTON (BUT, xsize=100, value='RESET', uvalue='RESET')
	tmp	= WIDGET_BUTTON (BUT, xsize=100, value='SAVE PLOT', uvalue='IMAGE')
	tmp	= WIDGET_BUTTON (BUT, xsize=100, value='QUIT', uvalue='QUIT')

	BUT	= WIDGET_BASE (CTRL, /col, frame=1, /align_left)
	l_plot	= WIDGET_DROPLIST (BUT, value=['line', 'line+dots', 'line+stars', 'dots', 'stars'], uvalue='STYLE', title='plot mean as:')
	l_line	= WIDGET_DROPLIST (BUT, value=['none', 'dashed', 'dotted', 'solid', '+'], uvalue='LINESTYLE', title='plot min/max:')

	BUT	= WIDGET_BASE (CTRL, /col)
	b_log	= CW_BGROUP (BUT, 'logarithmic plot', /nonexcl, uvalue='LOG_PLOT', set_value=log_plot)
	b_zero	= CW_BGROUP (BUT, 'extend to zero', /nonexcl, uvalue='ZERO', set_value=show_zero)
	b_line	= CW_BGROUP (BUT, 'plot zero line', /nonexcl, uvalue='LINE', set_value=show_line)

	tmp	= WIDGET_BASE (BASE, /row)
	BUT	= WIDGET_BASE (tmp, /col)

	BASE	= WIDGET_BASE (APP, /row)

	PLOTS	= WIDGET_BASE (BASE, /col)
	tmp	= WIDGET_BASE (PLOTS, /row)
	d_prof	= WIDGET_DRAW (tmp, xsize=plot_width, ysize=plot_height, retain=2)

	range = get_val_range (a_range)
	tmp	= WIDGET_LABEL (PLOTS, value='horizontal axis:')
	BUT	= WIDGET_BASE (PLOTS, /row)
	sa_min	= CW_FSLIDER (BUT, xsize=sl_width, title='minimum value', uvalue='A_MIN', /double, /edit, min=range[0], max=range[1], drag=1, value=a_range[0])
	CTRL	= WIDGET_BASE (BUT, /col, frame=0)
	sa_set	= WIDGET_BUTTON (CTRL, value='<= MINMAX =>', uvalue='A_MINMAX', sensitive=0)
	tmp	= WIDGET_BUTTON (CTRL, value='<= CENTER =>', uvalue='A_CENTER', sensitive=(a_min lt 0.0))
	sa_max	= CW_FSLIDER (BUT, xsize=sl_width, title='maximum value', uvalue='A_MAX', /double, /edit, min=range[0], max=range[1], drag=1, value=a_range[1])

	range = get_val_range (v_range)
	tmp	= WIDGET_LABEL (PLOTS, value='vertical axis:')
	BUT	= WIDGET_BASE (PLOTS, /row)
	sv_min	= CW_FSLIDER (BUT, xsize=sl_width, title='minimum value', uvalue='V_MIN', /double, /edit, min=range[0], max=range[1], drag=1, value=v_range[0])
	CTRL	= WIDGET_BASE (BUT, /col, frame=0)
	sv_set	= WIDGET_BUTTON (CTRL, value='<= MINMAX =>', uvalue='V_MINMAX', sensitive=0)
	sv_fr	= WIDGET_BUTTON (CTRL, value='<= COUPLE =>', uvalue='V_COUPLE', sensitive=0)
	sv_max	= CW_FSLIDER (BUT, xsize=sl_width, title='maximum value', uvalue='V_MAX', /double, /edit, min=range[0], max=range[1], drag=1, value=v_range[1])

	BASE	= WIDGET_BASE (APP, /row)


	WIDGET_CONTROL, MOTHER, /REALIZE
	WIDGET_CONTROL, d_prof, GET_VALUE = win

	WIDGET_CONTROL, BASE

	pc_axis_profile_reset

	if (keyword_set (log)) then log_plot = (log ne 0) and (a_min gt 0.0)
	if (keyword_set (min)) then a_range[0] = min
	if (keyword_set (max)) then a_range[1] = max
	WIDGET_CONTROL, b_log, SET_VALUE = log_plot
	WIDGET_CONTROL, b_log, sensitive = (a_min gt 0.0)
	WIDGET_CONTROL, sa_min, SET_VALUE = a_range[0]
	WIDGET_CONTROL, sa_max, SET_VALUE = a_range[1]

	XMANAGER, "pc_axis_profile", MOTHER, /no_block

	pc_axis_profile_draw

end

