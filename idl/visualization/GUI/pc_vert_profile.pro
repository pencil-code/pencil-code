;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_vert_profile.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;;  $Id$
;;;
;;;  Description:
;;;    Plots a vertical profile of the given 3D quantity.
;;;


; Event handling of vertical profile window
pro pc_vert_profile_event, event

	common vert_prof_common, z, t, prof_name, prof_mean, prof_min, prof_max, x_range, x_min, x_max, x_label, z_range, z_min, z_max, z_label, file_name
	common vert_prof_GUI_common, win, l_plot, l_line, plot_style, line_style, b_zero, b_line, b_log, show_zero, show_line, log_plot, sx_set, sz_set, sz_fr, z_coupled, sx_max, sx_min, sz_max, sz_min

	WIDGET_CONTROL, WIDGET_INFO (event.top, /CHILD)
	WIDGET_CONTROL, event.id, GET_UVALUE = eventval

	quit = -1
	DRAW_PROF = 0

	SWITCH eventval of
	'X_MIN': begin
		WIDGET_CONTROL, sx_min, GET_VALUE = val_min
		if (val_min gt x_range[1]) then begin
			val_min = x_range[1]
			WIDGET_CONTROL, sx_min, SET_VALUE = val_min
		end
		x_range[0] = val_min
		DRAW_PROF = 1
		break
	end
	'X_MAX': begin
		WIDGET_CONTROL, sx_max, GET_VALUE = val_max
		if (val_max lt x_range[0]) then begin
			val_max = x_range[0]
			WIDGET_CONTROL, sx_max, SET_VALUE = val_max
		end
		x_range[1] = val_max
		DRAW_PROF = 1
		break
	end
	'X_CENTER': begin
		x_width = x_range[1] - x_range[0]
		x_range = ([-x_width / 2.0, x_width / 2.0] < x_max) > x_min
		range_min = min (abs (x_range))
		x_range = (x_range < range_min) > (-range_min)
		x_width = x_range[1] - x_range[0]
		x_range = [-x_width / 2.0, x_width / 2.0]
		WIDGET_CONTROL, sx_min, SET_VALUE = x_range[0]
		WIDGET_CONTROL, sx_max, SET_VALUE = x_range[1]
		DRAW_PROF = 1
		break
	end
	'X_MINMAX': begin
		x_range = [x_min, x_max]
		WIDGET_CONTROL, sx_min, SET_VALUE = x_range[0]
		WIDGET_CONTROL, sx_max, SET_VALUE = x_range[1]
		DRAW_PROF = 1
		break
	end
	'Z_MINMAX': begin
		z_range = [z_min, z_max]
		WIDGET_CONTROL, sz_min, SET_VALUE = z_range[0]
		WIDGET_CONTROL, sz_max, SET_VALUE = z_range[1]
		DRAW_PROF = 1
		break
	end
	'Z_MIN': begin
		WIDGET_CONTROL, sz_min, GET_VALUE = val_min
		val_max = z_range[1]
		if (z_coupled lt 0) then begin
			if (val_min gt val_max) then begin
				val_min = val_max
				WIDGET_CONTROL, sz_min, SET_VALUE = val_min
			end
		end else begin
			if (val_min gt z_max - z_coupled) then begin
				val_min = z_max - z_coupled
				WIDGET_CONTROL, sz_min, SET_VALUE = val_min
			end
			val_max = val_min + z_coupled
			WIDGET_CONTROL, sz_max, SET_VALUE = val_max
			z_range[1] = val_max
		end
		z_range[0] = val_min
		DRAW_PROF = 1
		break
	end
	'Z_MAX': begin
		WIDGET_CONTROL, sz_max, GET_VALUE = val_max
		val_min = z_range[0]
		if (z_coupled lt 0) then begin
			if (val_max lt val_min) then begin
				val_max = val_min
				WIDGET_CONTROL, sz_max, SET_VALUE = val_max
			end
		end else begin
			if (val_max lt z_min + z_coupled) then begin
				val_max = z_min + z_coupled
				WIDGET_CONTROL, sz_max, SET_VALUE = val_max
			end
			val_min = val_max - z_coupled
			WIDGET_CONTROL, sz_min, SET_VALUE = val_min
			z_range[0] = val_min
		end
		z_range[1] = val_max
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
		pc_vert_profile_reset
		DRAW_PROF = 1
		break
	end
	'IMAGE': begin
		WIDGET_CONTROL, event.id, SENSITIVE = 0
		pc_save_image, file_name+".png", window=win
		save, z, t, prof_name, prof_mean, prof_min, prof_max, x_range, x_label, z_range, z_label, file=file_name+".xdr"
		WIDGET_CONTROL, event.id, SENSITIVE = 1
		break
	end
	'Z_COUPLE': begin
		WIDGET_CONTROL, sz_fr, set_value='<= RELEASE =>', set_uvalue='Z_RELEASE'
		z_coupled = z_range[1] - z_range[0]
		break
	end
	'Z_RELEASE': begin
		WIDGET_CONTROL, sz_fr, set_value='<= COUPLE =>', set_uvalue='Z_COUPLE'
		z_coupled = -1
		break
	end
	'QUIT': begin
		quit = event.top
		break
	end
	endswitch

	if (DRAW_PROF) then pc_vert_profile_draw

	WIDGET_CONTROL, sx_set, sensitive=((x_range[0] gt x_min) or (x_range[1] lt x_max))
	WIDGET_CONTROL, sz_set, sensitive=(((z_range[0] gt z_min) or (z_range[1] lt z_max)) and (z_coupled le -1.0))
	WIDGET_CONTROL, sz_fr, sensitive=((z_range[0] gt z_min) or (z_range[1] lt z_max))

	WIDGET_CONTROL, WIDGET_INFO (event.top, /CHILD)

	if (quit ge 0) then WIDGET_CONTROL, quit, /DESTROY

	return
end


; Draw the timeseries plots
pro pc_vert_profile_draw

	common vert_prof_common, z, t, prof_name, prof_mean, prof_min, prof_max, x_range, x_min, x_max, x_label, z_range, z_min, z_max, z_label, file_name
	common vert_prof_GUI_common, win, l_plot, l_line, plot_style, line_style, b_zero, b_line, b_log, show_zero, show_line, log_plot, sx_set, sz_set, sz_fr, z_coupled, sx_max, sx_min, sz_max, sz_min

	wset, win

	; increase charsize
	normal_charsize = !P.CHARSIZE
	if (normal_charsize le 0.0) then normal_charsize = 1.0
	!P.CHARSIZE = 1.25 * normal_charsize

	; setup horizontal plot range
	range = get_val_range (x_range)
	if (show_zero and (range[0] gt 0.0)) then range[0] = 0.0
	if (show_zero and (range[1] lt 0.0)) then range[1] = 0.0

	; plot profile mean value
	if (plot_style ge 3) then psym = 3 else psym = 0
	plot, prof_mean, z, xlog=log_plot, xs=3, ys=1, xr=range, yr=get_val_range (z_range), psym=psym, title=prof_name, xtitle=x_label, ytitle=z_label
	if (show_line) then oplot, [0.0, 0.0], minmax (z), linestyle=1, color=20020
	if (plot_style le 2) then oplot, prof_mean, z
	if (plot_style gt 0) then begin
		psym = 3
		color = 200
		if ((plot_style eq 2) or (plot_style eq 4)) then psym = 2
		if (plot_style ge 3) then color = -1
		oplot, prof_mean, z, psym=psym, color=color
	end

	; plot profile minimum and maximum values
	if (line_style gt 0) then begin
		linestyle = 0
		if (line_style le 3) then linestyle = 3 - line_style
		oplot, prof_min, z, linestyle=linestyle
		oplot, prof_max, z, linestyle=linestyle
	end

	; reset charsize
	!P.CHARSIZE = normal_charsize
end


; Reset to defaults
pro pc_vert_profile_reset

	common vert_prof_common, z, t, prof_name, prof_mean, prof_min, prof_max, x_range, x_min, x_max, x_label, z_range, z_min, z_max, z_label
	common vert_prof_GUI_common, win, l_plot, l_line, plot_style, line_style, b_zero, b_line, b_log, show_zero, show_line, log_plot, sx_set, sz_set, sz_fr, z_coupled, sx_max, sx_min, sz_max, sz_min

	; initial GUI settings
	x_range = [x_min, x_max]
	z_range = [z_min, z_max]
	z_coupled = -1
	plot_style = 1
	line_style = 1
	show_zero = 0
	show_line = 1
	log_plot = 0

	if (win ge 0) then begin
		WIDGET_CONTROL, l_plot, SET_DROPLIST_SELECT = plot_style
		WIDGET_CONTROL, l_line, SET_DROPLIST_SELECT = line_style
		WIDGET_CONTROL, b_zero, SET_VALUE = show_zero
		WIDGET_CONTROL, b_line, SET_VALUE = show_line
		WIDGET_CONTROL, b_log, SET_VALUE = log_plot
		WIDGET_CONTROL, b_log, sensitive = (x_min gt 0.0)
		WIDGET_CONTROL, sx_min, SET_VALUE = x_range[0]
		WIDGET_CONTROL, sx_max, SET_VALUE = x_range[1]
		WIDGET_CONTROL, sz_min, SET_VALUE = z_range[0]
		WIDGET_CONTROL, sz_max, SET_VALUE = z_range[1]
		pc_vert_profile_draw
	end
end


; Calculate and draw a vertical profile of the given 3D data
;
; data:        3D data cube (can be including ghost cells)
; coord:       coordinates for the vertical position of data in the cube,
;              asumed to be in the center of the data cube (eg. without ghost cells).
;              If omitted, the index numbers are used as coordinates.
; title:       title string for the plot
; min:         initial minimum value for data display
; max:         initial maximum value for data display
; log:         set this to use a logarithmic scale for data display
; horiz_label: label string for the horizontal axis
; vert_label:  label string for the vertical axis
; file_label:  label string for filenames (special characters will be filtered)
; time:        timestamp of the displayed data
;
pro pc_vert_profile, data, coord=coord, title=title, horiz_label=horiz_label, vert_label=vert_label, min=min, max=max, log=log, file_label=file_label, time=time

	common vert_prof_common, z, t, prof_name, prof_mean, prof_min, prof_max, x_range, x_min, x_max, x_label, z_range, z_min, z_max, z_label, file_name
	common vert_prof_GUI_common, win, l_plot, l_line, plot_style, line_style, b_zero, b_line, b_log, show_zero, show_line, log_plot, sx_set, sz_set, sz_fr, z_coupled, sx_max, sx_min, sz_max, sz_min

	num = (size (data))[3]
	if (n_elements (coord) eq 0) then coord = findgen (num)
	num_coord = size (coord, /n_elements)
	z = reform (coord, num_coord)

	t = "N/A"
	if (keyword_set (time)) then t = time

	prof_name = ""
	if (keyword_set (title)) then prof_name = title

	if (keyword_set (file_label)) then file_name = file_label else file_name = prof_name
	if (file_name eq "") then file_name = "vert_profile"
	pos = stregex (file_name, '[ \/]+', length=len)
	while (pos gt 0) do begin
		file_name = strmid (file_name, 0, pos) + "_" + strmid (file_name, pos+len)
		pos = stregex (file_name, '[ \/]+', length=len)
	end
	file_label = file_name

	; if the data contains ghost cells, but not the coordinates, center the plot:
	start_pos = (num - num_coord) / 2

	prof_mean = dblarr (num)
	prof_min = dblarr (num)
	prof_max = dblarr (num)

	for pos = 0, num_coord - 1 do begin
		prof_mean[pos] = mean (data[*,*,start_pos + pos], /double)
		tmp = minmax (data[*,*,start_pos + pos])
		prof_min[pos] = tmp[0]
		prof_max[pos] = tmp[1]
	end

	x_min = min (prof_min)
	x_max = max (prof_max)
	x_range = [x_min, x_max]
	if (keyword_set (horiz_label)) then x_label = horiz_label else x_label = ""

	; extend vertical profile plot range by half a grid distance
	if (num_coord le 1) then begin
		z_range = [coord[0]-1, coord[0]+1]
	end else begin
		z_range = [2*coord[0]-coord[1], 2*coord[num_coord-1]-coord[num_coord-2]]
	end
	z_min = z_range[0]
	z_max = z_range[1]
	if (keyword_set (vert_label)) then z_label = vert_label else z_label = ""

	; GUI default settings
	win = -1
	plot_width = 600
	plot_height = 500
	sl_width = (plot_width - 100) / 2
	pc_vert_profile_reset

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
	l_line	= WIDGET_DROPLIST (BUT, value=['none', 'dashed', 'dotted', 'solid'], uvalue='LINESTYLE', title='plot min/max:')

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

	range = get_val_range (x_range)
	tmp	= WIDGET_LABEL (PLOTS, value='horizontal axis:')
	BUT	= WIDGET_BASE (PLOTS, /row)
	sx_min	= CW_FSLIDER (BUT, xsize=sl_width, title='minimum value', uvalue='X_MIN', /double, /edit, min=range[0], max=range[1], drag=1, value=x_range[0])
	CTRL	= WIDGET_BASE (BUT, /col, frame=0)
	sx_set	= WIDGET_BUTTON (CTRL, value='<= MINMAX =>', uvalue='X_MINMAX', sensitive=0)
	tmp	= WIDGET_BUTTON (CTRL, value='<= CENTER =>', uvalue='X_CENTER', sensitive=(x_min lt 0.0))
	sx_max	= CW_FSLIDER (BUT, xsize=sl_width, title='maximum value', uvalue='X_MAX', /double, /edit, min=range[0], max=range[1], drag=1, value=x_range[1])

	range = get_val_range (z_range)
	tmp	= WIDGET_LABEL (PLOTS, value='vertical axis:')
	BUT	= WIDGET_BASE (PLOTS, /row)
	sz_min	= CW_FSLIDER (BUT, xsize=sl_width, title='minimum value', uvalue='Z_MIN', /double, /edit, min=range[0], max=range[1], drag=1, value=z_range[0])
	CTRL	= WIDGET_BASE (BUT, /col, frame=0)
	sz_set	= WIDGET_BUTTON (CTRL, value='<= MINMAX =>', uvalue='Z_MINMAX', sensitive=0)
	sz_fr	= WIDGET_BUTTON (CTRL, value='<= COUPLE =>', uvalue='Z_COUPLE', sensitive=0)
	sz_max	= CW_FSLIDER (BUT, xsize=sl_width, title='maximum value', uvalue='Z_MAX', /double, /edit, min=range[0], max=range[1], drag=1, value=z_range[1])

	BASE	= WIDGET_BASE (APP, /row)


	WIDGET_CONTROL, MOTHER, /REALIZE
	WIDGET_CONTROL, d_prof, GET_VALUE = win

	WIDGET_CONTROL, BASE

	pc_vert_profile_reset

	if (keyword_set (log)) then log_plot = (log ne 0) and (x_min gt 0.0)
	if (keyword_set (min)) then x_range[0] = min
	if (keyword_set (max)) then x_range[1] = max
	WIDGET_CONTROL, b_log, SET_VALUE = log_plot
	WIDGET_CONTROL, b_log, sensitive = (x_min gt 0.0)
	WIDGET_CONTROL, sx_min, SET_VALUE = x_range[0]
	WIDGET_CONTROL, sx_max, SET_VALUE = x_range[1]

	XMANAGER, "pc_vert_profile", MOTHER, /no_block

	pc_vert_profile_draw

end

