;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_radial_profile.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;;  $Id$
;;;
;;;  Description:
;;;    Plots a radial profile of the given 3D quantity.
;;;


; Event handling of radial profile window
pro pc_radial_profile_event, event

	common radial_prof_common, coords, axis_dir, mid, num, t, prof_name, prof, a_range, a_min, a_max, a_label, r_range, r_min, r_max, r_label, file_name
	common radial_prof_GUI_common, win, l_plot, plot_style, b_zero, b_line, b_log, show_zero, show_line, log_plot, sa_set, sr_set, sr_fr, r_coupled, sa_max, sa_min, sr_max, sr_min

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
	'R_MINMAX': begin
		r_range = [r_min, r_max]
		WIDGET_CONTROL, sr_min, SET_VALUE = r_range[0]
		WIDGET_CONTROL, sr_max, SET_VALUE = r_range[1]
		DRAW_PROF = 1
		break
	end
	'R_MIN': begin
		WIDGET_CONTROL, sr_min, GET_VALUE = val_min
		val_max = r_range[1]
		if (r_coupled lt 0) then begin
			if (val_min gt val_max) then begin
				val_min = val_max
				WIDGET_CONTROL, sr_min, SET_VALUE = val_min
			end
		end else begin
			if (val_min gt r_max - r_coupled) then begin
				val_min = r_max - r_coupled
				WIDGET_CONTROL, sr_min, SET_VALUE = val_min
			end
			val_max = val_min + r_coupled
			WIDGET_CONTROL, sr_max, SET_VALUE = val_max
			r_range[1] = val_max
		end
		r_range[0] = val_min
		DRAW_PROF = 1
		break
	end
	'R_MAX': begin
		WIDGET_CONTROL, sr_max, GET_VALUE = val_max
		val_min = r_range[0]
		if (r_coupled lt 0) then begin
			if (val_max lt val_min) then begin
				val_max = val_min
				WIDGET_CONTROL, sr_max, SET_VALUE = val_max
			end
		end else begin
			if (val_max lt r_min + r_coupled) then begin
				val_max = r_min + r_coupled
				WIDGET_CONTROL, sr_max, SET_VALUE = val_max
			end
			val_min = val_max - r_coupled
			WIDGET_CONTROL, sr_min, SET_VALUE = val_min
			r_range[0] = val_min
		end
		r_range[1] = val_max
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
		pc_radial_profile_reset
		DRAW_PROF = 1
		break
	end
	'IMAGE': begin
		WIDGET_CONTROL, event.id, SENSITIVE = 0
		pc_save_image, file_name+"_radial-profile.png", window=win
		save, coords, num, t, prof_name, prof, a_range, a_min, a_max, a_label, r_range, r_min, r_max, r_label, file_name, file=file_name+"_radial-profile.xdr"
		WIDGET_CONTROL, event.id, SENSITIVE = 1
		break
	end
	'R_COUPLE': begin
		WIDGET_CONTROL, sr_fr, set_value='<= RELEASE =>', set_uvalue='R_RELEASE'
		r_coupled = r_range[1] - r_range[0]
		break
	end
	'R_RELEASE': begin
		WIDGET_CONTROL, sr_fr, set_value='<= COUPLE =>', set_uvalue='R_COUPLE'
		r_coupled = -1
		break
	end
	'QUIT': begin
		quit = event.top
		break
	end
	endswitch

	if (DRAW_PROF) then pc_radial_profile_draw

	WIDGET_CONTROL, sa_set, sensitive=((a_range[0] gt a_min) or (a_range[1] lt a_max))
	WIDGET_CONTROL, sr_set, sensitive=(((r_range[0] gt r_min) or (r_range[1] lt r_max)) and (r_coupled le -1.0))
	WIDGET_CONTROL, sr_fr, sensitive=((r_range[0] gt r_min) or (r_range[1] lt r_max))

	WIDGET_CONTROL, WIDGET_INFO (event.top, /CHILD)

	if (quit ge 0) then WIDGET_CONTROL, quit, /DESTROY

	return
end


; Draw the profile
pro pc_radial_profile_draw

	common radial_prof_common, coords, axis_dir, mid, num, t, prof_name, prof, a_range, a_min, a_max, a_label, r_range, r_min, r_max, r_label, file_name
	common radial_prof_GUI_common, win, l_plot, plot_style, b_zero, b_line, b_log, show_zero, show_line, log_plot, sa_set, sr_set, sr_fr, r_coupled, sa_max, sa_min, sr_max, sr_min

	wset, win

	; increase charsize
	normal_charsize = !P.CHARSIZE
	if (normal_charsize le 0.0) then normal_charsize = 1.0
	!P.CHARSIZE = 1.25 * normal_charsize

	; setup horizontal plot range
	range = get_val_range (a_range)
	if (show_zero and (range[0] gt 0.0)) then range[0] = 0.0
	if (show_zero and (range[1] lt 0.0)) then range[1] = 0.0

	; plot profile mean value
	if (plot_style ge 3) then psym = 3 else psym = 0
	plot, coords.r, prof, ylog=log_plot, xs=3, ys=1, xr=get_val_range (r_range), yr=range, psym=psym, title=prof_name, xtitle=a_label, ytitle=r_label
	if (show_line) then oplot, [0.0, 0.0], minmax (prof), linestyle=1, color=20020
	if (plot_style le 2) then oplot, coords.r, prof
	if (plot_style gt 0) then begin
		psym = 3
		color = 200
		if ((plot_style eq 2) or (plot_style eq 4)) then psym = 2
		if (plot_style ge 3) then color = -1
		oplot, coords.r, prof, psym=psym, color=color
	end

	; reset charsize
	!P.CHARSIZE = normal_charsize
end


; Reset to defaults
pro pc_radial_profile_reset

	common radial_prof_common, coords, axis_dir, mid, num, t, prof_name, prof, a_range, a_min, a_max, a_label, r_range, r_min, r_max, r_label, file_name
	common radial_prof_GUI_common, win, l_plot, plot_style, b_zero, b_line, b_log, show_zero, show_line, log_plot, sa_set, sr_set, sr_fr, r_coupled, sa_max, sa_min, sr_max, sr_min

	; initial GUI settings
	a_range = [a_min, a_max]
	r_range = [r_min, r_max]
	r_coupled = -1
	if (num gt 1) then plot_style = 2 else plot_style = 4
	if (num ge 128) then plot_style = 1
	show_zero = 1
	show_line = 1
	log_plot = 0

	if (win ge 0) then begin
		WIDGET_CONTROL, l_plot, SET_DROPLIST_SELECT = plot_style
		WIDGET_CONTROL, b_zero, SET_VALUE = show_zero
		WIDGET_CONTROL, b_line, SET_VALUE = show_line
		WIDGET_CONTROL, b_log, SET_VALUE = log_plot
		WIDGET_CONTROL, b_log, sensitive = (a_min gt 0.0)
		WIDGET_CONTROL, sa_min, SET_VALUE = a_range[0]
		WIDGET_CONTROL, sa_max, SET_VALUE = a_range[1]
		WIDGET_CONTROL, sr_min, SET_VALUE = r_range[0]
		WIDGET_CONTROL, sr_max, SET_VALUE = r_range[1]
		pc_radial_profile_draw
	end
end


; Calculate and draw a radial profile of the given 3D data
;
; data:        3D data cube (can be including ghost cells)
; coord:       coordinate structures for the position of data in the cube,
;              asumed to be in the center of the data cube (eg. without ghost cells).
;              If omitted, the index numbers are used as coordinates.
; anchor:      anchor point for spherical shell integration.
;              If omitted, the midpoint of the box is used as anchor point.
; title:       title string for the plot
; min:         initial minimum value for data display
; max:         initial maximum value for data display
; log:         set this to use a logarithmic scale for data display
; horiz_label: label string for the axis (horizontal)
; vert_label:  label string for the value (vertical)
; file_label:  label string for filenames (special characters will be filtered)
; time:        timestamp of the displayed data
;
pro pc_radial_profile, data, coord=coord, anchor=anchor, title=title, horiz_label=horiz_label, vert_label=vert_label, min=min, max=max, log=log, file_label=file_label, time=time

	common radial_prof_common, coords, axis_dir, mid, num, t, prof_name, prof, a_range, a_min, a_max, a_label, r_range, r_min, r_max, r_label, file_name
	common radial_prof_GUI_common, win, l_plot, plot_style, b_zero, b_line, b_log, show_zero, show_line, log_plot, sa_set, sr_set, sr_fr, r_coupled, sa_max, sa_min, sr_max, sr_min

	if (n_elements (anchor) eq 0) then anchor = [ mean (coord.x), mean (coord.y), mean (coord.z) ]
	mid = anchor

	dimensionality = size (data, /n_dimensions)
	dims = size (data, /dimensions)
	nx = dims[0]
	if (dimensionality ge 2) then ny = dims[1]
	if (dimensionality ge 3) then nz = dims[2]

	; find radial direction
	axis_dir = 0
	if (max (coord.(axis_dir)) gt max (coord.y)) then axis_dir = 1
	if (max (coord.(axis_dir)) gt max (coord.z)) then axis_dir = 2
	r = coord.(axis_dir) - anchor[axis_dir]
	good = where (r ge 0.0, num)
	if (num le 0) then begin
		r = reverse (-r)
		good = where (r ge 0.0, num)
		if (num le 0) then begin
			print, "ERROR: no valid coordinates available in radial direction"
			return
		end
	end

	; center coordinates around anchor point
	coords = { x:coord.x-anchor[0], y:coord.y-anchor[1], z:coord.z-anchor[2], r:r[good] }

	t = "N/A"
	if (keyword_set (time)) then t = time

	prof_name = ""
	if (keyword_set (title)) then prof_name = title

	if (keyword_set (horiz_label)) then a_label = horiz_label else a_label = ""
	if (keyword_set (vert_label)) then r_label = vert_label else r_label = ""

	if (keyword_set (file_label)) then file_name = file_label else file_name = prof_name
	if (file_name eq "") then file_name = "radial_profile"
	pos = stregex (file_name, '[ \/]+', length=len)
	while (pos gt 0) do begin
		file_name = strmid (file_name, 0, pos) + "_" + strmid (file_name, pos+len)
		pos = stregex (file_name, '[ \/]+', length=len)
	end
	file_label = file_name

	; extend radial profile plot range by half a grid distance
	if (num le 1) then begin
		r_range = [coords.r[0]-1, coords.r[0]+1]
	end else begin
		r_range = [2*coords.r[0]-coords.r[1], 2*coords.r[num-1]-coords.r[num-2]]
	end
	r_min = r_range[0]
	r_max = r_range[1]

	; compute radial profile
	prof = dblarr (num)
	x_2 = coord.x^2
	y_2 = coord.y^2
	z_2 = coord.z^2
	cutoff = r_range[1]
	for pos_z = 0, nz - 1 do begin
		for pos_y = 0, ny - 1 do begin
			r = sqrt (x_2 + y_2[pos_y] + z_2[pos_z])
			for pos_x = 0, nx - 1 do begin
				if (r[pos_x] le cutoff) then begin
					pos = pc_find_index (r[pos_x], coords.r, /round)
					if (pos le num-1) then prof[pos] += data[pos_x,pos_y,pos_z]
				end
			end
		end
	end
	prof /= double (nx*ny*nz)

	a_range = minmax (prof)
	a_min = a_range[0]
	a_max = a_range[1]

	; GUI default settings
	win = -1
	plot_width = 600
	plot_height = 500
	sl_width = (plot_width - 100) / 2
	pc_radial_profile_reset

	if (prof_name eq "") then add = "" else add = " of "
	MOTHER	= WIDGET_BASE (title="PC radial profile analysis"+add+prof_name)
	APP	= WIDGET_BASE (MOTHER, /col)

	BASE	= WIDGET_BASE (APP, /row)

	CTRL	= WIDGET_BASE (BASE, /row)

	BUT	= WIDGET_BASE (CTRL, /col)
	tmp	= WIDGET_BUTTON (BUT, xsize=100, value='RESET', uvalue='RESET')
	tmp	= WIDGET_BUTTON (BUT, xsize=100, value='SAVE PLOT', uvalue='IMAGE')
	tmp	= WIDGET_BUTTON (BUT, xsize=100, value='QUIT', uvalue='QUIT')

	BUT	= WIDGET_BASE (CTRL, /col, frame=1, /align_left)
	l_plot	= WIDGET_DROPLIST (BUT, value=['line', 'line+dots', 'line+stars', 'dots', 'stars'], uvalue='STYLE', title='plot mean as:')

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
	tmp	= WIDGET_LABEL (PLOTS, value='radial axis:')
	BUT	= WIDGET_BASE (PLOTS, /row)
	sa_min	= CW_FSLIDER (BUT, xsize=sl_width, title='minimum value', uvalue='A_MIN', /double, /edit, min=range[0], max=range[1], drag=1, value=a_range[0])
	CTRL	= WIDGET_BASE (BUT, /col, frame=0)
	sa_set	= WIDGET_BUTTON (CTRL, value='<= MINMAX =>', uvalue='A_MINMAX', sensitive=0)
	tmp	= WIDGET_BUTTON (CTRL, value='<= CENTER =>', uvalue='A_CENTER', sensitive=(a_min lt 0.0))
	sa_max	= CW_FSLIDER (BUT, xsize=sl_width, title='maximum value', uvalue='A_MAX', /double, /edit, min=range[0], max=range[1], drag=1, value=a_range[1])

	range = get_val_range (r_range)
	tmp	= WIDGET_LABEL (PLOTS, value='vertical axis:')
	BUT	= WIDGET_BASE (PLOTS, /row)
	sr_min	= CW_FSLIDER (BUT, xsize=sl_width, title='minimum value', uvalue='R_MIN', /double, /edit, min=range[0], max=range[1], drag=1, value=r_range[0])
	CTRL	= WIDGET_BASE (BUT, /col, frame=0)
	sr_set	= WIDGET_BUTTON (CTRL, value='<= MINMAX =>', uvalue='R_MINMAX', sensitive=0)
	sr_fr	= WIDGET_BUTTON (CTRL, value='<= COUPLE =>', uvalue='R_COUPLE', sensitive=0)
	sr_max	= CW_FSLIDER (BUT, xsize=sl_width, title='maximum value', uvalue='R_MAX', /double, /edit, min=range[0], max=range[1], drag=1, value=r_range[1])

	BASE	= WIDGET_BASE (APP, /row)


	WIDGET_CONTROL, MOTHER, /REALIZE
	WIDGET_CONTROL, d_prof, GET_VALUE = win

	WIDGET_CONTROL, BASE

	pc_radial_profile_reset

	if (keyword_set (log)) then log_plot = (log ne 0) and (a_min gt 0.0)
	if (keyword_set (min)) then a_range[0] = min
	if (keyword_set (max)) then a_range[1] = max
	WIDGET_CONTROL, b_log, SET_VALUE = log_plot
	WIDGET_CONTROL, b_log, sensitive = (a_min gt 0.0)
	WIDGET_CONTROL, sa_min, SET_VALUE = a_range[0]
	WIDGET_CONTROL, sa_max, SET_VALUE = a_range[1]

	XMANAGER, "pc_radial_profile", MOTHER, /no_block

	pc_radial_profile_draw

end

