;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_show_ts.pro     ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  $Id$
;;;
;;;  Description:
;;;   GUI for investigation and analysis of the timeseries.
;;;
;;;  To do:
;;;   Add more comments


; Event handling of visualisation window
pro pc_show_ts_event, event

	common timeseries_common, time_start, time_end, ts, units, run_par, start_par, orig_dim, lvx_min, lvx_max, lvy_min, lvy_max, rvx_min, rvx_max, rvy_min, rvy_max, l_plot, r_plot, l_xy, r_xy, l_sx, l_sy, r_sx, r_sy, plot_style
	common timeseries_gui_common, l_x, l_y, r_x, r_y, ls_min, ls_max, rs_min, rs_max, ls_fr, rs_fr, ls_xy, rs_xy, l_coupled, r_coupled, lx_range, ly_range, rx_range, ry_range, s_line

	WIDGET_CONTROL, WIDGET_INFO (event.top, /CHILD)
	WIDGET_CONTROL, event.id, GET_UVALUE = eventval

	quit = -1
	L_DRAW_TS = 0
	R_DRAW_TS = 0

	SWITCH eventval of
	'ANALYZE': begin
		pc_show_ts_analyze
		break
	end
	'L_XY': begin
		if (l_xy ne event.index) then begin
			l_xy = event.index
			if (l_xy eq 0) then val_min=[lvx_min,lx_range] else val_min=[lvy_min,ly_range]
			if (l_xy eq 0) then val_max=[lvx_max,lx_range] else val_max=[lvy_max,ly_range]
			WIDGET_CONTROL, ls_min, SET_VALUE = val_min
			WIDGET_CONTROL, ls_max, SET_VALUE = val_max
			if (l_coupled ge 0) then if (l_xy eq 0) then l_coupled = lvx_max-lvx_min else l_coupled = lvy_max-lvy_min
		end
		break
	end
	'R_XY': begin
		if (r_xy ne event.index) then begin
			r_xy = event.index
			if (r_xy eq 0) then val_min=[rvx_min,rx_range] else val_min=[rvy_min,ry_range]
			if (r_xy eq 0) then val_max=[rvx_max,rx_range] else val_max=[rvy_max,ry_range]
			WIDGET_CONTROL, rs_min, SET_VALUE = val_min
			WIDGET_CONTROL, rs_max, SET_VALUE = val_max
			if (r_coupled ge 0) then if (r_xy eq 0) then r_coupled = rvx_max-rvx_min else r_coupled = rvy_max-rvy_min
		end
		break
	end
	'LS_MIN': begin
		WIDGET_CONTROL, ls_min, GET_VALUE = val_min
		if (l_coupled lt 0) then begin
			if (l_xy eq 0) then val_max = lvx_max else val_max = lvy_max
			if (val_min gt val_max) then begin
				val_min = val_max
				WIDGET_CONTROL, ls_min, SET_VALUE = val_min
			end
		end else begin
			if (l_xy eq 0) then val_max = lx_range[1] else val_max = ly_range[1]
			if (val_min gt val_max-l_coupled) then begin
				val_min = val_max-l_coupled
				WIDGET_CONTROL, ls_min, SET_VALUE = val_min
			end
			val_max = val_min+l_coupled
			WIDGET_CONTROL, ls_max, SET_VALUE = val_max
			if (l_xy eq 0) then lvx_max = val_max else lvy_max = val_max
		end
		if (l_xy eq 0) then lvx_min = val_min else lvy_min = val_min
		L_DRAW_TS = 1
		break
	end
	'LS_MAX': begin
		WIDGET_CONTROL, ls_max, GET_VALUE = val_max
		if (l_coupled lt 0) then begin
			if (l_xy eq 0) then val_min = lvx_min else val_min = lvy_min
			if (val_max lt val_min) then begin
				val_max = val_min
				WIDGET_CONTROL, ls_max, SET_VALUE = val_max
			end
		end else begin
			if (l_xy eq 0) then val_min = lx_range[0] else val_min = ly_range[0]
			if (val_max lt val_min+l_coupled) then begin
				val_max = val_min+l_coupled
				WIDGET_CONTROL, ls_max, SET_VALUE = val_max
			end
			val_min = val_max-l_coupled
			WIDGET_CONTROL, ls_min, SET_VALUE = val_min
			if (l_xy eq 0) then lvx_min = val_min else lvy_min = val_min
		end
		if (l_xy eq 0) then lvx_max = val_max else lvy_max = val_max
		L_DRAW_TS = 1
		break
	end
	'RS_MIN': begin
		WIDGET_CONTROL, rs_min, GET_VALUE = val_min
		if (r_coupled lt 0) then begin
			if (r_xy eq 0) then val_max = rvx_max else val_max = rvy_max
			if (val_min gt val_max) then begin
				val_min = val_max
				WIDGET_CONTROL, rs_min, SET_VALUE = val_min
			end
		end else begin
			if (r_xy eq 0) then val_max = rx_range[1] else val_max = ry_range[1]
			if (val_min gt val_max-r_coupled) then begin
				val_min = val_max-r_coupled
				WIDGET_CONTROL, rs_min, SET_VALUE = val_min
			end
			val_max = val_min+r_coupled
			WIDGET_CONTROL, rs_max, SET_VALUE = val_max
			if (r_xy eq 0) then rvx_max = val_max else rvy_max = val_max
		end
		if (r_xy eq 0) then rvx_min = val_min else rvy_min = val_min
		R_DRAW_TS = 1
		break
	end
	'RS_MAX': begin
		WIDGET_CONTROL, rs_max, GET_VALUE = val_max
		if (r_coupled lt 0) then begin
			if (r_xy eq 0) then val_min = rvx_min else val_min = rvy_min
			if (val_max lt val_min) then begin
				val_max = val_min
				WIDGET_CONTROL, rs_max, SET_VALUE = val_max
			end
		end else begin
			if (r_xy eq 0) then val_min = rx_range[0] else val_min = ry_range[0]
			if (val_max lt val_min+r_coupled) then begin
				val_max = val_min+r_coupled
				WIDGET_CONTROL, rs_max, SET_VALUE = val_max
			end
			val_min = val_max-r_coupled
			WIDGET_CONTROL, rs_min, SET_VALUE = val_min
			if (r_xy eq 0) then rvx_min = val_min else rvy_min = val_min
		end
		if (r_xy) eq 0 then rvx_max = val_max else rvy_max = val_max
		R_DRAW_TS = 1
		break
	end
	'L_X': begin
		if (l_sx ne event.index) then begin
			l_sx = event.index
			lx_range = minmax (ts.(l_sx))
			lvx_min = lx_range[0]
			lvx_max = lx_range[1]
			if (l_xy eq 0) then begin
				WIDGET_CONTROL, ls_min, SET_VALUE = [lvx_min,lx_range]
				WIDGET_CONTROL, ls_max, SET_VALUE = [lvx_max,lx_range]
			end
			L_DRAW_TS = 1
		end
		break
	end
	'L_Y': begin
		if (l_sy ne event.index) then begin
			l_sy = event.index
			ly_range = minmax (ts.(l_sy))
			lvy_min = ly_range[0]
			lvy_max = ly_range[1]
			if (l_xy eq 1) then begin
				WIDGET_CONTROL, ls_min, SET_VALUE = [lvy_min,ly_range]
				WIDGET_CONTROL, ls_max, SET_VALUE = [lvy_max,ly_range]
			end
			L_DRAW_TS = 1
		end
		break
	end
	'R_X': begin
		if (r_sx ne event.index) then begin
			r_sx = event.index
			rx_range = minmax (ts.(r_sx))
			rvx_min = rx_range[0]
			rvx_max = rx_range[1]
			if (r_xy eq 0) then begin
				WIDGET_CONTROL, rs_min, SET_VALUE = [rvx_min,rx_range]
				WIDGET_CONTROL, rs_max, SET_VALUE = [rvx_max,rx_range]
			end
			R_DRAW_TS = 1
		end
		break
	end
	'R_Y': begin
		if (r_sy ne event.index) then begin
			r_sy = event.index
			ry_range = minmax (ts.(r_sy))
			rvy_min = ry_range[0]
			rvy_max = ry_range[1]
			if (r_xy eq 1) then begin
				WIDGET_CONTROL, rs_min, SET_VALUE = [rvy_min,ry_range]
				WIDGET_CONTROL, rs_max, SET_VALUE = [rvy_max,ry_range]
			end
			R_DRAW_TS = 1
		end
		break
	end
	'STYLE': begin
		if (plot_style ne event.index) then begin
			plot_style = event.index
			L_DRAW_TS = 1
			R_DRAW_TS = 1
		end
		break
	end
	'RESET': begin
		; reset_ts_GUI
		break
	end
	'REFRESH': begin
		pc_read_ts, obj=ts, datadir=datadir, /quiet
		lx_range = minmax (ts.(l_sx))
		ly_range = minmax (ts.(l_sy))
		rx_range = minmax (ts.(r_sx))
		ry_range = minmax (ts.(r_sy))
		if (l_xy eq 0) then begin
			WIDGET_CONTROL, ls_min, SET_VALUE = [lvx_min,lx_range]
			WIDGET_CONTROL, ls_max, SET_VALUE = [lvx_max,lx_range]
		end else begin
			WIDGET_CONTROL, ls_min, SET_VALUE = [lvy_min,ly_range]
			WIDGET_CONTROL, ls_max, SET_VALUE = [lvy_max,ly_range]
		end
		if (r_xy eq 0) then begin
			WIDGET_CONTROL, rs_min, SET_VALUE = [rvx_min,rx_range]
			WIDGET_CONTROL, rs_max, SET_VALUE = [rvx_max,rx_range]
		end else begin
			WIDGET_CONTROL, rs_min, SET_VALUE = [rvy_min,ry_range]
			WIDGET_CONTROL, rs_max, SET_VALUE = [rvy_max,ry_range]
		end
		pc_show_ts_draw, 1, 1
		break
	end
	'L_COUPLE': begin
		WIDGET_CONTROL, ls_fr, set_value='<= RELEASE =>', set_uvalue='L_RELEASE'
		if (l_xy eq 0) then l_coupled = lvx_max - lvx_min else l_coupled = lvy_max - lvy_min
		break
	end
	'L_RELEASE': begin
		WIDGET_CONTROL, ls_fr, set_value='<= COUPLE =>', set_uvalue='L_COUPLE'
		l_coupled = -1
		break
	end
	'R_COUPLE': begin
		WIDGET_CONTROL, rs_fr, set_value='<= RELEASE =>', set_uvalue='R_RELEASE'
		if (r_xy eq 0) then r_coupled = rvx_max - rvx_min else r_coupled = rvy_max - rvy_min
		break
	end
	'R_RELEASE': begin
		WIDGET_CONTROL, rs_fr, set_value='<= COUPLE =>', set_uvalue='R_COUPLE'
		r_coupled = -1
		break
	end
	'QUIT': begin
		quit = event.top
		break
	end
	endswitch

	pc_show_ts_draw, L_DRAW_TS, R_DRAW_TS

	WIDGET_CONTROL, WIDGET_INFO (event.top, /CHILD)

	if (quit ge 0) then WIDGET_CONTROL, quit, /DESTROY

	return
end


; Draw the timeseries plots
pro pc_show_ts_draw, l_draw, r_draw

	common timeseries_common, time_start, time_end, ts, units, run_par, start_par, orig_dim, lvx_min, lvx_max, lvy_min, lvy_max, rvx_min, rvx_max, rvy_min, rvy_max, l_plot, r_plot, l_xy, r_xy, l_sx, l_sy, r_sx, r_sy, plot_style

	if (l_draw ne 0) then begin
		wset, l_plot
		xr = [lvx_min,lvx_max]
		if (lvx_min eq lvx_max) then xr[1] = lvx_min * (1.0 + 1.e-14)
		if (plot_style le 2) then plot, ts.(l_sx), ts.(l_sy), xr=xr, yr=[lvy_min*0.95,lvy_max*1.05], /xs, /ys
		if (plot_style eq 1) then oplot, ts.(l_sx), ts.(l_sy), psym=3, color=200
		if (plot_style eq 2) then oplot, ts.(l_sx), ts.(l_sy), psym=2, color=200
		if (plot_style eq 3) then plot, ts.(l_sx), ts.(l_sy), xr=xr, yr=[lvy_min*0.95,lvy_max*1.05], /xs, /ys, psym=3
		if (plot_style eq 4) then plot, ts.(l_sx), ts.(l_sy), xr=xr, yr=[lvy_min*0.95,lvy_max*1.05], /xs, /ys, psym=2
	end

	if (r_draw ne 0) then begin
		wset, r_plot
		xr = [rvx_min,rvx_max]
		if (rvx_min eq rvx_max) then xr[1] = rvx_min * (1.0 + 1.e-14)
		if (plot_style le 2) then plot, ts.(r_sx), ts.(r_sy), xr=xr, yr=[rvy_min*0.95,rvy_max*1.05], /xs, /ys
		if (plot_style eq 1) then oplot, ts.(r_sx), ts.(r_sy), psym=3, color=200
		if (plot_style eq 2) then oplot, ts.(r_sx), ts.(r_sy), psym=2, color=200
		if (plot_style eq 3) then plot, ts.(r_sx), ts.(r_sy), xr=xr, yr=[rvy_min*0.95,rvy_max*1.05], /xs, /ys, psym=3
		if (plot_style eq 4) then plot, ts.(r_sx), ts.(r_sy), xr=xr, yr=[rvy_min*0.95,rvy_max*1.05], /xs, /ys, psym=2
	end
end


; Analyze the timeseries plots
pro pc_show_ts_analyze

	common timeseries_common, time_start, time_end, ts, units, run_par, start_par, orig_dim, lvx_min, lvx_max, lvy_min, lvy_max, rvx_min, rvx_max, rvy_min, rvy_max, l_plot, r_plot, l_xy, r_xy, l_sx, l_sy, r_sx, r_sy, plot_style

	charsize = 1.25
	old_multi = !P.MULTI
	old_x_margin = !X.margin
	!X.margin[0] += 3
	x_margin_both = (!X.margin > max (old_x_margin))

	window, 11, xsize=1000, ysize=400, title='timestep analysis', retain=2
	!P.MULTI = [0, 2, 1]

	print, "starting values:"
	print, "dt    :", ts.dt[0]
	plot, ts.dt, title = 'dt', xc=charsize, yc=charsize, /yl

	tags = tag_names (ts)
	if (any (strcmp (tags, 't', /fold_case))) then begin
		time = ts.t
	endif else begin
		time = ts.it
	endelse
	x_minmax = minmax (time > time_start)
	if (time_end gt 0) then x_minmax = minmax (x_minmax < time_end)
	y_minmax = minmax (ts.dt)
	if (any (strcmp (tags, 'dtu', /fold_case)))       then y_minmax = minmax ([y_minmax, ts.dtu])
	if (any (strcmp (tags, 'dtv', /fold_case)))       then y_minmax = minmax ([y_minmax, ts.dtv])
	if (any (strcmp (tags, 'dtnu', /fold_case)))      then y_minmax = minmax ([y_minmax, ts.dtnu])
	if (any (strcmp (tags, 'dtb', /fold_case)))       then y_minmax = minmax ([y_minmax, ts.dtb])
	if (any (strcmp (tags, 'dteta', /fold_case)))     then y_minmax = minmax ([y_minmax, ts.dteta])
	if (any (strcmp (tags, 'dtc', /fold_case)))       then y_minmax = minmax ([y_minmax, ts.dtc])
	if (any (strcmp (tags, 'dtchi', /fold_case)))     then y_minmax = minmax ([y_minmax, ts.dtchi])
	if (any (strcmp (tags, 'dtchi2', /fold_case)))    then y_minmax = minmax ([y_minmax, ts.dtchi2])
	if (any (strcmp (tags, 'dtspitzer', /fold_case))) then y_minmax = minmax ([y_minmax, ts.dtspitzer])
	if (any (strcmp (tags, 'dtd', /fold_case)))       then y_minmax = minmax ([y_minmax, ts.dtd])

	time *= units.time
	ts.dt *= units.time
	x_minmax *= units.time
	y_minmax *= units.time

	plot, time, ts.dt, title = 'dt(t) u{-t} v{-p} nu{.v} b{.r} eta{-g} c{.y} chi{-.b} chi2{-.o} d{-l} [s]', xrange=x_minmax, /xs, xc=charsize, yc=charsize, yrange=y_minmax, /yl
	if (any (strcmp (tags, 'dtu', /fold_case))) then begin
		oplot, time, ts.dtu*units.time, linestyle=2, color=11061000
		print, "dtu      :", ts.dtu[0]
	end
	if (any (strcmp (tags, 'dtv', /fold_case))) then begin
		oplot, time, ts.dtv*units.time, linestyle=2, color=128255200
		print, "dtv      :", ts.dtv[0]
	end
	if (any (strcmp (tags, 'dtnu', /fold_case))) then begin
		oplot, time, ts.dtnu*units.time, linestyle=1, color=128000128
		print, "dtnu     :", ts.dtnu[0]
	end
	if (any (strcmp (tags, 'dtb', /fold_case))) then begin
		oplot, time, ts.dtb*units.time, linestyle=1, color=200
		print, "dtb      :", ts.dtb[0]
	end
	if (any (strcmp (tags, 'dteta', /fold_case))) then begin
		oplot, time, ts.dteta*units.time, linestyle=2, color=220200200
		print, "dteta    :", ts.dteta[0]
	end
	if (any (strcmp (tags, 'dtc', /fold_case))) then begin
		oplot, time, ts.dtc*units.time, linestyle=1, color=61695
		print, "dtc      :", ts.dtc[0]
	end
	if (any (strcmp (tags, 'dtchi', /fold_case))) then begin
		oplot, time, ts.dtchi*units.time, linestyle=3, color=115100200
		print, "dtchi    :", ts.dtchi[0]
	end
	if (any (strcmp (tags, 'dtchi2', /fold_case))) then begin
		oplot, time, ts.dtchi2*units.time, linestyle=3, color=41215
		print, "dtchi2   :", ts.dtchi2[0]
	end
	if (any (strcmp (tags, 'dtspitzer', /fold_case))) then begin
		oplot, time, ts.dtspitzer*units.time, linestyle=3, color=41215000
		print, "dtspitzer:", ts.dtspitzer[0]
	end
	if (any (strcmp (tags, 'dtd', /fold_case))) then begin
		oplot, time, ts.dtd*units.time, linestyle=2, color=16737000
		print, "dtc      :", ts.dtd[0]
	end

	window, 12, xsize=1000, ysize=800, title='time series analysis', retain=2
	!P.MULTI = [0, 2, 2, 0, 0]

	max_subplots = 4
	num_subplots = 0

	N_grid = orig_dim.nx * orig_dim.ny * orig_dim.nz
	volume = start_par.Lxyz[0] * start_par.Lxyz[1] * start_par.Lxyz[2]

	if (any (strcmp (tags, 'eem', /fold_case)) and any (strcmp (tags, 'ekintot', /fold_case)) and any (strcmp (tags, 'ethm', /fold_case)) and (num_subplots lt max_subplots)) then begin
		num_subplots += 1
		energy_int = (ts.eem*N_grid + ts.ekintot) * units.mass / units.velocity^2
		energy_therm = (ts.ethm + ts.ekintot/volume) * units.mass / (units.velocity^2 * units.length^3)
		plot, time, energy_int, title = 'Energy {w} and E/V {r} conservation', xrange=x_minmax, /xs, xmar=x_margin_both, xc=charsize, yc=charsize, ytitle='<E> [J]', ys=10, /noerase
		plot, time, energy_therm, color=200, xrange=x_minmax, xs=5, xmar=x_margin_both, xc=charsize, yc=charsize, ys=6, /noerase
		axis, xc=charsize, yc=charsize, yaxis=1, yrange=!Y.CRANGE, /ys, ytitle='<E/V> [J/kg]'
		plot, time, energy_int, linestyle=2, xrange=x_minmax, xs=5, xmar=x_margin_both, xc=charsize, yc=charsize, ys=6, /noerase
		!P.MULTI = [max_subplots-num_subplots, 2, 2, 0, 0]
	end else if (any (strcmp (tags, 'eem', /fold_case)) and any (strcmp (tags, 'ekintot', /fold_case)) and any (strcmp (tags, 'totmass', /fold_case)) and (num_subplots lt max_subplots)) then begin
		num_subplots += 1
		mass = ts.totmass * units.mass / units.default_mass
		energy = (ts.eem*N_grid + ts.ekintot) / ts.totmass / units.velocity^2
		plot, time, energy, title = 'Energy {w} and mass {r} conservation', xrange=x_minmax, /xs, xmar=x_margin_both, xc=charsize, yc=charsize, ytitle='<E>_m [J/kg]', ys=10, /noerase
		plot, time, mass, color=200, xrange=x_minmax, xs=5, xmar=x_margin_both, xc=charsize, yc=charsize, ys=6, /noerase
		axis, xc=charsize, yc=charsize, yaxis=1, yrange=!Y.CRANGE, /ys, ytitle='total mass ['+units.default_mass_str+']'
		plot, time, energy, linestyle=2, xrange=x_minmax, xs=5, xmar=x_margin_both, xc=charsize, yc=charsize, ys=6, /noerase
		!P.MULTI = [max_subplots-num_subplots, 2, 2, 0, 0]
	end else if (any (strcmp (tags, 'eem', /fold_case)) and any (strcmp (tags, 'ekintot', /fold_case)) and (num_subplots lt max_subplots)) then begin
		num_subplots += 1
		energy = (ts.eem*N + ts.ekintot) * units.mass / units.velocity^2
		plot, time, energy, title = 'Energy conservation', xrange=x_minmax, /xs, xc=charsize, yc=charsize, ytitle='<E> [J]'
	end else if (any (strcmp (tags, 'ethm', /fold_case)) and any (strcmp (tags, 'ekintot', /fold_case)) and (num_subplots lt max_subplots)) then begin
		num_subplots += 1
		energy = (ts.ethm + ts.ekintot/volume) * units.mass / (units.velocity^2 * units.length^3)
		plot, time, energy, title = 'Energy conservation', xrange=x_minmax, /xs, xc=charsize, yc=charsize, ytitle='<E/V> [J/m^3]'
	end
	if (any (strcmp (tags, 'TTmax', /fold_case)) and any (strcmp (tags, 'rhomin', /fold_case)) and (num_subplots lt max_subplots)) then begin
		num_subplots += 1
		Temp_max = ts.TTmax * units.temperature
		rho_min = ts.rhomin * units.density / units.default_density
		plot, time, Temp_max, title = 'Maximum temperature {w} and minimum density {.r}', xrange=x_minmax, /xs, xmar=x_margin_both, xc=charsize, yc=charsize, ytitle='maximum temperature [K]', /yl, ys=10, /noerase
		plot, time, rho_min, color=200, xrange=x_minmax, xs=5, xmar=x_margin_both, xc=charsize, yc=charsize, /yl, ys=6, /noerase
		axis, xc=charsize, yc=charsize, yaxis=1, yrange=10.^(!Y.CRANGE), /ys, /yl, ytitle='minimum density ['+units.default_density_str+']'
		plot, time, Temp_max, linestyle=2, xrange=x_minmax, xs=5, xmar=x_margin_both, xc=charsize, yc=charsize, /yl, ys=6, /noerase
		!P.MULTI = [max_subplots-num_subplots, 2, 2, 0, 0]
	end else if (any (strcmp (tags, 'TTm', /fold_case)) and any (strcmp (tags, 'rhomin', /fold_case)) and (num_subplots lt max_subplots)) then begin
		num_subplots += 1
		Temp_mean = ts.TTm * units.temperature
		rho_min = ts.rhomin * units.density / units.default_density
		plot, time, Temp_mean, title = 'Mean temperature {w} and minimum density {.r}', xrange=x_minmax, /xs, xmar=x_margin_both, xc=charsize, yc=charsize, ytitle='<T> [K]', /yl, ys=10, /noerase
		plot, time, rho_min, color=200, xrange=x_minmax, xs=5, xmar=x_margin_both, xc=charsize, yc=charsize, /yl, ys=6, /noerase
		axis, xc=charsize, yc=charsize, yaxis=1, yrange=10.^(!Y.CRANGE), /ys, /yl, ytitle='minimum density ['+units.default_density_str+']'
		plot, time, Temp_mean, linestyle=2, xrange=x_minmax, xs=5, xmar=x_margin_both, xc=charsize, yc=charsize, /yl, ys=6, /noerase
		!P.MULTI = [max_subplots-num_subplots, 2, 2, 0, 0]
	end else if (any (strcmp (tags, 'TTm', /fold_case)) and any (strcmp (tags, 'TTmax', /fold_case)) and (num_subplots lt max_subplots)) then begin
		num_subplots += 1
		Temp_max = ts.TTmax * units.temperature
		Temp_mean = ts.TTm * units.temperature
		yrange = [ min (Temp_mean), max (Temp_max) ]
		plot, time, Temp_max, title = 'Maximum temperature {w} and mean temperature {.r}', xrange=x_minmax, /xs, xc=charsize, yc=charsize, ytitle='maximum and mean temperature [K]', yrange=yrange, /yl
		oplot, time, Temp_mean, color=200
		oplot, time, Temp_max, linestyle=2
	end else if (any (strcmp (tags, 'TTmax', /fold_case)) and (num_subplots lt max_subplots)) then begin
		num_subplots += 1
		Temp_max = ts.TTmax * units.temperature
		plot, time, Temp_max, title = 'Maximum temperature [K]', xrange=x_minmax, /xs, xc=charsize, yc=charsize, /yl
	end else if (any (strcmp (tags, 'TTm', /fold_case)) and (num_subplots lt max_subplots)) then begin
		num_subplots += 1
		Temp_mean = ts.TTm * units.temperature
		plot, time, Temp_mean, title = 'Mean temperature [K]', xrange=x_minmax, /xs, xc=charsize, yc=charsize, /yl
	end else if (any (strcmp (tags, 'rhomin', /fold_case)) and (num_subplots lt max_subplots)) then begin
		num_subplots += 1
		rho_min = ts.rhomin * units.density / units.default_density
		plot, time, rho_min, title = 'rho_min(t) ['+units.default_density_str+']', xrange=x_minmax, /xs, xc=charsize, yc=charsize, /yl
	end
	if (any (strcmp (tags, 'j2m', /fold_case)) and any (strcmp (tags, 'visc_heatm', /fold_case)) and any (tag_names (run_par) eq "ETA") and (num_subplots lt max_subplots)) then begin
		num_subplots += 1
		HR_ohm = run_par.eta * start_par.mu0 * ts.j2m * units.density * units.velocity^3 / units.length
		visc_heat_mean = ts.visc_heatm * units.density * units.velocity^3 / units.length
		yrange = [ min ([HR_ohm, visc_heat_mean]), max ([HR_ohm, visc_heat_mean]) ]
		plot, time, HR_ohm, title = 'Mean Ohmic heating rate {w} and viscous heating rate {.r}', xrange=x_minmax, /xs, xc=charsize, yc=charsize, ytitle='heating rates [W/m^3]', yrange=yrange, /yl
		oplot, time, visc_heat_mean, color=200
		oplot, time, HR_ohm, linestyle=2
	end else if (any (strcmp (tags, 'j2m', /fold_case)) and any (tag_names (run_par) eq "ETA") and (num_subplots lt max_subplots)) then begin
		num_subplots += 1
		mu0_SI = 4.0 * !Pi * 1.e-7
		HR_ohm = run_par.eta * start_par.mu0 * ts.j2m * units.density * units.velocity^3 / units.length
		j_abs = sqrt (ts.j2m) * units.velocity * sqrt (start_par.mu0 / mu0_SI * units.density) / units.length
		plot, time, HR_ohm, title = 'Mean Ohmic heating rate {w} and mean current density {.r}', xrange=x_minmax, /xs, xmar=x_margin_both, xc=charsize, yc=charsize, ytitle='HR = <eta*mu0*j^2> [W/m^3]', /yl, ys=10, /noerase
		plot, time, j_abs, color=200, xrange=x_minmax, xs=5, xmar=x_margin_both, xc=charsize, yc=charsize, /yl, ys=6, /noerase
		axis, xc=charsize, yc=charsize, yaxis=1, yrange=10.^(!Y.CRANGE), /ys, /yl, ytitle='sqrt(<j^2>) [A/m^2]'
		plot, time, HR_ohm, linestyle=2, xrange=x_minmax, xs=5, xmar=x_margin_both, xc=charsize, yc=charsize, /yl, ys=6, /noerase
		!P.MULTI = [max_subplots-num_subplots, 2, 2, 0, 0]
	end else if (any (strcmp (tags, 'visc_heatm', /fold_case)) and (num_subplots lt max_subplots)) then begin
		num_subplots += 1
		visc_heat_mean = ts.visc_heatm * units.density * units.velocity^3 / units.length
		plot, time, visc_heat_mean, title = 'Mean viscous heating rate', xrange=x_minmax, /xs, xc=charsize, yc=charsize, ytitle='heating rate [W/m^3]', /yl
	end
	if (any (strcmp (tags, 'umax', /fold_case)) and (num_subplots lt max_subplots)) then begin
		num_subplots += 1
		u_max = ts.umax * units.velocity / units.default_velocity
		u_title = 'u_max(t){w}'
		if (any (strcmp (tags, 'urms', /fold_case))) then begin
			u_title += ' u_rms{.r}'
		end else if (any (strcmp (tags, 'u2m', /fold_case))) then begin
			u_title += ' sqrt(<u^2>){.-b}'
		end
		plot, time, u_max, title = u_title+' ['+units.default_velocity_str+']', xrange=x_minmax, /xs, xc=charsize, yc=charsize
		if (any (strcmp (tags, 'urms', /fold_case))) then begin
			urms = ts.urms * units.velocity / units.default_velocity
			oplot, time, urms, linestyle=1, color=200
		end else if (any (strcmp (tags, 'u2m', /fold_case))) then begin
			u2m = sqrt (ts.u2m) * units.velocity / units.default_velocity
			oplot, time, u2m, linestyle=3, color=115100200
		end
	end
	if (any (strcmp (tags, 'totmass', /fold_case)) and (num_subplots lt max_subplots)) then begin
		num_subplots += 1
		mass = ts.totmass * units.mass / units.default_mass
		plot, time, mass, title = 'Mass conservation', xrange=x_minmax, /xs, xc=charsize, yc=charsize
	end

	!X.margin = old_x_margin
	!P.MULTI = old_multi
end


; Show timeseries analysis window
pro pc_show_ts, object=time_series, unit=unit, param=param, run_param=run_param, start_time=start_time, end_time=end_time, datadir=datadir

	common timeseries_common, time_start, time_end, ts, units, run_par, start_par, orig_dim, lvx_min, lvx_max, lvy_min, lvy_max, rvx_min, rvx_max, rvy_min, rvy_max, l_plot, r_plot, l_xy, r_xy, l_sx, l_sy, r_sx, r_sy, plot_style
	common timeseries_gui_common, l_x, l_y, r_x, r_y, ls_min, ls_max, rs_min, rs_max, ls_fr, rs_fr, ls_xy, rs_xy, l_coupled, r_coupled, lx_range, ly_range, rx_range, ry_range

	; GUI settings
	@pc_gui_settings
	col_width = 220
	plot_width = 2 * col_width
	plot_height = plot_width
	sl_width = col_width - 52

	if (not keyword_set (datadir)) then datadir = pc_get_datadir()
	pc_read_dim, obj=orig_dim, datadir=datadir, /quiet
	if (not keyword_set (unit)) then pc_units, obj=unit, datadir=datadir, param=param, dim=orig_dim, /quiet
	if (not any (strcmp (tag_names (unit), "default_length", /fold_case))) then unit = create_struct (unit, display_units)
	units = unit

	if (keyword_set (time_series)) then ts = time_series
	if (n_elements (ts) le 0) then pc_read_ts, obj=ts, datadir=datadir, /quiet
	time_series = ts

	if (not keyword_set (param)) then pc_read_param, obj=param, datadir=datadir, /quiet
	if (not keyword_set (run_param)) then pc_read_param, obj=run_param, datadir=datadir, /param2, /quiet
	if (not keyword_set (start_time)) then start_time = min (ts.t) * unit.time
	if (not keyword_set (end_time)) then end_time = max (ts.t) * unit.time

	time_start = start_time
	time_end = end_time
	run_par = run_param
	start_par = param

	plots = tag_names (ts)
	num_plots = n_elements (plots)

	l_sx = 0
	l_sy = 2 < (num_plots-1)
	r_sx = 1 < (num_plots-1)
	r_sy = 2 < (num_plots-1)
	l_xy = 0
	r_xy = 0
	lx_range = minmax (ts.(l_sx))
	ly_range = minmax (ts.(l_sy))
	rx_range = minmax (ts.(r_sx))
	ry_range = minmax (ts.(r_sy))
	lvx_min = lx_range[0]
	lvx_max = lx_range[1]
	lvy_min = ly_range[0]
	lvy_max = ly_range[1]
	rvx_min = rx_range[0]
	rvx_max = rx_range[1]
	rvy_min = ry_range[0]
	rvy_max = ry_range[1]
	l_coupled = -1
	r_coupled = -1
	plot_style = 0

	MOTHER	= WIDGET_BASE (title='PC timeseries analysis')
	APP	= WIDGET_BASE (MOTHER, /col)

	BASE	= WIDGET_BASE (APP, /row)

	tmp	= WIDGET_BASE (BASE, /row)
	BUT	= WIDGET_BASE (tmp, /col)
	L_X	= WIDGET_DROPLIST (BUT, xsize=plot_width-60, value=plots, uvalue='L_X', title='LEFT plot:')
	L_Y	= WIDGET_LIST (BUT, value=plots, uvalue='L_Y', ysize=(num_plots<12)>4); , /multiple
	WIDGET_CONTROL, L_X, SET_DROPLIST_SELECT = l_sx
	WIDGET_CONTROL, L_Y, SET_LIST_SELECT = l_sy

	CTRL	= WIDGET_BASE (BASE, /col)
	BUT	= WIDGET_BASE (CTRL, /col, frame=1, /align_center)
	tmp	= WIDGET_BUTTON (BUT, xsize=100, value='RESET', uvalue='RESET', sensitive=0)
	tmp	= WIDGET_BUTTON (BUT, xsize=100, value='REFRESH', uvalue='REFRESH')
	tmp	= WIDGET_BUTTON (BUT, xsize=100, value='ANALYZE', uvalue='ANALYZE')
	tmp	= WIDGET_BUTTON (BUT, xsize=100, value='QUIT', uvalue='QUIT')
	BUT	= WIDGET_BASE (CTRL, /col, /align_center)
	tmp	= WIDGET_LABEL (CTRL, value='plotting style:', frame=0)
	tmp	= WIDGET_DROPLIST (CTRL, value=['line', 'line+dots', 'line+stars', 'dots', 'stars'], uvalue='STYLE')
	WIDGET_CONTROL, tmp, SET_DROPLIST_SELECT = plot_style

	tmp	= WIDGET_BASE (BASE, /row)
	BUT	= WIDGET_BASE (tmp, /col)
	R_X	= WIDGET_DROPLIST (BUT, xsize=plot_width-60, value=plots, uvalue='R_X', title='RIGHT plot:')
	R_Y	= WIDGET_LIST (BUT, value=plots, uvalue='R_Y', ysize=(num_plots<12)>4) ; , /multiple
	WIDGET_CONTROL, R_X, SET_DROPLIST_SELECT = r_sx
	WIDGET_CONTROL, R_Y, SET_LIST_SELECT = r_sy

	BASE	= WIDGET_BASE (APP, /row)

	tmp	= WIDGET_BASE (BASE, /col)
	PLOTS	= WIDGET_BASE (tmp, /row)
	dplot_l	= WIDGET_DRAW (PLOTS, xsize=plot_width, ysize=plot_height, retain=2)

	tmp	= WIDGET_BASE (BASE, /col)
	PLOTS	= WIDGET_BASE (tmp, /row)
	dplot_r	= WIDGET_DRAW (PLOTS, xsize=plot_width, ysize=plot_height, retain=2)

	BASE	= WIDGET_BASE (APP, /row)

	BUT	= WIDGET_BASE (BASE, /row)
	ls_min	= CW_FSLIDER (BUT, xsize=sl_width, title='minimum value', uvalue='LS_MIN', /double, /edit, min=lx_range[0], max=lx_range[1], drag=1, value=lvx_min)
	CTRL	= WIDGET_BASE (BUT, /col, frame=0)
	tmp	= WIDGET_DROPLIST (CTRL, value=['X', 'Y'], uvalue='L_XY', title='axis:')
	ls_fr	= WIDGET_BUTTON (CTRL, value='<= COUPLE =>', uvalue='L_COUPLE')
	ls_max	= CW_FSLIDER (BUT, xsize=sl_width, title='maximum value', uvalue='LS_MAX', /double, /edit, min=lx_range[0], max=lx_range[1], drag=1, value=lvx_max)

	BUT	= WIDGET_BASE (BASE, /row)
	rs_min	= CW_FSLIDER (BUT, xsize=sl_width, title='minimum value', uvalue='RS_MIN', /double, /edit, min=rx_range[0], max=rx_range[1], drag=1, value=rvx_min)
	CTRL	= WIDGET_BASE (BUT, /col)
	tmp	= WIDGET_DROPLIST (CTRL, value=['X', 'Y'], uvalue='R_XY', title='axis:')
	rs_fr	= WIDGET_BUTTON (CTRL, value='<= COUPLE =>', uvalue='R_COUPLE')
	rs_max	= CW_FSLIDER (BUT, xsize=sl_width, title='maximum value', uvalue='RS_MAX', /double, /edit, min=rx_range[0], max=rx_range[1], drag=1, value=rvx_max)


	WIDGET_CONTROL, MOTHER, /REALIZE
	WIDGET_CONTROL, dplot_l, GET_VALUE = l_plot
	WIDGET_CONTROL, dplot_r, GET_VALUE = r_plot

	WIDGET_CONTROL, BASE

	XMANAGER, "pc_show_ts", MOTHER, /no_block

	pc_show_ts_draw, 1, 1

end

