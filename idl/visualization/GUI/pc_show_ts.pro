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
pro timeseries_event, event

	common timeseries_common, time_start, time_end, ts, units, run_par, start_par
;	common timeseries_gui_common, 

	WIDGET_CONTROL, WIDGET_INFO (event.top, /CHILD)
	WIDGET_CONTROL, event.id, GET_UVALUE = eventval

	quit = -1
	DRAW_TS = 0

	SWITCH eventval of
	'SCALE': begin
		abs_scale = event.select
		break
	end
	'SL_MIN': begin
		WIDGET_CONTROL, sl_min, GET_VALUE = val_min
		if (val_min gt val_max) then begin
			val_min = val_max
			WIDGET_CONTROL, sl_min, SET_VALUE = val_min
		end
		pos_b[selected_cube,sub_aver] = val_min
		DRAW_TS = 1
		break
	end
	'SL_MAX': begin
		WIDGET_CONTROL, sl_max, GET_VALUE = val_max
		if (val_max lt val_min) then begin
			val_max = val_min
			WIDGET_CONTROL, sl_max, SET_VALUE = val_max
		end
		pos_t[selected_cube,sub_aver] = val_max
		DRAW_TS = 1
		break
	end
	'VAR_LEFT': begin
		if (selected_cube ne event.index) then begin
			prepare_cube, event.index
			DRAW_TS = 1
		end
		break
	end
	'VAR_RIGHT': begin
		if (selected_cube ne event.index) then begin
			prepare_cube, event.index
			DRAW_TS = 1
		end
		break
	end
	'RESET': begin
		reset_GUI
		break
	end
	'IMAGE': begin
		WIDGET_CONTROL, image, SENSITIVE = 0
		WIDGET_CONTROL, sl_min, SENSITIVE = 0
		WIDGET_CONTROL, sl_max, SENSITIVE = 0
		save_images, "PNG"
		WIDGET_CONTROL, image, SENSITIVE = 1
		WIDGET_CONTROL, sl_min, SENSITIVE = 1
		WIDGET_CONTROL, sl_max, SENSITIVE = 1
		break
	end
	'QUIT': begin
		quit = event.top
		break
	end
	endswitch

	if (DRAW_TS) then draw_timeseries

	WIDGET_CONTROL, WIDGET_INFO (event.top, /CHILD)

	if (quit ge 0) then WIDGET_CONTROL, quit, /DESTROY

	return
end


; Draw the timeseries plots
pro draw_timeseries

	common timeseries_common, time_start, time_end, ts, units, run_par, start_par
;	common timeseries_gui_common, 

	charsize = 1.25
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

	if (any (strcmp (tags, 'eem', /fold_case)) and any (strcmp (tags, 'ethm', /fold_case)) and any (strcmp (tags, 'ekintot', /fold_case)) and any (strcmp (tags, 'totmass', /fold_case)) and (num_subplots lt max_subplots)) then begin
		num_subplots += 1
		mass = ts.totmass * units.mass / units.default_mass
		energy = (ts.eem + ts.ekintot/ts.totmass) * units.mass / units.velocity^2
		plot, time, energy, title = 'Energy {w} and mass {r} conservation', xrange=x_minmax, /xs, xmar=x_margin_both, xc=charsize, yc=charsize, ytitle='<E> [J]', ys=10, /noerase
		plot, time, mass, color=200, xrange=x_minmax, xs=5, xmar=x_margin_both, xc=charsize, yc=charsize, ys=6, /noerase
		axis, xc=charsize, yc=charsize, yaxis=1, yrange=!Y.CRANGE, /ys, ytitle='total mass ['+units.default_mass_str+']'
		plot, time, energy, linestyle=2, xrange=x_minmax, xs=5, xmar=x_margin_both, xc=charsize, yc=charsize, ys=6, /noerase
		!P.MULTI = [max_subplots-num_subplots, 2, 2, 0, 0]
	end else if (any (strcmp (tags, 'totmass', /fold_case)) and (num_subplots lt max_subplots)) then begin
		num_subplots += 1
		mass = ts.totmass * units.mass / units.default_mass
		plot, time, mass, title = 'Mass conservation', xrange=x_minmax, /xs, xc=charsize, yc=charsize
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
	!X.margin = old_x_margin
end


; Show timeseries analysis window
pro pc_show_ts, object=time_series, units=units_struct, param=param, run_param=run_param, start_time=start_time, end_time=end_time, datadir=datadir

	common timeseries_common, time_start, time_end, ts, units, run_par, start_par

	if (not keyword_set (datadir)) then datadir = pc_get_datadir()

	if (keyword_set (units_struct)) then units = units_struct
	if (n_elements (units) le 0) then begin
		pc_units, obj=unit, datadir=datadir, /quiet
		units = { velocity:unit.velocity, time:unit.time, temperature:unit.temperature, length:unit.length, density:unit.density, mass:unit.density*unit.length^3, magnetic_field:unit.magnetic_field, default_length:1, default_time:1, default_velocity:1, default_density:1, default_mass:1, default_magnetic_field:1, default_length_str:'m', default_time_str:'s', default_velocity_str:'m/s', default_density_str:'kg/m^3', default_mass_str:'kg', default_magnetic_field_str:'Tesla' }
	end
	units_struct = units

	if (keyword_set (time_series)) then ts = time_series
	if (n_elements (ts) le 0) then pc_read_ts, obj=ts, datadir=datadir, /quiet
	time_series = ts

	if (not keyword_set (param)) then pc_read_param, obj=param, datadir=datadir, /quiet
	if (not keyword_set (run_param)) then pc_read_param, obj=run_param, datadir=datadir, /param2, /quiet
	if (not keyword_set (start_time)) then start_time = min (ts.t) * units.time
	if (not keyword_set (end_time)) then end_time = max (ts.t) * units.time

	time_start = start_time
	time_end = end_time
	run_par = run_param
	start_par = param

	draw_timeseries
end

