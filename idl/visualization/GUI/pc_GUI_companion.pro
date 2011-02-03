;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_GUI_companion.pro     ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  $Id$
;;;
;;;  Description:
;;;   Framework for precalculation and comparision of output in pencil units.
;;;   Companion procedures needed by 'pc_GUI.pro'.
;;;
;;;  To do:
;;;   Add more comments


; Prepares the varset
pro prepare_varset, num, units, coords, varset, overset

	common varset_common, set, overplot, oversets, unit, coord, varsets, varfiles, sources

	unit = units
	coord = { x:coords.x*unit.length, y:coords.y*unit.length, z:coords.z*unit.length }

	varfiles = { title:"-", loaded:0, number:0, precalc_done:0 }
	varfiles = replicate (varfiles, num)

	varsets = replicate (varset, num)
	oversets = replicate (overset, num)
end


; Precalculates a data set and loads data, if necessary
pro precalc, i, number=number, varfile=varfile, show_aver=show_aver, vars=vars

	common varset_common, set, overplot, oversets, unit, coord, varsets, varfiles, sources

	; Default settings
	default, show_aver, 0
	default, number, i

	if (varfiles[i].number le 0) then varfiles[i].number = number

	if (varfiles[i].loaded eq 0) then begin
		default, varfile, "var.dat"
		if (n_elements (vars) eq 0) then begin
			print, 'Reading: ', varfile, ' ... please wait!'
			pc_read_var, varfile=varfile, object=vars, /quiet
			sources = tag_names (vars)
		end
		varfiles[i].title = varfile
		varfiles[i].loaded = 1
		varfiles[i].precalc_done = 0
	end

	if (varfiles[i].precalc_done eq 0) then begin
		print, 'Calculation started...'
		precalc_data, number, vars
		varfiles[i].precalc_done = 1
		print, 'Ready.'
	end

	if (show_aver) then draw_averages, number
end


; Precalculates a data set
pro precalc_data, i, vars

	common varset_common, set, overplot, oversets, unit, coord, varsets, varfiles, sources

	tags = tag_names (varsets[i])

	; Compute all desired quantities from available source data
	if (any (strcmp (sources, 'uu', /fold_case))) then begin
		if (any (strcmp (tags, 'u_abs', /fold_case))) then begin
			; Absolute velocity
			varsets[i].u_abs = sqrt (dot2 (vars.uu)) * unit.velocity / unit.default_velocity
		end
		if (any (strcmp (tags, 'u_x', /fold_case))) then begin
			; Velocity x-component
			varsets[i].u_x = vars.uu[*,*,*,0] * unit.velocity / unit.default_velocity
		end
		if (any (strcmp (tags, 'u_y', /fold_case))) then begin
			; Velocity y-component
			varsets[i].u_y = vars.uu[*,*,*,1] * unit.velocity / unit.default_velocity
		end
		if (any (strcmp (tags, 'u_z', /fold_case))) then begin
			; Velocity z-component
			varsets[i].u_z = vars.uu[*,*,*,2] * unit.velocity / unit.default_velocity
		end
	end
	if (any (strcmp (tags, 'Temp', /fold_case))) then begin
		; Temperature
		if (any (strcmp (sources, 'lnTT', /fold_case))) then begin
			varsets[i].Temp = exp (vars.lnTT) * unit.temperature
		end else if (any (strcmp (sources, 'TT', /fold_case))) then begin
			varsets[i].Temp = vars.TT * unit.temperature
		end
	end
	if (any (strcmp (sources, 'aa', /fold_case))) then begin
		if (any (strcmp (tags, 'Ax', /fold_case))) then begin
			; Magnetic vector potential x-component
			varsets[i].Ax = vars.aa[*,*,*,0]
		end
		if (any (strcmp (tags, 'Ay', /fold_case))) then begin
			; Magnetic vector potential y-component
			varsets[i].Ay = vars.aa[*,*,*,1]
		end
		if (any (strcmp (tags, 'Az', /fold_case))) then begin
			; Magnetic vector potential z-component
			varsets[i].Az = vars.aa[*,*,*,2]
		end
		; Magnetic field
		bb = curl (vars.aa) / unit.length
		if (any (strcmp (tags, 'bx', /fold_case))) then begin
			; Magnetic field x-component
			varsets[i].bx = bb[*,*,*,0]
		end
		if (any (strcmp (tags, 'by', /fold_case))) then begin
			; Magnetic field y-component
			varsets[i].by = bb[*,*,*,1]
		end
		if (any (strcmp (tags, 'bz', /fold_case))) then begin
			; Magnetic field z-component
			varsets[i].bz = bb[*,*,*,2]
		end
		if (any (strcmp (tags, 'rho_mag', /fold_case))) then begin
			; Magnetic energy density
			varsets[i].rho_mag = dot2 (bb)
		end
		if (any (strcmp (tags, 'j', /fold_case))) then begin
			; Current density
			varsets[i].j = sqrt (sqrt (dot2 (curlcurl (vars.aa))) / unit.length^2)
		end
	end
	if (any (strcmp (tags, 'ln_rho', /fold_case))) then begin
		; Natural logarithmic density
		if (any (strcmp (sources, 'lnrho', /fold_case))) then begin
			varsets[i].ln_rho = alog (exp (vars.lnrho) * unit.density / unit.default_density)
		end else if (any (strcmp (sources, 'rho', /fold_case))) then begin
			varsets[i].ln_rho = alog (vars.rho * unit.density / unit.default_density)
		end
		if (any (strcmp (tags, 'rho_u_z', /fold_case)) and any (strcmp (sources, 'uu', /fold_case))) then begin
			; Vertical component of impulse density
			varsets[i].rho_u_z = exp (varsets[i].ln_rho) * vars.uu[*,*,*,2] * unit.velocity / unit.default_velocity
		end
	end else if (any (strcmp (tags, 'log_rho', /fold_case))) then begin
		; Logarithmic density
		if (any (strcmp (sources, 'lnrho', /fold_case))) then begin
			varsets[i].log_rho = alog10 (exp (vars.lnrho) * unit.density / unit.default_density)
		end else if (any (strcmp (sources, 'rho', /fold_case))) then begin
			varsets[i].log_rho = alog10 (vars.rho * unit.density / unit.default_density)
		end
		if (any (strcmp (tags, 'rho_u_z', /fold_case)) and any (strcmp (sources, 'uu', /fold_case))) then begin
			; Vertical component of impulse density
			varsets[i].rho_u_z = 10.0^(varsets[i].log_rho) * vars.uu[*,*,*,2] * unit.velocity / unit.default_velocity
		end
	end else if (any (strcmp (tags, 'rho', /fold_case))) then begin
		; Density
		if (any (strcmp (sources, 'lnrho', /fold_case))) then begin
			varsets[i].rho = exp (vars.lnrho) * unit.density / unit.default_density
		end else if (any (strcmp (sources, 'rho', /fold_case))) then begin
			varsets[i].rho = vars.rho * unit.density / unit.default_density
		end
		if (any (strcmp (tags, 'rho_u_z', /fold_case)) and any (strcmp (sources, 'uu', /fold_case))) then begin
			; Vertical component of impulse density
			varsets[i].rho_u_z = varsets[i].rho * vars.uu[*,*,*,2] * unit.velocity / unit.default_velocity
		end
	end

	over_tags = tag_names (oversets[i])
	if (any (strcmp (sources, 'uu', /fold_case))) then begin
		if (any (strcmp (over_tags, 'u', /fold_case))) then begin
			; Velocity overplot
			oversets[i].u = float (vars.uu * unit.velocity / unit.default_velocity)
		end
	end
	if (any (strcmp (sources, 'aa', /fold_case))) then begin
		if (any (strcmp (over_tags, 'b', /fold_case))) then begin
			; Magnetic field overplot
			oversets[i].b = float (bb)
		end
		if (any (strcmp (over_tags, 'a_contour', /fold_case))) then begin
			; Magnetic field lines overplot
			oversets[i].a_contour = float (vars.aa)
		end
	end
end


; Show timeseries analysis window
pro show_timeseries, ts, tags, unit, start_time=start_time

	if (n_elements (ts) gt 0) then begin

		default, start_time, 0
		add_title = ''
		if (start_time > 0) then add_title = ' (starting at the frist selected snapshot)'

		window, 11, xsize=1000, ysize=400, title='timestep analysis'+add_title, retain=2
		!P.MULTI = [0, 2, 1]

		print, "starting values:"
		print, "dt    :", ts.dt[0]
		plot, ts.dt, title = 'dt', /yl

		tags = tag_names (ts)
		x_minmax = minmax (ts.t > start_time)
		y_minmax = minmax (ts.dt)
		if (any (strcmp (tags, 'dtu', /fold_case)))    then y_minmax = minmax ([y_minmax, ts.dtu])
		if (any (strcmp (tags, 'dtv', /fold_case)))    then y_minmax = minmax ([y_minmax, ts.dtv])
		if (any (strcmp (tags, 'dtnu', /fold_case)))   then y_minmax = minmax ([y_minmax, ts.dtnu])
		if (any (strcmp (tags, 'dtb', /fold_case)))    then y_minmax = minmax ([y_minmax, ts.dtb])
		if (any (strcmp (tags, 'dteta', /fold_case)))  then y_minmax = minmax ([y_minmax, ts.dteta])
		if (any (strcmp (tags, 'dtc', /fold_case)))    then y_minmax = minmax ([y_minmax, ts.dtc])
		if (any (strcmp (tags, 'dtchi', /fold_case)))  then y_minmax = minmax ([y_minmax, ts.dtchi])
		if (any (strcmp (tags, 'dtchi2', /fold_case))) then y_minmax = minmax ([y_minmax, ts.dtchi2])

		ts.t *= unit.time
		ts.dt *= unit.time
		x_minmax *= unit.time
		y_minmax *= unit.time

		plot, ts.t, ts.dt, title = 'dt(tt) u{-t} v{-p} nu{.v} b{.r} eta{-g} c{.y} chi{-.b} chi2{-.o} [s]', xrange=x_minmax, /xs, yrange=y_minmax, /yl
		if (any (strcmp (tags, 'dtu', /fold_case))) then begin
			oplot, ts.t, ts.dtu*unit.time, linestyle=2, color=11061000
			print, "dtu   :", ts.dtu[0]
		end
		if (any (strcmp (tags, 'dtv', /fold_case))) then begin
			oplot, ts.t, ts.dtv*unit.time, linestyle=2, color=128255200
			print, "dtv   :", ts.dtv[0]
		end
		if (any (strcmp (tags, 'dtnu', /fold_case))) then begin
			oplot, ts.t, ts.dtnu*unit.time, linestyle=1, color=128000128
			print, "dtnu  :", ts.dtnu[0]
		end
		if (any (strcmp (tags, 'dtb', /fold_case))) then begin
			oplot, ts.t, ts.dtb*unit.time, linestyle=1, color=200
			print, "dtb   :", ts.dtb[0]
		end
		if (any (strcmp (tags, 'dteta', /fold_case))) then begin
			oplot, ts.t, ts.dteta*unit.time, linestyle=2, color=220200200
			print, "dteta :", ts.dteta[0]
		end
		if (any (strcmp (tags, 'dtc', /fold_case))) then begin
			oplot, ts.t, ts.dtc*unit.time, linestyle=1, color=61695
			print, "dtc   :", ts.dtc[0]
		end
		if (any (strcmp (tags, 'dtchi', /fold_case))) then begin
			oplot, ts.t, ts.dtchi*unit.time, linestyle=3, color=115100200
			print, "dtchi :", ts.dtchi[0]
		end
		if (any (strcmp (tags, 'dtchi2', /fold_case))) then begin
			oplot, ts.t, ts.dtchi2*unit.time, linestyle=3, color=41215
			print, "dtchi2:", ts.dtchi2[0]
		end

		window, 12, xsize=1000, ysize=800, title='time series analysis'+add_title, retain=2
		!P.MULTI = [0, 2, 2]

		max_subplots = 4
		num_subplots = 0

		if (any (strcmp (tags, 'eem', /fold_case)) and any (strcmp (tags, 'ethm', /fold_case)) and any (strcmp (tags, 'ekintot', /fold_case)) and any (strcmp (tags, 'totmass', /fold_case)) and (num_subplots lt max_subplots)) then begin
			num_subplots += 1
			mass = ts.totmass * unit.mass / unit.default_mass
			energy = (ts.eem + ts.ekintot/ts.totmass) * unit.mass / unit.velocity^2
			plot, ts.t, energy, title = 'Mass and energy conservation', xrange=x_minmax, xtitle='mean energy [J]', /xs
			oplot, ts.t, mass*mean (energy)/mean (mass)
			axis, yaxis=1, yrange=!Y.CRANGE*mean (energy)/mean (mass), /ys, ytitle='total mass ['+unit.default_mass_str+']'
		end
		if (any (strcmp (tags, 'TTmax', /fold_case)) and (num_subplots lt max_subplots)) then begin
			num_subplots += 1
			Temp_max = ts.TTmax * unit.temperature
			plot, ts.t, Temp_max, title = 'Temp_max(tt) [K]', xrange=x_minmax, /xs, /yl
		end
		if (any (strcmp (tags, 'umax', /fold_case)) and (num_subplots lt max_subplots)) then begin
			num_subplots += 1
			u_max = ts.umax * unit.velocity / unit.default_velocity
			plot, ts.t, u_max, title = 'u_max(tt) ['+unit.default_velocity_str+']', xrange=x_minmax, /xs
		end
		if (any (strcmp (tags, 'rhomin', /fold_case)) and (num_subplots lt max_subplots)) then begin
			num_subplots += 1
			rho_min = ts.rhomin * unit.density / unit.default_density
			plot, ts.t, rho_min, title = 'rho_min(tt) ['+unit.default_density_str+']', xrange=x_minmax, /xs, /yl
		end
	end
end


; Dummy routine
pro pc_GUI_companion

	pc_GUI_companion_loaded = 1
end

