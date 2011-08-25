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
pro prepare_varset, num, units, coords, varset, overset, dir, params, run_params

	common varset_common, set, overplot, oversets, unit, coord, varsets, varfiles, datadir, sources, param, run_param

	datadir = dir

	unit = units
	coord = coords
	param = params
	run_param = run_params

	varfiles = { title:"-", time:-1.0, loaded:0, number:-1, precalc_done:0 }
	varfiles = replicate (varfiles, num)

	varsets = replicate (varset, num)
	oversets = replicate (overset, num)
end


; Precalculates a data set and loads data, if necessary
pro precalc, i, number=number, varfile=varfile, datadir=dir, dim=dim, grid=grid, param=par, run_param=run_par, varcontent=varcontent, allprocs=allprocs, show_aver=show_aver, time=time

	common varset_common, set, overplot, oversets, unit, coord, varsets, varfiles, datadir, sources, param, run_param

	; Default settings
	default, show_aver, 0
	default, number, i
	default, dir, pc_get_datadir()
	default, datadir, dir
	if (keyword_set (par)) then param = par
	if (keyword_set (run_par)) then run_param = run_par

	if (varfiles[i].number le 0) then varfiles[i].number = number

	if (varfiles[i].loaded eq 0) then begin
		default, varfile, "var.dat"
		if (n_elements (vars) eq 0) then begin
			print, 'Reading: ', varfile, ' ... please wait!'
			if (i eq 0) then begin
				pc_read_var, varfile=varfile, object=vars, datadir=datadir, dim=dim, grid=grid, param=param, par2=run_param, varcontent=varcontent, allprocs=allprocs, /nostats
			endif else begin
				pc_read_var, varfile=varfile, object=vars, datadir=datadir, dim=dim, grid=grid, param=param, par2=run_param, varcontent=varcontent, allprocs=allprocs, /quiet
			endelse
			sources = tag_names (vars)
			precalc_data, number, vars
			print, 'Ready.'
		end
		varfiles[i].title = varfile
		varfiles[i].loaded = 1
		varfiles[i].precalc_done = 1
		varfiles[i].time = vars.t
		time = vars.t
		vars = 0
	end

	if (show_aver) then draw_averages, number
	if (keyword_set (par)) then par = param
	if (keyword_set (run_par)) then run_par = run_param
end


; Precalculates a data set
pro precalc_data, i, vars

	common varset_common, set, overplot, oversets, unit, coord, varsets, varfiles, datadir, sources, param, run_param

	; First and last physical value, excluding ghost cells
	l1 = coord.l1
	l2 = coord.l2
	m1 = coord.m1
	m2 = coord.m2
	n1 = coord.n1
	n2 = coord.n2

	; Target size of local reduced data block
	tx = coord.nx
	ty = coord.ny
	tz = coord.nz

	tags = tag_names (varsets[i])

	; Compute all desired quantities from available source data
	if (any (strcmp (sources, 'uu', /fold_case))) then begin
		if (any (strcmp (tags, 'u_abs', /fold_case))) then begin
			; Absolute velocity
			varsets[i].u_abs = sqrt (congrid (dot2 (vars.uu[l1:l2,m1:m2,n1:n2,*]), tx, ty, tz, /center, /interp)) * unit.velocity / unit.default_velocity
		end
		if (any (strcmp (tags, 'u_x', /fold_case))) then begin
			; Velocity x-component
			varsets[i].u_x = congrid (reform (vars.uu[l1:l2,m1:m2,n1:n2,0]), tx, ty, tz, /center, /interp) * unit.velocity / unit.default_velocity
		end
		if (any (strcmp (tags, 'u_y', /fold_case))) then begin
			; Velocity y-component
			varsets[i].u_y = congrid (reform (vars.uu[l1:l2,m1:m2,n1:n2,1]), tx, ty, tz, /center, /interp) * unit.velocity / unit.default_velocity
		end
		if (any (strcmp (tags, 'u_z', /fold_case))) then begin
			; Velocity z-component
			varsets[i].u_z = congrid (reform (vars.uu[l1:l2,m1:m2,n1:n2,2]), tx, ty, tz, /center, /interp) * unit.velocity / unit.default_velocity
		end
	end
	if (any (strcmp (tags, 'Temp', /fold_case))) then begin
		; Temperature
		if (any (strcmp (sources, 'lnTT', /fold_case))) then begin
			varsets[i].Temp = exp (congrid (vars.lnTT[l1:l2,m1:m2,n1:n2], tx, ty, tz, /center, /interp)) * unit.temperature
		end else if (any (strcmp (sources, 'TT', /fold_case))) then begin
			varsets[i].Temp = congrid (vars.TT[l1:l2,m1:m2,n1:n2], tx, ty, tz, /center, /interp) * unit.temperature
		end
	end
	if (any (strcmp (tags, 'ln_rho', /fold_case))) then begin
		; Natural logarithmic density
		if (any (strcmp (sources, 'lnrho', /fold_case))) then begin
			varsets[i].ln_rho = alog (exp (congrid (vars.lnrho[l1:l2,m1:m2,n1:n2], tx, ty, tz, /center, /interp)) * unit.density / unit.default_density)
		end else if (any (strcmp (sources, 'rho', /fold_case))) then begin
			varsets[i].ln_rho = alog (congrid (vars.rho[l1:l2,m1:m2,n1:n2], tx, ty, tz, /center, /interp) * unit.density / unit.default_density)
		end
	end else if (any (strcmp (tags, 'log_rho', /fold_case))) then begin
		; Logarithmic density
		if (any (strcmp (sources, 'lnrho', /fold_case))) then begin
			varsets[i].log_rho = alog10 (exp (congrid (vars.lnrho[l1:l2,m1:m2,n1:n2], tx, ty, tz, /center, /interp)) * unit.density / unit.default_density)
		end else if (any (strcmp (sources, 'rho', /fold_case))) then begin
			varsets[i].log_rho = alog10 (congrid (vars.rho[l1:l2,m1:m2,n1:n2], tx, ty, tz, /center, /interp) * unit.density / unit.default_density)
		end
	end else if (any (strcmp (tags, 'rho', /fold_case))) then begin
		; Density
		if (any (strcmp (sources, 'lnrho', /fold_case))) then begin
			varsets[i].rho = exp (congrid (vars.lnrho[l1:l2,m1:m2,n1:n2], tx, ty, tz, /center, /interp)) * unit.density / unit.default_density
		end else if (any (strcmp (sources, 'rho', /fold_case))) then begin
			varsets[i].rho = congrid (vars.rho[l1:l2,m1:m2,n1:n2], tx, ty, tz, /center, /interp) * unit.density / unit.default_density
		end
	end
	if (any (strcmp (tags, 'rho_u_z', /fold_case)) and any (strcmp (sources, 'uu', /fold_case))) then begin
		; Vertical component of impulse density
		if (any (strcmp (sources, 'lnrho', /fold_case))) then begin
			varsets[i].rho_u_z = exp (congrid (vars.lnrho[l1:l2,m1:m2,n1:n2], tx, ty, tz, /center, /interp)) * congrid (vars.uu[l1:l2,m1:m2,n1:n2,2], tx, ty, tz, /center, /interp) * unit.density*unit.velocity / (unit.default_density*unit.default_velocity)
		end else if (any (strcmp (sources, 'rho', /fold_case))) then begin
			varsets[i].rho_u_z = congrid (vars.rho[l1:l2,m1:m2,n1:n2], tx, ty, tz, /center, /interp) * congrid (vars.uu[l1:l2,m1:m2,n1:n2,2], tx, ty, tz, /center, /interp) * unit.density*unit.velocity / (unit.default_density*unit.default_velocity)
		endif
	end
	if (any (strcmp (sources, 'aa', /fold_case))) then begin
		if (any (strcmp (tags, 'Ax', /fold_case))) then begin
			; Magnetic vector potential x-component
			varsets[i].Ax = congrid (reform (vars.aa[l1:l2,m1:m2,n1:n2,0]), tx, ty, tz, /center, /interp) * unit.magnetic_field
		end
		if (any (strcmp (tags, 'Ay', /fold_case))) then begin
			; Magnetic vector potential y-component
			varsets[i].Ay = congrid (reform (vars.aa[l1:l2,m1:m2,n1:n2,1]), tx, ty, tz, /center, /interp) * unit.magnetic_field
		end
		if (any (strcmp (tags, 'Az', /fold_case))) then begin
			; Magnetic vector potential z-component
			varsets[i].Az = congrid (reform (vars.aa[l1:l2,m1:m2,n1:n2,2]), tx, ty, tz, /center, /interp) * unit.magnetic_field
		end
		; Magnetic field
		bb = curl (vars.aa) * unit.magnetic_field
		if (any (strcmp (tags, 'bx', /fold_case))) then begin
			; Magnetic field x-component
			varsets[i].bx = congrid (reform (bb[l1:l2,m1:m2,n1:n2,0]), tx, ty, tz, /center, /interp) / unit.default_magnetic_field
		end
		if (any (strcmp (tags, 'by', /fold_case))) then begin
			; Magnetic field y-component
			varsets[i].by = congrid (reform (bb[l1:l2,m1:m2,n1:n2,1]), tx, ty, tz, /center, /interp) / unit.default_magnetic_field
		end
		if (any (strcmp (tags, 'bz', /fold_case))) then begin
			; Magnetic field z-component
			varsets[i].bz = congrid (reform (bb[l1:l2,m1:m2,n1:n2,2]), tx, ty, tz, /center, /interp) / unit.default_magnetic_field
		end
		if (any (strcmp (tags, 'rho_mag', /fold_case))) then begin
			; Magnetic energy density
			varsets[i].rho_mag = congrid (dot2 (bb[l1:l2,m1:m2,n1:n2,*]), tx, ty, tz, /center, /interp)
		end
		mu0_SI = 4.0 * !Pi * 1.e-7
		if (any (strcmp (tags, 'HR_ohm', /fold_case))) then begin
			; Ohming heating rate [W / m^3] = [kg/m^3] * [m/s]^3 / [m]
			varsets[i].HR_ohm = run_param.eta * param.mu0 * congrid (dot2 ((curlcurl (vars.aa))[l1:l2,m1:m2,n1:n2,*]), tx, ty, tz, /center, /interp) * unit.density * unit.velocity^3 / unit.length
		end
		if (any (strcmp (tags, 'j', /fold_case))) then begin
			; Current density [A / m^2]
			if (any (strcmp (tags, 'HR_ohm', /fold_case))) then begin
				varsets[i].j = sqrt (varsets[i].HR_ohm / (run_param.eta * mu0_SI * unit.velocity * unit.length))
			end else begin
				varsets[i].j = sqrt (congrid (dot2 ((curlcurl (vars.aa))[l1:l2,m1:m2,n1:n2,*]), tx, ty, tz, /center, /interp)) * unit.velocity * sqrt (param.mu0 / mu0_SI * unit.density) / unit.length
			end
		end
	end

	over_tags = tag_names (oversets[i])
	if (any (strcmp (sources, 'uu', /fold_case))) then begin
		if (any (strcmp (over_tags, 'u', /fold_case))) then begin
			; Velocity overplot
			oversets[i].u[*,*,*,0] = float (congrid (reform (vars.uu[l1:l2,m1:m2,n1:n2,0]), tx, ty, tz, /center, /interp) * unit.velocity / unit.default_velocity)
			oversets[i].u[*,*,*,1] = float (congrid (reform (vars.uu[l1:l2,m1:m2,n1:n2,1]), tx, ty, tz, /center, /interp) * unit.velocity / unit.default_velocity)
			oversets[i].u[*,*,*,2] = float (congrid (reform (vars.uu[l1:l2,m1:m2,n1:n2,2]), tx, ty, tz, /center, /interp) * unit.velocity / unit.default_velocity)
		end
	end
	if (any (strcmp (sources, 'aa', /fold_case))) then begin
		if (any (strcmp (over_tags, 'b', /fold_case))) then begin
			; Magnetic field overplot
			oversets[i].b[*,*,*,0] = float (congrid (reform (bb[l1:l2,m1:m2,n1:n2,0]), tx, ty, tz, /center, /interp) / unit.default_magnetic_field)
			oversets[i].b[*,*,*,1] = float (congrid (reform (bb[l1:l2,m1:m2,n1:n2,1]), tx, ty, tz, /center, /interp) / unit.default_magnetic_field)
			oversets[i].b[*,*,*,2] = float (congrid (reform (bb[l1:l2,m1:m2,n1:n2,2]), tx, ty, tz, /center, /interp) / unit.default_magnetic_field)
		end
		bb = 0
		if (any (strcmp (over_tags, 'a_contour', /fold_case))) then begin
			; Magnetic field lines overplot
			oversets[i].a_contour[*,*,*,0] = float (congrid (reform (vars.aa[l1:l2,m1:m2,n1:n2,0]), tx, ty, tz, /center, /interp) * unit.magnetic_field)
			oversets[i].a_contour[*,*,*,1] = float (congrid (reform (vars.aa[l1:l2,m1:m2,n1:n2,1]), tx, ty, tz, /center, /interp) * unit.magnetic_field)
			oversets[i].a_contour[*,*,*,2] = float (congrid (reform (vars.aa[l1:l2,m1:m2,n1:n2,2]), tx, ty, tz, /center, /interp) * unit.magnetic_field)
		end
	end
end


; Show timeseries analysis window
pro show_timeseries, ts, tags, unit, param, run_param, start_time=start_time, end_time=end_time

	if (n_elements (ts) gt 0) then begin

		default, start_time, 0
		default, end_time, 0
		add_title = ''
		if (start_time > 0) then add_title = ' (starting at the frist selected snapshot)'
		if (end_time > 0) then add_title = ' (ending at the last selected snapshot)'
		if ((start_time > 0) and (end_time > 0)) then add_title = ' (showing only selected snapshots)'

		window, 11, xsize=1000, ysize=400, title='timestep analysis'+add_title, retain=2
		!P.MULTI = [0, 2, 1]

		print, "starting values:"
		print, "dt    :", ts.dt[0]
		plot, ts.dt, title = 'dt', /yl

		tags = tag_names (ts)
		if (any (strcmp (tags, 't', /fold_case))) then begin
			time = ts.t
		endif else begin
			time = ts.it
		endelse
		x_minmax = minmax (time > start_time)
		if (end_time > 0) then x_minmax = minmax (x_minmax < end_time)
		y_minmax = minmax (ts.dt)
		if (any (strcmp (tags, 'dtu', /fold_case)))    then y_minmax = minmax ([y_minmax, ts.dtu])
		if (any (strcmp (tags, 'dtv', /fold_case)))    then y_minmax = minmax ([y_minmax, ts.dtv])
		if (any (strcmp (tags, 'dtnu', /fold_case)))   then y_minmax = minmax ([y_minmax, ts.dtnu])
		if (any (strcmp (tags, 'dtb', /fold_case)))    then y_minmax = minmax ([y_minmax, ts.dtb])
		if (any (strcmp (tags, 'dteta', /fold_case)))  then y_minmax = minmax ([y_minmax, ts.dteta])
		if (any (strcmp (tags, 'dtc', /fold_case)))    then y_minmax = minmax ([y_minmax, ts.dtc])
		if (any (strcmp (tags, 'dtchi', /fold_case)))  then y_minmax = minmax ([y_minmax, ts.dtchi])
		if (any (strcmp (tags, 'dtchi2', /fold_case))) then y_minmax = minmax ([y_minmax, ts.dtchi2])
		if (any (strcmp (tags, 'dtd', /fold_case)))    then y_minmax = minmax ([y_minmax, ts.dtd])

		time *= unit.time
		ts.dt *= unit.time
		x_minmax *= unit.time
		y_minmax *= unit.time

		plot, time, ts.dt, title = 'dt(t) u{-t} v{-p} nu{.v} b{.r} eta{-g} c{.y} chi{-.b} chi2{-.o} d{-l} [s]', xrange=x_minmax, /xs, yrange=y_minmax, /yl
		if (any (strcmp (tags, 'dtu', /fold_case))) then begin
			oplot, time, ts.dtu*unit.time, linestyle=2, color=11061000
			print, "dtu   :", ts.dtu[0]
		end
		if (any (strcmp (tags, 'dtv', /fold_case))) then begin
			oplot, time, ts.dtv*unit.time, linestyle=2, color=128255200
			print, "dtv   :", ts.dtv[0]
		end
		if (any (strcmp (tags, 'dtnu', /fold_case))) then begin
			oplot, time, ts.dtnu*unit.time, linestyle=1, color=128000128
			print, "dtnu  :", ts.dtnu[0]
		end
		if (any (strcmp (tags, 'dtb', /fold_case))) then begin
			oplot, time, ts.dtb*unit.time, linestyle=1, color=200
			print, "dtb   :", ts.dtb[0]
		end
		if (any (strcmp (tags, 'dteta', /fold_case))) then begin
			oplot, time, ts.dteta*unit.time, linestyle=2, color=220200200
			print, "dteta :", ts.dteta[0]
		end
		if (any (strcmp (tags, 'dtc', /fold_case))) then begin
			oplot, time, ts.dtc*unit.time, linestyle=1, color=61695
			print, "dtc   :", ts.dtc[0]
		end
		if (any (strcmp (tags, 'dtchi', /fold_case))) then begin
			oplot, time, ts.dtchi*unit.time, linestyle=3, color=115100200
			print, "dtchi :", ts.dtchi[0]
		end
		if (any (strcmp (tags, 'dtchi2', /fold_case))) then begin
			oplot, time, ts.dtchi2*unit.time, linestyle=3, color=41215
			print, "dtchi2:", ts.dtchi2[0]
		end
		if (any (strcmp (tags, 'dtd', /fold_case))) then begin
			oplot, time, ts.dtd*unit.time, linestyle=2, color=16737000
			print, "dtc   :", ts.dtd[0]
		end

		window, 12, xsize=1000, ysize=800, title='time series analysis'+add_title, retain=2
		!P.MULTI = [0, 2, 2]

		max_subplots = 4
		num_subplots = 0

		if (any (strcmp (tags, 'eem', /fold_case)) and any (strcmp (tags, 'ethm', /fold_case)) and any (strcmp (tags, 'ekintot', /fold_case)) and any (strcmp (tags, 'totmass', /fold_case)) and (num_subplots lt max_subplots)) then begin
			num_subplots += 1
			mass = ts.totmass * unit.mass / unit.default_mass
			energy = (ts.eem + ts.ekintot/ts.totmass) * unit.mass / unit.velocity^2
			plot, time, energy, linestyle=2, title = 'Mass {.r} and energy {-w} conservation', xrange=x_minmax, /xs, ys=8
			oplot, time, mass*mean (energy)/mean (mass), linestyle=1, color=200
			axis, yaxis=0, yrange=!Y.CRANGE, /ys, ytitle='mean energy [J]'
			axis, yaxis=1, yrange=!Y.CRANGE*mean (energy)/mean (mass), /ys, ytitle='total mass ['+unit.default_mass_str+']'
		end else if (any (strcmp (tags, 'totmass', /fold_case)) and (num_subplots lt max_subplots)) then begin
			num_subplots += 1
			mass = ts.totmass * unit.mass / unit.default_mass
			plot, time, mass, title = 'Mass conservation', xrange=x_minmax, /xs, ys=8
		end
		if (any (strcmp (tags, 'TTmax', /fold_case)) and any (strcmp (tags, 'j2m', /fold_case)) and (num_subplots lt max_subplots)) then begin
			num_subplots += 1
			Temp_max = ts.TTmax * unit.temperature
			HR_ohm = run_param.eta * param.mu0 * ts.j2m * unit.density * unit.velocity^3 / unit.length
			plot, time, Temp_max, title = 'Maximum temperature [K] {-w} and average Ohmic heating rate [W/m^3] {.r}', xrange=x_minmax, /xs, /yl, ys=8
			oplot, time, HR_ohm*mean (Temp_max)/mean (HR_ohm), linestyle=1, color=200
			axis, yaxis=0, yrange=!Y.CRANGE, /ys, ytitle='Maximum temperature [K]'
			axis, yaxis=1, yrange=!Y.CRANGE*mean (Temp_max)/mean (HR_ohm), /ys, ytitle='HR_ohm [W/m^3]'
		end else if (any (strcmp (tags, 'TTmax', /fold_case)) and (num_subplots lt max_subplots)) then begin
			num_subplots += 1
			Temp_max = ts.TTmax * unit.temperature
			plot, time, Temp_max, title = 'Maximum temperature [K]', xrange=x_minmax, /xs, /yl
		end else if (any (strcmp (tags, 'j2m', /fold_case)) and (num_subplots lt max_subplots)) then begin
			num_subplots += 1
			HR_ohm = run_param.eta * param.mu0 * ts.j2m * unit.density * unit.velocity^3 / unit.length
			j_abs = sqrt (ts.j2m) * unit.velocity * sqrt (param.mu0 / mu0_SI * unit.density) / unit.length
			plot, time, HR_ohm, title = 'Ohmic heating rate [W/m^2] {-w} and mean current density [A/m^2] {.r}', xrange=x_minmax, /xs, /yl, ys=8
			oplot, time, j_abs*mean (HR_ohm)/mean (j_abs), linestyle=1, color=200
			axis, yaxis=0, yrange=!Y.CRANGE, /ys, ytitle='Ohmic heating rate [W/m^3]'
			axis, yaxis=1, yrange=!Y.CRANGE*mean (HR_ohm)/mean (j_abs), /ys, ytitle='current density [A/m^2]'
		end
		if (any (strcmp (tags, 'umax', /fold_case)) and (num_subplots lt max_subplots)) then begin
			num_subplots += 1
			u_max = ts.umax * unit.velocity / unit.default_velocity
			u_title = 'u_max(t){-w}'
			if (any (strcmp (tags, 'urms', /fold_case))) then begin
				u_title += ' u_rms{.r}'
			end else if (any (strcmp (tags, 'u2m', /fold_case))) then u_title += ' <u^2>^0.5{.-b}'
			plot, time, u_max, title = u_title+' ['+unit.default_velocity_str+']', xrange=x_minmax, /xs
			if (any (strcmp (tags, 'urms', /fold_case))) then begin
				urms = ts.urms * unit.velocity / unit.default_velocity
				oplot, time, urms, linestyle=1, color=200
			end else if (any (strcmp (tags, 'u2m', /fold_case))) then begin
				u2m = sqrt (ts.u2m) * unit.velocity / unit.default_velocity
				oplot, time, u2m, linestyle=3, color=115100200
			end
		end
		if (any (strcmp (tags, 'rhomin', /fold_case)) and (num_subplots lt max_subplots)) then begin
			num_subplots += 1
			rho_min = ts.rhomin * unit.density / unit.default_density
			plot, time, rho_min, title = 'rho_min(t) ['+unit.default_density_str+']', xrange=x_minmax, /xs, /yl
		end
	end
end


; Dummy routine
pro pc_gui_companion

	pc_gui_companion_loaded = 1
end

